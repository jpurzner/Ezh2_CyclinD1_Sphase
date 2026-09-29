# Spinning up AWS compute for the ChIP × bursting work

Written 2026-09-29. Assumes you already have an AWS account with data in S3. Nothing here has
been run — I have no credentials on this machine, and launching billable infrastructure is yours
to do.

---

## 0. First: you probably need less than you think

Sizing the two workloads separately, because they are very different jobs:

| workload | what it does | honest requirement |
|---|---|---|
| **Burst imputation** (public ChIP × Larsson kinetics) | download processed BED/bigWig, aggregate to gene level, fit a regression | **8 vCPU / 16–32 GB**. Embarrassingly parallel, IO-bound. The matrix it produces is a few MB. |
| **Your own GNP/MB ChIP** from FASTQ | adapter trim → align → dedup → peak call → signal | **16–32 vCPU, 64 GB, fast scratch**. This is the real compute. |

If the plan is only the first one, a `c6i.2xlarge` (~$0.34/hr) is ample and a big instance is
wasted money. If you are also aligning your own reads, that justifies the larger box — but then
size it for the alignment, not the analysis.

**Split them.** They have different shapes and different lifetimes: alignment is a few hours of
heavy batch, the analysis is interactive and iterative.

## 1. Decide region first — it is the one hard-to-undo choice

**Put compute where your data already is.** S3→EC2 within a region is free and fast; cross-region
is charged and slow.

ENCODE does not pull you anywhere: `s3://encode-public` is in **us-west-2**, but it is part of the
**AWS Open Data Sponsorship Program**, so the bucket owner absorbs egress and reading it from
another region costs you nothing. You pay only in latency. *(Worth confirming against current AWS
billing docs before a large transfer — sponsorship terms do change.)*

So: if your lab data is in `ca-central-1`, stay there. Do not move to us-west-2 for ENCODE's sake.

## 2. Instance choice

| use | type | vCPU / RAM | ~on-demand | notes |
|---|---|---|---|---|
| analysis only | `c6i.2xlarge` | 8 / 16 GB | ~$0.34/hr | plenty for Phase 1 + 2 |
| analysis, roomier | `m6i.2xlarge` | 8 / 32 GB | ~$0.38/hr | if you hold matrices in memory |
| **+ your own alignment** | `c6id.8xlarge` | 32 / 64 GB | ~$1.62/hr | the `d` gives **local NVMe scratch**, which matters a lot for BAM IO |
| cheapest for batch | any of the above **on Spot** | — | **60–90% off** | fine for interruptible work; not for a long interactive session |

`c6id`/`m6id` (with `d`) include instance-store NVMe. For alignment that is worth more than extra
cores — sorting BAMs on EBS is often the bottleneck.

Rough cost: a full 8-hour day on `c6i.2xlarge` is **under $3**. Even the 32-core box for a day is
~$13 on demand, ~$4 on Spot. This is not an expensive problem.

## 3. Storage

* **Root EBS `gp3`, 100–200 GB.** gp3 gives 3000 IOPS / 125 MB/s baseline at no extra charge —
  do not use `gp2`. ~$0.08/GB-month.
* **Scratch**: instance-store NVMe if you chose a `d` type (free, ephemeral, fast). Otherwise add
  a second gp3 volume.
* **Keep durable outputs in S3**, not on EBS. The instance is cattle; the bucket is the record.
* **Do not** put the 28 GB of bigWigs on an expensive volume you keep. Pull, aggregate, push the
  few-MB matrix to S3, delete the instance.

## 4. Access — use a role, not keys

Do **not** copy access keys onto the instance.

1. Create an IAM role, e.g. `chip-analysis-ec2`.
2. Attach a policy granting read on ENCODE (public, no auth needed) and read/write on **only your
   project prefix**:

```json
{ "Version": "2012-10-17", "Statement": [
  { "Effect": "Allow",
    "Action": ["s3:GetObject","s3:PutObject","s3:ListBucket"],
    "Resource": ["arn:aws:s3:::YOUR-BUCKET","arn:aws:s3:::YOUR-BUCKET/chip-bursting/*"] }
]}
```

3. Attach the role at launch (*Advanced details → IAM instance profile*).
4. Add an **S3 gateway VPC endpoint** — it is free, keeps S3 traffic off the NAT gateway, and NAT
   egress charges are the usual way a small analysis instance produces a surprising bill.

## 5. Launch

Console → EC2 → Launch instance:

* **AMI**: Ubuntu 24.04 LTS (x86_64). Avoid ARM/Graviton unless you want to check every
  bioinformatics binary — the savings are not worth the debugging here.
* **Type**: from §2.
* **Key pair**: existing, or use **EC2 Instance Connect / SSM Session Manager** and skip SSH keys.
* **Network**: your existing VPC; private subnet + SSM is cleaner than a public IP.
* **Storage**: 150 GB gp3.
* **IAM instance profile**: the role from §4.
* **Tags**: `Project=ezh2-bursting`, `Owner=<you>` — makes the bill legible later.

Then:

```bash
sudo apt-get update && sudo apt-get install -y python3-pip bedtools awscli
pip3 install --break-system-packages pyBigWig numpy pandas scikit-learn pyarrow
```

## 6. Set a budget alarm before you start

Billing → Budgets → monthly cost budget, e.g. $50, alert at 80%. Two minutes, and it is the
difference between noticing a forgotten instance in a day versus at month end.

The classic waste here is **a stopped-but-not-terminated instance still holding a large EBS
volume**. Stopping does not stop EBS charges. Terminate when done.

## 7. What to actually run, in order

1. **Phase 1** (`~400 MB`, minutes): pull ChIP-Atlas peak BEDs for the chosen antigens, filter to
   **wild-type experiments only** (`protein_noise_bursting_chromatin.md` §2c — 37% of the mESC
   Suz12 set is KO-flagged and it changed the answer), build gene × feature matrices.
2. **Train mESC → test GM12878.** This is the decision point and needs no lab data. If held-out
   cross-cell-type R² is near zero, stop; the mapping does not transfer and your own ChIP will not
   rescue it.
3. **Phase 2** only if step 2 justifies it: continuous RPM from bigWigs, matching Ochiai's
   measurement (promoter −1000/+100 and gene body separately).
4. Push the derived matrices to S3, terminate, continue modelling locally — they are a few MB.

## 8. If you are also processing your own ChIP

Separate job, separate instance, and the sequence is: reference + index to instance-store scratch
→ align → sort/dedup → peak call and signal → push BAM/bigWig/peaks to S3 → terminate. Use Spot
with checkpointing, or just accept on-demand for a few hours since the total is small.

Worth considering **AWS Batch** if this becomes routine — it queues per-sample jobs and terminates
instances automatically, which removes the single biggest cost risk (leaving something running).

---

## Open questions that change the above

1. **Which region is your lab data in?** Compute goes there.
2. **Is your own FASTQ→BAM alignment in scope**, or only the public analysis? This changes the
   instance by ~4× and is the difference between "small box" and "real box".
3. **Existing VPC/IAM, or starting fresh?** Affects whether §4–5 are a 10-minute or 45-minute job.
