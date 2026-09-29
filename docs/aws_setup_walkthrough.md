# AWS setup: ChIP × bursting analysis + GNP/MB alignment

**Your setup:** account exists, lab data in S3, **region `us-east-1`**, and **both** workloads in
scope — your own GNP/MB ChIP from FASTQ *and* the public-data burst-imputation analysis.

Nothing here has been run. I have no credentials on this machine, and launching billable
infrastructure is yours to do.

---

## 0. Use two instances, not one

They are different jobs with different lifetimes, and conflating them is how you end up paying
for 32 idle cores while you stare at a regression.

| | workload | instance | us-east-1 on-demand | lifetime |
|---|---|---|---|---|
| **A** | your ChIP: FASTQ → BAM → peaks/signal | `c6id.8xlarge` (32 vCPU, 64 GB, **1.9 TB NVMe**) | ~$1.61/hr | hours, then terminate |
| **B** | burst imputation + modelling | `c6i.2xlarge` (8 vCPU, 16 GB) | ~$0.34/hr | days, intermittent, stop when idle |

**Whole thing should land under $50.** Alignment of ~20 ChIP samples is ~2–3 h wall-clock
(~$4 on-demand, ~$1.50 on Spot); the analysis box at a few hours a day for a week is ~$10–15.

`us-east-1` is the cheapest region and has the widest instance availability, so you are already in
the right place. ENCODE's `s3://encode-public` lives in **us-west-2**, so those reads are
cross-region — but it is an **AWS Open Data** bucket, so the sponsor absorbs egress and it costs
you nothing beyond latency. *(Re-check against current billing docs before a very large pull;
sponsorship terms do change.)* Do **not** move compute to us-west-2 for ENCODE's sake — stay with
your own data.

---

## 1. Instance A — alignment

### Why `c6id` specifically

The `d` is the point: **1.9 TB of local NVMe**. BAM sorting and dedup are IO-bound, and doing them
on EBS is usually the bottleneck — you would need expensive provisioned IOPS to match what the
instance store gives you free.

32 vCPU is sized for **4 concurrent samples at 8 threads each**, not for one 32-thread job —
bowtie2/BWA scale poorly past ~16 threads, so parallelism across samples wins.

### Space budget

| | |
|---|---|
| mm10 + bowtie2 index | ~5 GB |
| FASTQ, 20 samples PE @ ~40M reads | ~80 GB |
| BAM + intermediates | ~150 GB |
| **total** | **~250 GB** → comfortable on 1.9 TB |

### The gotcha: instance store is not mounted for you

It is also **ephemeral** — wiped on stop *or* terminate. Anything you want to keep goes to S3.

```bash
lsblk                                    # find the NVMe, usually nvme1n1
sudo mkfs.xfs /dev/nvme1n1
sudo mkdir -p /scratch && sudo mount /dev/nvme1n1 /scratch
sudo chown $USER:$USER /scratch
df -h /scratch
```

### Setup

```bash
sudo apt-get update
sudo apt-get install -y bowtie2 samtools bedtools macs2 fastp pigz awscli
# or conda/mamba if you want pinned versions
```

### Shape of the run

```bash
cd /scratch
aws s3 cp s3://YOUR-BUCKET/refs/mm10_bowtie2.tar.gz . && tar xzf mm10_bowtie2.tar.gz
aws s3 sync s3://YOUR-BUCKET/gnp-chip/fastq/ ./fastq/

ls fastq/*_R1.fastq.gz | sed 's/_R1.*//' | xargs -P 4 -I{} bash -c '
  s=$(basename {});
  fastp -i {}_R1.fastq.gz -I {}_R2.fastq.gz -o {}_1.fq.gz -O {}_2.fq.gz -w 4 -j $s.json -h $s.html
  bowtie2 -p 8 -x mm10 -1 {}_1.fq.gz -2 {}_2.fq.gz \
    | samtools sort -@ 4 -m 2G -T /scratch/tmp_$s -o $s.bam
  samtools index $s.bam
'
aws s3 sync . s3://YOUR-BUCKET/gnp-chip/bam/ --exclude "*" --include "*.bam*"
```

`-P 4` × `-p 8` = 32 threads. Adjust if your sample count differs.

**Spot is a good fit here** (60–90% off) as long as you checkpoint per sample — each finished BAM
goes to S3 immediately, so an interruption costs you one sample, not the run.

**Terminate when done.** Stopping still bills the EBS root volume, and the NVMe contents are gone
either way.

---

## 2. Instance B — the analysis box

`c6i.2xlarge`, 150 GB gp3 root (**gp3, not gp2** — 3000 IOPS / 125 MB/s baseline included).

```bash
sudo apt-get update && sudo apt-get install -y python3-pip bedtools awscli
pip3 install --break-system-packages pyBigWig numpy pandas scikit-learn pyarrow
```

For the bulk ChIP-Atlas / ENCODE pull, use parallel transfers — a single-stream `curl` is what
made this slow locally (measured 4.3 MB/s):

```bash
# s5cmd is much faster than aws s3 cp for many objects
curl -sL https://github.com/peak/s5cmd/releases/latest/download/s5cmd_Linux-64bit.tar.gz \
  | tar xz -C /usr/local/bin s5cmd
```

**Stop it, don't terminate it**, between sessions — you keep the environment and pay only for the
150 GB gp3 (~$12/month, and less if you shrink it).

---

## 3. Access — IAM role, never keys on the box

1. IAM → Roles → create `chip-analysis-ec2`, trusted entity **EC2**.
2. Inline policy, scoped to your prefixes rather than the whole bucket:

```json
{ "Version": "2012-10-17", "Statement": [
  { "Effect": "Allow", "Action": ["s3:ListBucket"],
    "Resource": "arn:aws:s3:::YOUR-BUCKET",
    "Condition": {"StringLike": {"s3:prefix": ["gnp-chip/*","chip-bursting/*","refs/*"]}}},
  { "Effect": "Allow", "Action": ["s3:GetObject","s3:PutObject"],
    "Resource": ["arn:aws:s3:::YOUR-BUCKET/gnp-chip/*",
                 "arn:aws:s3:::YOUR-BUCKET/chip-bursting/*",
                 "arn:aws:s3:::YOUR-BUCKET/refs/*"]}
]}
```

3. Attach at launch under *Advanced details → IAM instance profile*.
4. **Add an S3 gateway VPC endpoint** (VPC → Endpoints → `com.amazonaws.us-east-1.s3`, type
   *Gateway*). It is **free**, and without it S3 traffic from a private subnet goes through the
   NAT gateway at ~$0.045/GB — which on a 28 GB pull is small, but NAT egress is the classic way a
   cheap analysis instance produces a surprising bill.

ENCODE needs no credentials at all — `s3://encode-public` is public. Use `--no-sign-request`:

```bash
aws s3 ls --no-sign-request s3://encode-public/ --region us-west-2
```

---

## 4. Before you launch: two minutes of insurance

* **Budget alarm.** Billing → Budgets → monthly cost budget, $50, alert at 80%.
* **Tags at launch:** `Project=ezh2-bursting`, `Owner=<you>`. Makes the bill readable later.
* Know the difference: **stop** keeps (and bills) the EBS volume; **terminate** releases it.
  Instance-store NVMe is lost on both.

Access without SSH keys is worth setting up once: attach `AmazonSSMManagedInstanceCore` to the
role and use **Session Manager** — no key pair, no open port 22, no public IP.

---

## 5. Order of work

1. **Instance A**: align your GNP/MB ChIP, push BAM/bigWig/peaks to S3, **terminate**.
2. **Instance B**, Phase 1 (~400 MB, minutes): ChIP-Atlas peak BEDs for the chosen antigens,
   **wild-type experiments only** — 37% of the mESC Suz12 set is KO-flagged and filtering changed
   the Ezh2 result from 0.96× to 0.51× (`protein_noise_bursting_chromatin.md` §2c). Build the
   gene × feature matrices.
3. **The decision point: train on mESC, test on GM12878.** Needs no lab data. If held-out
   cross-cell-type R² is near zero, stop — the mapping does not transfer and your own ChIP will
   not rescue it (`protein_noise_burst_imputation.md` §5).
4. **Phase 2** only if step 3 justifies it: continuous RPM from bigWigs (~28 GB), promoter
   −1000/+100 and gene body separately, matching Ochiai.
5. Apply to your GNP/MB ChIP, prioritising the **±vismodegib / ±EZH2i** contrast — the part a
   borrowed burst size structurally cannot give you.
6. Push the derived matrices to S3 and pull them down; they are a few MB and the rest of the
   modelling is comfortable locally.

**If this becomes routine, move step 1 to AWS Batch.** It queues per-sample jobs and terminates
instances on completion, which eliminates the single largest cost risk — leaving something
running.
