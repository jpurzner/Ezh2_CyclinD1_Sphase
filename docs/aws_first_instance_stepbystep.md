# Launching your first EC2 instance — step by step

Companion to `aws_setup_walkthrough.md`, written for a first-timer. Region is **us-east-1**
throughout.

**Start with the small box.** Learn the mechanics on the $0.34/hr analysis instance, not the
$1.61/hr alignment one. Everything transfers, and mistakes are cheap. The alignment instance is
§7 and only differs in four fields.

AWS changes console wording occasionally — if a label reads slightly differently, the meaning
below should still identify it.

---

## 0. Two things before you touch EC2

### 0a. Set a budget alarm (2 minutes, saves real money)

The most common beginner cost event is a forgotten running instance. This catches it.

1. Top-right account menu → **Billing and Cost Management**
2. Left sidebar → **Budgets** → **Create budget**
3. **Use a template** → **Monthly cost budget**
4. Amount: **$50**. Email: yours. → **Create budget**

### 0b. Check your region — this trips up nearly everyone

Top-right of the console, next to your account name, is a region selector. **Set it to
`US East (N. Virginia) us-east-1`.**

EC2 is region-scoped: an instance launched in the wrong region is invisible from every other
region, and people routinely "lose" instances and keep paying for them. If a resource you made
seems to have vanished, check the region selector first.

---

## 1. Create the IAM role (do this before launching)

This is what lets the instance read and write your S3 buckets **without putting access keys on
the machine**. It also enables browser-based login, so you need no SSH key at all.

1. Search bar → **IAM** → left sidebar **Roles** → **Create role**
2. **Trusted entity type**: *AWS service*. **Use case**: **EC2**. → **Next**
3. On the permissions page, search and tick **`AmazonSSMManagedInstanceCore`**
   *(this is what enables browser login — do not skip it)* → **Next**
4. **Role name**: `chip-analysis-ec2` → **Create role**

Now add your S3 access:

5. Open the role you just made → **Permissions** tab → **Add permissions** → **Create inline
   policy** → **JSON** tab
6. Paste this, **replacing `YOUR-BUCKET`** with your actual bucket name:

```json
{ "Version": "2012-10-17", "Statement": [
  { "Effect": "Allow", "Action": ["s3:ListBucket"],
    "Resource": "arn:aws:s3:::YOUR-BUCKET" },
  { "Effect": "Allow", "Action": ["s3:GetObject","s3:PutObject"],
    "Resource": "arn:aws:s3:::YOUR-BUCKET/*" }
]}
```

7. **Next** → name it `s3-project-access` → **Create policy**

> Scoping this to specific prefixes is better practice, but bucket-wide is fine to start and one
> less thing to get wrong on your first attempt. Tighten it once things work.

---

## 2. Launch the instance

Search bar → **EC2** → **Instances** → orange **Launch instances** button.

Work down the form:

### Name and tags
**Name:** `bursting-analysis`

Click **Add additional tags** → **Add tag**: Key `Project`, Value `ezh2-bursting`. This is how the
bill becomes readable later.

### Application and OS Images (AMI)
Select **Ubuntu**. Confirm the description says **Ubuntu Server 24.04 LTS**, **64-bit (x86)**, and
is marked *Free tier eligible*.

> **Not Arm/Graviton.** The x86 option avoids checking whether every bioinformatics binary has an
> Arm build. The savings are not worth the debugging.

### Instance type
Click the dropdown, type `c6i.2xlarge`, select it. (8 vCPU, 16 GiB.)

> It will warn this is not free-tier. That is expected — ~$0.34/hr, and you will stop it when idle.

### Key pair (login)
Select **Proceed without a key pair (Not recommended)**.

> The warning is aimed at people with no other access route. You have one: the SSM role from §1
> gives browser-based login, which is *more* secure than a key file — nothing to lose, nothing
> to leak, and no open SSH port.

### Network settings
Click **Edit**.

- **VPC**: your default VPC
- **Subnet**: *No preference*
- **Auto-assign public IP**: **Enable** ← required for the browser-login agent to reach AWS
- **Firewall (security groups)**: **Create security group**
  - Name: `chip-analysis-sg`
  - **Delete every inbound rule**, including the SSH one offered by default. You want **zero
    inbound rules.**

> Zero inbound is correct and is the main security win here. Browser login works over an outbound
> connection the instance makes, so nothing needs to be open to the internet. The default "allow
> SSH from 0.0.0.0/0" exposes port 22 to the entire internet and should never be accepted.

### Configure storage
Change the root volume from the default **8 GiB** to **150 GiB**, type **gp3**.

> **This is the single most common beginner mistake.** 8 GiB fills almost immediately once you
> start downloading data, and you discover it hours in. Growing a volume later is possible but
> fiddly. 150 GiB of gp3 is about $12/month, prorated — pennies for a few days.

### Advanced details
Expand it. You only need one field, near the top:

- **IAM instance profile**: select **`chip-analysis-ec2`**

> Easy to miss, and without it the instance cannot reach S3 *and* browser login will not work.
> It can be attached later, but it is one click here.

### Launch
Right-hand **Summary** panel → **Launch instance**.

---

## 3. Connect

1. **Instances** → wait for **Instance state** `Running` and **Status check** `2/2 checks passed`
   *(the second one takes 2–3 minutes; browser login will not work until it passes)*
2. Tick your instance → **Connect** button (top right)
3. Choose the **Session Manager** tab → **Connect**

A terminal opens in the browser. If Session Manager is greyed out or fails, see §6.

First, become a normal user with a home directory:

```bash
sudo su - ubuntu
```

> Session Manager drops you in as `ssm-user` in `/usr/bin`, which is disorienting. This puts you
> where you expect to be.

---

## 4. Verify before installing anything

```bash
df -h /                       # want ~150G, not 8G
nproc                         # want 8
free -g                       # want ~15
aws sts get-caller-identity   # should show .../chip-analysis-ec2/i-xxxx
aws s3 ls s3://YOUR-BUCKET/   # should list your bucket
```

If `get-caller-identity` errors or `s3 ls` says *Access Denied*, the role is not attached or the
policy is wrong — fix that before going further. Everything downstream depends on it.

Then set up:

```bash
sudo apt-get update
sudo apt-get install -y python3-pip bedtools unzip
# AWS CLI v2 -- `apt-get install awscli` may install v1; use the official installer
curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip -q awscliv2.zip && sudo ./aws/install && aws --version   # expect 2.x
pip3 install --break-system-packages pyBigWig numpy pandas scikit-learn pyarrow

# faster bulk S3 transfers than the aws cli
curl -sL https://github.com/peak/s5cmd/releases/latest/download/s5cmd_Linux-64bit.tar.gz \
  | sudo tar xz -C /usr/local/bin s5cmd

# confirm the public ENCODE bucket is reachable with no credentials
aws s3 ls --no-sign-request s3://encode-public/ --region us-west-2 | head
```

---

## 5. Shutting down — the part that actually costs money

**Stop** — keeps the disk and everything on it, no compute charges. You still pay ~$12/month for
the 150 GiB volume. Use this between working sessions.
> Instances → tick → **Instance state** → **Stop instance**

**Terminate** — deletes the instance and its disk. **Irreversible.** Use when finished entirely.
> Instances → tick → **Instance state** → **Terminate instance**

Rules of thumb:
- Stop it **every time you finish for the day**. Nothing survives in memory anyway.
- A stopped instance's public IP changes when you restart it. Harmless with Session Manager.
- After terminating, check **Elastic Block Store → Volumes** for anything left in `available`
  state. Orphaned volumes bill forever and are invisible on the Instances page.

Check spend: **Billing and Cost Management → Cost Explorer**. It lags ~24 h.

---

## 6. If something goes wrong

| symptom | cause | fix |
|---|---|---|
| **Connect** → Session Manager greyed out | role missing, or agent not up yet | Wait for `2/2 checks passed`. Then: tick instance → **Actions → Security → Modify IAM role** → attach `chip-analysis-ec2` → reboot. |
| Session Manager "not connected" | no route to the internet | Instance needs a public IP in a public subnet. Check **Auto-assign public IP** was Enable. |
| `aws s3 ls` → *Access Denied* | policy bucket name wrong | IAM → Roles → your role → check the inline policy has your real bucket name |
| `No space left on device` | root volume left at 8 GiB | Volumes → select → **Modify volume** → 150 GiB, then `sudo growpart /dev/nvme0n1 1 && sudo resize2fs /dev/nvme0n1p1` |
| Instance not in the list | wrong region | Set region selector to **us-east-1** |
| `pip3 install` → *externally-managed-environment* | Ubuntu 24.04 policy | use `--break-system-packages` as above, or a venv |

---

## 7. The alignment instance, once you are comfortable

Same procedure, four changes:

| field | value |
|---|---|
| Name | `gnp-chip-align` |
| Instance type | **`c6id.8xlarge`** (32 vCPU, 64 GiB, 1.9 TB NVMe) |
| Root storage | 100 GiB gp3 *(data goes on the NVMe, not here)* |
| everything else | identical — same role, zero inbound, public IP enabled |

**The NVMe drive is not mounted for you.** After connecting:

```bash
lsblk                                   # find the ~1.7T disk, usually nvme1n1
sudo mkfs.xfs /dev/nvme1n1
sudo mkdir -p /scratch && sudo mount /dev/nvme1n1 /scratch
sudo chown ubuntu:ubuntu /scratch
df -h /scratch                          # confirm ~1.7T
```

> **`/scratch` is erased when the instance stops — not just when it terminates.** This surprises
> people badly. Push every finished BAM to S3 as it completes, and never leave the only copy of
> anything there.

Because of that, for this instance **terminate rather than stop** when the run is done — there is
nothing on it worth keeping, and a stopped 32-core instance still bills its EBS volume.

---

## Appendix — the CLI equivalent

Once you have done it through the console once, this is the same launch in one command, useful for
repeatability:

```bash
aws ec2 run-instances \
  --region us-east-1 \
  --image-id resolve:ssm:/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id \
  --instance-type c6i.2xlarge \
  --iam-instance-profile Name=chip-analysis-ec2 \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":150,"VolumeType":"gp3"}}]' \
  --metadata-options 'HttpTokens=required' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=bursting-analysis},{Key=Project,Value=ezh2-bursting}]'
```

Verify it against a console launch before trusting it — the AMI resolver path in particular is
worth confirming returns a current image.
