# Installing the AWS CLI

## On the EC2 instance (Ubuntu) — what you need

**Do not use `apt-get install awscli`.** The Ubuntu archive package has historically been AWS CLI
**v1**, which is deprecated and behaves differently from v2 (different output defaults, no SSO
support, missing subcommands). Use Amazon's official v2 installer:

```bash
sudo apt-get update
sudo apt-get install -y unzip
curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip -q awscliv2.zip
sudo ./aws/install
aws --version        # must say 2.x -- if it says 1.x, an apt version is shadowing it
```

`x86_64` is correct for the `c6i.2xlarge` / `c6id.8xlarge` instances chosen — both Intel. On a
Graviton (ARM) instance it would be `awscli-exe-linux-aarch64.zip`.

If you already ran `apt-get install awscli` and got v1, remove it first so it does not shadow v2:

```bash
sudo apt-get remove -y awscli && hash -r
```

### The full setup line for each instance

**Analysis box** (`c6i.2xlarge`):

```bash
sudo apt-get update
sudo apt-get install -y python3-pip bedtools unzip
pip3 install --break-system-packages pyBigWig numpy pandas scikit-learn pyarrow

curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip -q awscliv2.zip && sudo ./aws/install

curl -sL https://github.com/peak/s5cmd/releases/latest/download/s5cmd_Linux-64bit.tar.gz \
  | sudo tar xz -C /usr/local/bin s5cmd
```

**Alignment box** (`c6id.8xlarge`):

```bash
sudo apt-get update
sudo apt-get install -y bowtie2 samtools bedtools macs2 fastp pigz unzip

curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip -q awscliv2.zip && sudo ./aws/install
```

### Do NOT run `aws configure` on EC2

This is the part people get wrong. The **IAM instance role** (`chip-analysis-ec2`) supplies
credentials to the instance automatically, rotated by AWS. Running `aws configure` there would
write long-lived keys to disk that *override* the role — strictly worse.

Verify the role is working instead:

```bash
aws sts get-caller-identity
```

Expect an ARN containing **`assumed-role/chip-analysis-ec2/i-0abc...`**. If you instead see an
IAM *user* ARN, something wrote static credentials and they are taking precedence — check
`~/.aws/credentials` and delete it.

Then confirm actual access:

```bash
aws s3 ls s3://YOUR-BUCKET/                                        # your data (needs the role)
aws s3 ls --no-sign-request s3://encode-public/ --region us-west-2 # public, no credentials
```

`--no-sign-request` is required for ENCODE: it is a public bucket, and signing the request with
your credentials will fail rather than succeed.

---

## On your Mac — already done

**`apt-get` does not exist on macOS** — it is Debian/Ubuntu's package manager. The macOS
equivalent is Homebrew, and the CLI is already installed here:

```
aws-cli/2.37.5 Python/3.14.7 Darwin/24.5.0 source/arm64
```

Installed with `/opt/homebrew/bin/brew install awscli`, already on your `PATH`.

### Why the explicit `/opt/homebrew/` path

You have **two Homebrew installations**:

| path | architecture | PATH position |
|---|---|---|
| `/usr/local/bin/brew` | **Intel, under Rosetta 2** | 4 — *first* |
| `/opt/homebrew/bin/brew` | **native ARM** | 18 |

Plain `brew install awscli` would have used the Intel one and produced an emulated x86 binary.
Because `/usr/local/bin` comes first, anything installed with the Intel brew shadows the ARM
version of the same tool. Prefer `/opt/homebrew/bin/brew` for new installs; to make it the
default, add to `~/.zshrc`:

```bash
eval "$(/opt/homebrew/bin/brew shellenv)"
```

### Credentials on the Mac — you do this, not me

I have not configured credentials and will not; access keys are yours to handle.

```bash
aws configure          # IAM access keys: ID, secret, region us-east-1, output json
aws configure sso      # better if your institution uses Identity Center -- short-lived tokens
```

`aws configure` writes plain-text keys to `~/.aws/credentials`. Prefer SSO where available.

### Useful before any instance exists

Find which region each bucket is in — this decides where compute should go:

```bash
for b in $(aws s3api list-buckets --query 'Buckets[].Name' --output text); do
  printf "%-40s %s\n" "$b" "$(aws s3api get-bucket-location --bucket $b --query LocationConstraint --output text)"
done
```

> `None` means **us-east-1** — an old API quirk, not an error.
