# SSH access and Claude Code on the EC2 instance

Two parts: real SSH **without opening any inbound port**, then Claude Code on the instance.

You already have an ed25519 key at `~/.ssh/id_ed25519` (created June 2025) — we reuse it rather
than making another.

---

## Part 1 — SSH over Session Manager

### Why not plain SSH

Plain SSH needs inbound port 22 open in the security group. Even restricted to your IP, that is a
public listener, and university IPs change. **SSH over SSM keeps zero inbound rules**: the
instance makes an *outbound* connection to AWS, and your SSH traffic is tunnelled through it.

You still get everything real SSH gives you — `scp`, `rsync`, VS Code Remote-SSH, port forwarding,
and a proper terminal for Claude Code.

> This revises the earlier "proceed without a key pair" advice. That was right when browser-only
> access was the plan. Now that you want SSH, importing your existing key is better.

### Step 1 — install the Session Manager plugin (on your Mac)

**You need to run this — it needs your password, which I will not handle.**

```bash
/opt/homebrew/bin/brew install --cask session-manager-plugin
```

Verify:

```bash
session-manager-plugin --version
```

### Step 2 — import your existing public key into EC2

This registers your *public* key with AWS. Your private key never leaves the Mac.

```bash
aws ec2 import-key-pair \
  --region us-east-1 \
  --key-name jpurzner-ed25519 \
  --public-key-material fileb://~/.ssh/id_ed25519.pub
```

Confirm:

```bash
aws ec2 describe-key-pairs --region us-east-1 --key-names jpurzner-ed25519
```

At launch, select **`jpurzner-ed25519`** in the *Key pair (login)* field instead of "Proceed
without a key pair". Everything else in `aws_first_instance_stepbystep.md` stays the same —
**still zero inbound rules.**

### Step 3 — configure SSH on your Mac

Create `~/.ssh/config` (you have none yet):

```
Host i-* mi-*
    User ubuntu
    IdentityFile ~/.ssh/id_ed25519
    ProxyCommand sh -c "aws ssm start-session --target %h --document-name AWS-StartSSHSession --parameters 'portNumber=%p'"
    ServerAliveInterval 60
```

```bash
chmod 600 ~/.ssh/config
```

### Step 4 — connect

Use the **instance ID**, not an IP:

```bash
ssh i-0abc123def456       # get the ID from the EC2 console
```

You should land at `ubuntu@ip-10-x-x-x`. Now `scp`, `rsync` and VS Code Remote-SSH all work:

```bash
scp results.csv i-0abc123def456:/home/ubuntu/
rsync -av ./analysis/ i-0abc123def456:/home/ubuntu/analysis/
```

### If it fails

| symptom | cause |
|---|---|
| `TargetNotConnected` | SSM agent not up, or the IAM role lacks `AmazonSSMManagedInstanceCore`. Wait for `2/2 checks passed`. |
| `session-manager-plugin not found` | Step 1 not done, or not on PATH |
| `Permission denied (publickey)` | key pair not selected at launch, or wrong `User` (Ubuntu AMIs use `ubuntu`, not `ec2-user`) |
| hangs with no output | your local AWS credentials are not configured or not us-east-1 |

---

## Part 2 — Claude Code on the instance

Requirements are met: Ubuntu 24.04 (needs 20.04+), 16–64 GB RAM (needs 4 GB+), x64.

### Install

```bash
curl -fsSL https://claude.ai/install.sh | bash
```

Then open a new shell (or `source ~/.bashrc`) and check:

```bash
claude --version        # e.g. 2.1.211 (Claude Code)
claude doctor           # read-only diagnostics
```

Installs to `~/.local/bin/claude`, **no sudo needed**, and **auto-updates in the background**.

> If `claude: command not found`, `~/.local/bin` is not on PATH:
> `echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc`

There is also a signed apt repository if you prefer package management — but it does **not**
auto-update, so for a transient instance the native installer is the better fit.

### Authenticate — the part that is awkward headless

You need a **Pro, Max, Team, Enterprise or Console** account. The free plan does not include
Claude Code.

**Interactive login works fine over SSH**, despite there being no browser on the server:

1. Run `claude`
2. It tries to open a browser and fails. **Press `c`** to copy the login URL.
3. Paste that URL into a browser **on your Mac** and sign in.
4. The browser will show a **login code** rather than redirecting — this is expected over SSH,
   because the browser cannot reach the server's local callback.
5. Paste that code back at the `Paste code here if prompted` prompt.
6. You should see `Login successful`.

Credentials are then stored at `~/.claude/.credentials.json`, mode `0600`.

**For unattended or scripted runs**, generate a long-lived token *on your Mac*:

```bash
claude setup-token          # run locally, opens a browser, prints a 1-year OAuth token
```

Then on the instance:

```bash
export CLAUDE_CODE_OAUTH_TOKEN='...'          # add to ~/.bashrc to persist
```

Caveat: that token can only make model requests — no Remote Control, no claude.ai connectors.
Locally configured MCP servers still work.

> **Do not** set `ANTHROPIC_API_KEY` unless you specifically want Console API billing. It takes
> precedence over your subscription, and if it belongs to an expired org you get confusing auth
> failures. `unset ANTHROPIC_API_KEY` and check `/status` if things look wrong.

### Two things to be aware of

**Claude Code on that box inherits the instance IAM role.** Anything it runs can reach whatever
`chip-analysis-ec2` can reach — which is exactly why that role was scoped to your project
prefixes rather than given blanket S3 access. Worth keeping it tight.

**Your credentials live on the instance.** Terminating destroys them, which is fine. But if you
ever make an AMI or snapshot of that instance, run `/logout` first so the token is not baked into
the image.

### Suggested layout

```bash
mkdir -p ~/work && cd ~/work
git clone <your repo>        # or: aws s3 sync s3://YOUR-BUCKET/... .
claude
```

Run `claude` from the project directory — it picks up `CLAUDE.md` and `.claude/` settings from
there, same as locally.

### Keeping sessions alive

An SSH drop kills your Claude Code session. Use `tmux`:

```bash
sudo apt-get install -y tmux
tmux new -s claude          # start
# detach: Ctrl-b then d
tmux attach -t claude       # reattach after reconnecting
```

Worth doing from the start on a remote box — network blips are routine and losing a long session
is annoying.
