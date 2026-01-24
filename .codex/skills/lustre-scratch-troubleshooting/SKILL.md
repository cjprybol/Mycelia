---
name: lustre-scratch-troubleshooting
description: Diagnose and fix Lustre scratch "disk full" errors on /global/scratch when filesystem space appears available. Use for failed downloads/writes in scratch, TMPDIR issues, OST/MDT saturation, quota checks, and re-striping temp paths to healthy pools.
---

# Lustre Scratch Troubleshooting

## Overview

Diagnose "disk full" errors on Lustre scratch when overall free space exists. Identify quota, OST/MDT, or striping issues and redirect writes to a healthy pool.

## Quick triage

Check global space, inodes, and quotas:

```bash
df -h /global/scratch
df -i /global/scratch
lfs quota -u $USER /global/scratch
lfs quota -g $(id -gn) /global/scratch
```

Inspect OST/MDT usage:

```bash
lfs df -h /global/scratch
lfs df -i /global/scratch
```

NOTE: If specific OSTs are 100% while filesystem_summary has space, writes can fail depending on striping.

## Diagnose striping and pools

Check where the failing path is striped:

```bash
lfs getstripe -d /global/scratch/users/$USER/workspace/tmp
lfs getstripe /path/to/failing/file
```

List available pools (cluster-specific names vary):

```bash
lfs pool_list /global/scratch
lfs pool_list lr6.ddn_hdd
lfs pool_list lr6.ddn_nvme2
```

## Redirect TMPDIR to a healthy pool

Create a temp directory with explicit striping that avoids full OSTs:

```bash
mkdir -p /global/scratch/users/$USER/workspace/tmp_hdd
lfs setstripe -p ddn_hdd -c 1 -S 1M /global/scratch/users/$USER/workspace/tmp_hdd
lfs getstripe -d /global/scratch/users/$USER/workspace/tmp_hdd
```

Use it for a single command:

```bash
TMPDIR=/global/scratch/users/$USER/workspace/tmp_hdd <command>
```

If available, use the helper script in this repo:

```bash
bin/reset_lustre_tmp.sh --pool ddn_hdd
```

Optionally swap the shared temp path to a symlink (only after verifying it is safe to move):

```bash
ls -la /global/scratch/users/$USER/workspace/tmp
backup="/global/scratch/users/$USER/workspace/tmp_old_$(date +%Y%m%d_%H%M%S)"
mv /global/scratch/users/$USER/workspace/tmp "$backup"
ln -s /global/scratch/users/$USER/workspace/tmp_hdd /global/scratch/users/$USER/workspace/tmp
```

NOTE: Pool names and OST IDs vary by cluster. Prefer pools with ample free space and avoid pools tied to full OSTs.

## Validate

Confirm write access and rerun the failing command:

```bash
touch $TMPDIR/._lustre_write_test
rm -f $TMPDIR/._lustre_write_test
```

If errors persist, recheck striping on the specific path and consult cluster admins for OST health.
