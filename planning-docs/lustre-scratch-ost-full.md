## Lustre scratch disk-full troubleshooting (tmp/OST)

### Summary
We hit a "write error (disk full?)" when downloading NCBI datasets into `/global/scratch`, even though `df -h` reported ample free space. Root cause was Lustre OST saturation: specific OSTs were 100% full, and the default tmp directory used a layout that placed small files onto the `ddn_nvme2` pool backed by those full OSTs. The fix was to rehome `TMPDIR` onto a `ddn_hdd`-striped directory and optionally symlink the old tmp path to it.

### Symptoms
- Errors like: `write error (disk full?)` for files under `/global/scratch/users/<user>/workspace/tmp/...`
- `df -h /global/scratch` showed large free space
- Inodes were not exhausted

### Findings
- User/group quotas on `/global/scratch` were reported as `quota = 0, limit = 0` (no enforced quota).
- `lfs df -h /global/scratch` showed some OSTs at 100% (notably in the `ddn_nvme2` pool).
- `lfs getstripe -d /global/scratch/users/<user>/workspace/tmp` showed composite layout:
  - first extent on `ddn_nvme2` (small files land here)
  - later extents on `ddn_hdd`
- Small writes failed because they targeted full OSTs in `ddn_nvme2`.

### Commands used
```bash
df -h /global/scratch
df -i /global/scratch
lfs quota -u $USER /global/scratch
lfs quota -g $(id -gn) /global/scratch

lfs df -h /global/scratch
lfs df -i /global/scratch

lfs getstripe -d /global/scratch/users/$USER/workspace/tmp
lfs pool_list /global/scratch
lfs pool_list lr6.ddn_nvme2
lfs pool_list lr6.ddn_hdd
```

### Fix applied
Create a tmp directory striped explicitly on `ddn_hdd`, then point `TMPDIR` to it or symlink the old tmp path:
```bash
mkdir -p /global/scratch/users/$USER/workspace/tmp_hdd
lfs setstripe -p ddn_hdd -c 1 -S 1M /global/scratch/users/$USER/workspace/tmp_hdd
lfs getstripe -d /global/scratch/users/$USER/workspace/tmp_hdd

# One-off usage
TMPDIR=/global/scratch/users/$USER/workspace/tmp_hdd <command>

# Optional: preserve existing path
backup="/global/scratch/users/$USER/workspace/tmp_old_$(date +%Y%m%d_%H%M%S)"
mv /global/scratch/users/$USER/workspace/tmp "$backup"
ln -s /global/scratch/users/$USER/workspace/tmp_hdd /global/scratch/users/$USER/workspace/tmp
```

### Helper script
Use `bin/reset_lustre_tmp.sh` to create a new striped tmp directory and repoint the `tmp` symlink:
```bash
bin/reset_lustre_tmp.sh --pool ddn_hdd
```

### Reusable prompt
```
We are seeing "disk full" errors on /global/scratch even though df shows free space.
Please:
1) Check lfs quota (user/group) and lfs df -h /global/scratch to find full OSTs.
2) Inspect striping of the failing path with lfs getstripe.
3) If the tmp path stripes to a full pool (e.g., ddn_nvme2), create a new tmp directory
   striped on a healthier pool (e.g., ddn_hdd) and set TMPDIR or symlink.
4) Confirm write access and rerun the failing command.
```
