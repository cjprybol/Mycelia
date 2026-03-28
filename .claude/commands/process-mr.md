# /process-mr — MR to PR Processing Pipeline

Process a merge request from the queue: rebase, create GitHub PR, run 5-persona
review, apply labels. This is the Refinery's primary workflow for the PR review
pipeline.

## Usage

```text
/process-mr [MR_ID]
```

If no MR_ID provided, claim the next ready MR from the queue.

## Instructions

### Step 1: Claim MR

```bash
export BEADS_DOLT_AUTO_START=false
export BEADS_DOLT_PORT=3308

if [ -n "$ARGUMENTS" ]; then
    MR_ID="$ARGUMENTS"
else
    # Get next ready MR
    MR_ID=$(gt refinery next mycelia 2>/dev/null | head -1 | awk '{print $1}')
    if [ -z "$MR_ID" ]; then
        echo "No MRs ready to process."
        # Exit gracefully
    fi
fi

gt refinery claim mycelia "$MR_ID"
echo "Claimed MR: $MR_ID"
```

### Step 2: Get MR details

```bash
export BEADS_DOLT_AUTO_START=false
export BEADS_DOLT_PORT=3308

bd show "$MR_ID" 2>/dev/null
```

Extract from the MR bead:

- The polecat branch name (from MR metadata or branch listing)
- The source issue ID
- The issue title (for the PR title)
- The source:td-xxx label (for cross-referencing)

### Step 3: Rebase onto master

```bash
git fetch origin master
git fetch origin <BRANCH_NAME>
git checkout <BRANCH_NAME>
git rebase origin/master
```

If rebase fails (conflict):

```bash
git rebase --abort
echo "CONFLICT: Rebase failed for $MR_ID. Spawning conflict-resolution polecat."
gt sling <ISSUE_ID> mycelia --context "Conflict resolution: rebase <BRANCH> onto master failed. Resolve conflicts and push."
gt refinery release mycelia "$MR_ID"
```

Stop processing this MR. The conflict polecat will re-submit via `gt done`.

If rebase succeeds:

```bash
git push origin <BRANCH_NAME> --force-with-lease
```

### Step 4: Create GitHub PR

```bash
PR_URL=$(gh pr create \
    --head "<BRANCH_NAME>" \
    --base master \
    --title "<ISSUE_TITLE>" \
    --body "$(cat <<'PR_BODY'
## Summary

<First 20 lines of bead description>

## Source

- GT Bead: <ISSUE_ID>
- Todo Bead: <SOURCE_TD_ID>
- MR: <MR_ID>

---
Gas Town Refinery · Mycelia rig
PR_BODY
)" \
    --label "polecat" \
    --label "review-pending")

PR_NUM=$(echo "$PR_URL" | grep -oE '[0-9]+$')
echo "Created PR #$PR_NUM: $PR_URL"
```

### Step 5: Run 5-persona review

```text
/review-pr <PR_NUM>
```

### Step 6: Handle review result

**If APPROVED** (all grades A/B):

- PR is labeled "auto-approved" by `/review-pr`
- Check auto-merge config:

  ```bash
  AUTO_MERGE=$(python3 -c "import json; print(json.loads(open('$HOME/gt/mycelia/settings/config.json').read()).get('auto_merge_approved', False))")
  if [ "$AUTO_MERGE" = "True" ]; then
      gh pr merge "$PR_NUM" --auto --squash
      echo "Auto-merge enabled for PR #$PR_NUM"
  else
      echo "PR #$PR_NUM labeled auto-approved — waiting for human merge"
  fi
  ```

- Send MERGED notification to witness:

  ```bash
  gt mail send mycelia/witness -s "MERGED <polecat-name>" -m "Branch: <BRANCH_NAME>
  Issue: <ISSUE_ID>
  PR: #$PR_NUM
  Merged-At: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  ```

- Run post-merge cleanup:

  ```bash
  gt mq post-merge mycelia "$MR_ID"
  ```

- Release the MR claim:

  ```bash
  gt refinery release mycelia "$MR_ID"
  ```

- Done. Report success.

**If CHANGES_REQUESTED** (any grade C/D/F):

- PR is labeled "changes-requested" by `/review-pr`
- Collect blocking issues from the review comment:

  ```bash
  BLOCKING=$(gh pr view "$PR_NUM" --comments --json comments \
    -q '.comments[-1].body' | grep -E 'CRITICAL|MAJOR')
  ```

- Dispatch a fix polecat:

  ```bash
  gt sling <ISSUE_ID> mycelia --context "Fix review issues on PR #$PR_NUM (branch <BRANCH_NAME>):
  $BLOCKING
  Push fixes to the SAME branch. Do NOT create a new PR."
  ```

- After fix polecat completes and pushes, re-run:

  ```text
  /review-pr <PR_NUM>
  ```

- **Max 3 fix iterations.** After 3 failed reviews:

  ```bash
  gh pr edit "$PR_NUM" --remove-label "changes-requested" \
      --add-label "needs-human-review"
  gh pr comment "$PR_NUM" --body "3 fix iterations exhausted. Escalating to human review."
  ```

### Step 7: Cleanup

```bash
echo "MR $MR_ID processed. PR #$PR_NUM"
echo "Review: <VERDICT>"
echo "Labels: $(gh pr view $PR_NUM --json labels -q '[.labels[].name] | join(", ")')"
```
