#!/usr/bin/env bash
# Merge the open release pull request if the Tests workflow passed on its head.
#
# This is the weekly release. release-please.yml runs it on its Monday cron, and
# on a dispatch with `merge_release_pr`. A maintainer can always merge the
# release pull request by hand instead; that is the same release, sooner.
#
# The gate is the one release.yml applies before it publishes: a completed,
# successful run of the workflow named `Tests` on the exact commit. Advisory
# checks (complexity monitoring) do not hold a release back here either. The
# required checks on `main` still apply to the merge itself, so a red required
# check refuses it even if this script were wrong.
#
# Checks still running are not a failure. The release pull request is rebuilt
# on every releasable merge, so its checks may be minutes old at cron time. The
# script warns and leaves the merge for a dispatch or the next Monday.
#
# Writes `merged=true` or `merged=false` to $GITHUB_OUTPUT.
#
# Usage: merge-release-pr.sh   (needs GH_TOKEN and GH_REPO)
set -euo pipefail

: "${GH_REPO:?GH_REPO must name the repository}"
OUT="${GITHUB_OUTPUT:-/dev/null}"
GATE_WORKFLOW="Tests"
LABEL="autorelease: pending"

prs=$(gh pr list --base main --state open --label "${LABEL}" --json number --jq '.[].number')
count=$(printf '%s' "${prs}" | grep -c . || true)
if [ "${count}" -eq 0 ]; then
  echo "No release pull request is open: nothing releasable since the last tag."
  echo "merged=false" >>"${OUT}"
  exit 0
fi
if [ "${count}" -gt 1 ]; then
  echo "::error::${count} open pull requests carry '${LABEL}' (${prs//$'\n'/, }). Expected one; merge the right one by hand."
  exit 1
fi
pr="${prs}"

sha=$(gh pr view "${pr}" --json headRefOid --jq .headRefOid)
read -r status conclusion < <(
  gh api "repos/${GH_REPO}/actions/runs?head_sha=${sha}&per_page=100" \
    --jq "[.workflow_runs[] | select(.name == \"${GATE_WORKFLOW}\")]
          | sort_by(.created_at) | last
          | \"\(.status // \"absent\") \(.conclusion // \"none\")\""
)
echo "Release pull request #${pr} at ${sha}: ${GATE_WORKFLOW} is ${status}/${conclusion}."

if [ "${status}" = "absent" ]; then
  echo "::error::No ${GATE_WORKFLOW} run exists on release pull request #${pr} at ${sha}, so nothing will ever pass. Not merged. Close and reopen it to start its checks, then merge it once green."
  exit 1
fi
if [ "${status}" != "completed" ]; then
  echo "::warning::${GATE_WORKFLOW} has not finished on release pull request #${pr} (${status}). Not merged. Merge it by hand once green, or run: gh workflow run release-please.yml -f merge_release_pr=true"
  echo "merged=false" >>"${OUT}"
  exit 0
fi
if [ "${conclusion}" != "success" ]; then
  echo "::error::${GATE_WORKFLOW} concluded '${conclusion}' on release pull request #${pr}. Not merged."
  exit 1
fi

# --match-head-commit: if release-please pushed a new head after the check
# above, the merge is refused rather than landing an untested commit.
gh pr merge "${pr}" --squash --match-head-commit "${sha}"
echo "Merged release pull request #${pr} at ${sha}."
echo "merged=true" >>"${OUT}"
