#!/usr/bin/env bash
# Approve the workflow runs GitHub parked at one commit of the release pull
# request, which release-please.yml pushes with this repository's GITHUB_TOKEN.
#
# A `pull_request` event on a pull request pushed with GITHUB_TOKEN creates its
# runs in an approval-required state: `status: completed, conclusion:
# action_required`, and zero check runs. The required checks on `main` read the
# pull request's check rollup, so the release pull request would sit BLOCKED
# until someone clicked "Approve workflows to run" -- after every push, because
# release-please force-pushes the branch on each releasable merge and an
# approval covers only the commit it was given on. Dispatching the check
# workflows instead does not help: a workflow_dispatch run puts the right check
# names on the commit, but it is not associated with the pull request, so the
# rollup stays empty. Approving the parked runs makes them run as ordinary pull
# request runs, and the required checks then gate the merge exactly as on any
# other pull request. GITHUB_TOKEN may approve them; the caller's job needs
# `actions: write`.
#
# Runs are found by the head SHA, never by branch or by count: a branch query
# also returns the runs parked at every superseded head (approving one of those
# does nothing), and the number of workflows that run on a pull request changes
# whenever one is added. The runs for one event are created together, but a
# title or body edit raises its own `edited` event a moment later, so the list
# is polled until it stops growing before anything is approved.
#
# Usage: approve-parked-runs.sh <head-sha>   (needs GH_TOKEN and GH_REPO)
set -euo pipefail

SHA="${1:?usage: approve-parked-runs.sh <head-sha>}"
: "${GH_REPO:?GH_REPO must name the repository}"

POLL_SECONDS=10
STABLE_POLLS=3 # the same non-empty run list this many polls in a row
MAX_POLLS=30   # five minutes

list_runs() {
  gh api "repos/${GH_REPO}/actions/runs?head_sha=${SHA}&per_page=100" \
    --jq '.workflow_runs[] | "\(.id)\t\(.conclusion // "")\t\(.name)"' | sort
}

current=''
previous=''
stable=0
for _ in $(seq "${MAX_POLLS}"); do
  current=$(list_runs)
  ids=$(printf '%s\n' "${current}" | cut -f1)
  if [ -n "${current}" ] && [ "${ids}" = "${previous}" ]; then
    stable=$((stable + 1))
    [ "${stable}" -ge "${STABLE_POLLS}" ] && break
  else
    stable=1
  fi
  previous="${ids}"
  sleep "${POLL_SECONDS}"
done

if [ -z "${current}" ]; then
  echo "::error::No workflow run appeared at ${SHA} within $((MAX_POLLS * POLL_SECONDS)) s, so the release pull request has no checks and stays BLOCKED."
  exit 1
fi
if [ "${stable}" -lt "${STABLE_POLLS}" ]; then
  echo "::warning::The run list at ${SHA} was still changing after $((MAX_POLLS * POLL_SECONDS)) s; approving what is parked now."
fi

parked=$(printf '%s\n' "${current}" | awk -F'\t' '$2 == "action_required"')
if [ -z "${parked}" ]; then
  echo "Nothing parked at ${SHA}."
  exit 0
fi

while IFS=$'\t' read -r id _ name; do
  gh api --method POST "repos/${GH_REPO}/actions/runs/${id}/approve" --silent
  echo "Approved ${name} (run ${id}) at ${SHA}."
done <<<"${parked}"
