#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# ifb_purge_synced_outputs.sh — routine IFB inode/space reclaim
#
# The ABCfold array on the IFB cluster leaves ~8k small per-sample files in
# results/abcfold/ (raw per-backend CIF/JSON/npz) plus the big results/metadata/
# arrays.h5. Once those are rsynced back here and postprocessing has consumed
# them, the cluster copy is pure archive and just eats the shared npf_abinitio
# inode quota. This script deletes a remote subtree ONLY after proving, with
# rsync's own file-by-file comparison, that every file in it is already present
# locally at the same size (timestamp / permission differences are ignored —
# they don't mean the content differs).
#
#   ./scripts/ifb_purge_synced_outputs.sh            # dry-run: report only
#   ./scripts/ifb_purge_synced_outputs.sh --purge    # actually delete on IFB
#
# Safe to run repeatedly. If local is NOT a complete mirror of a subtree, that
# subtree is skipped and the missing/differing files are listed — rsync them
# down first:  rsync -av <IFB>:<remote>/results/<sub>/ results/<sub>/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

IFB_HOST="${IFB_HOST:-core.cluster.france-bioinformatique.fr}"
REMOTE_DIR="${REMOTE_DIR:-/shared/projects/npf_abinitio/ABCfold_ifb_AGO1_MiRs_complexes}"
LOCAL_DIR="${LOCAL_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# subtrees (relative to the repo root) that are safe to purge from IFB once
# fully mirrored here. Add more as the pipeline grows.
SUBTREES=(results/abcfold results/metadata)

PURGE=0
[[ "${1:-}" == "--purge" ]] && PURGE=1

cd "$LOCAL_DIR"
echo "local  : $LOCAL_DIR"
echo "remote : $IFB_HOST:$REMOTE_DIR"
echo "mode   : $([[ $PURGE == 1 ]] && echo 'PURGE (will delete on IFB)' || echo 'dry-run (report only)')"
echo

remote_count() { ssh "$IFB_HOST" "find '$REMOTE_DIR' -xdev 2>/dev/null | wc -l"; }
before=$(remote_count)
echo "IFB checkout file count (start): $before"
echo

purged_any=0
for sub in "${SUBTREES[@]}"; do
  echo "── $sub ──────────────────────────────────────────────"
  if ! ssh "$IFB_HOST" "test -d '$REMOTE_DIR/$sub'"; then
    echo "  not present on IFB — nothing to do"; echo; continue
  fi
  mkdir -p "$sub"

  # rsync dry-run, IFB -> local. Flag only files whose CONTENT differs or that
  # are missing locally: itemize flags with s (size), c (checksum) or + (new).
  # (bash 3.2 — no mapfile; read into an array.)
  bad=()
  while IFS= read -r line; do
    [ -n "$line" ] && bad+=("$line")
  done < <(
    rsync -ni -a --no-perms --no-owner --no-group \
      "$IFB_HOST:$REMOTE_DIR/$sub/" "$sub/" 2>/dev/null \
    | awk '/^>f/ { f=substr($0,3,2); if (f ~ /[sc+]/) print substr($0,12) }'
  )

  if (( ${#bad[@]} > 0 )); then
    echo "  ✗ NOT a complete local mirror — ${#bad[@]} file(s) missing/differing:"
    printf '      %s\n' "${bad[@]:0:15}"
    (( ${#bad[@]} > 15 )) && echo "      … and $(( ${#bad[@]} - 15 )) more"
    echo "  → pull them first:  rsync -av $IFB_HOST:$REMOTE_DIR/$sub/ $sub/"
    echo "  → skipping (nothing deleted)"; echo; continue
  fi

  rc=$(ssh "$IFB_HOST" "find '$REMOTE_DIR/$sub' -xdev | wc -l")
  rsz=$(ssh "$IFB_HOST" "du -sh '$REMOTE_DIR/$sub' | cut -f1")
  echo "  ✓ fully mirrored locally ($rc files, $rsz on IFB)"
  if (( PURGE == 1 )); then
    ssh "$IFB_HOST" "rm -rf '$REMOTE_DIR/$sub'"
    echo "  ✓ deleted on IFB"
    purged_any=1
  else
    echo "  (dry-run) would delete: $IFB_HOST:$REMOTE_DIR/$sub"
  fi
  echo
done

if (( PURGE == 1 && purged_any == 1 )); then
  after=$(remote_count)
  echo "IFB checkout file count: $before → $after  (freed $(( before - after )) inodes)"
fi
