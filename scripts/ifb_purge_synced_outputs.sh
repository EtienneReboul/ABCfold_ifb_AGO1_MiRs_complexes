#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# ifb_purge_synced_outputs.sh — routine IFB inode/space reclaim
#
# Two independent reclaim steps on the IFB cluster:
#
#  1. results/abcfold/ + results/metadata/  — the ABCfold array leaves ~8k
#     small per-sample files (raw per-backend CIF/JSON/npz) plus the big
#     arrays.h5. Once rsynced back and consumed by postprocessing the cluster
#     copy is pure archive. Deleted ONLY after rsync proves every file is
#     already present locally at the same size (timestamp/permission-only
#     differences are ignored — they don't mean the content differs).
#
#  2. conda package cache  — `micromamba clean --all` on $MAMBA_ROOT_PREFIX:
#     drops downloaded tarballs, the repodata index, and extracted packages
#     that no environment links. NEVER touches envs/ (the 6 backend envs +
#     metadata-compress that abcfold needs stay intact).
#
#   ./scripts/ifb_purge_synced_outputs.sh            # dry-run: report only
#   ./scripts/ifb_purge_synced_outputs.sh --purge    # actually reclaim
#
# Safe to run repeatedly. If local is NOT a complete mirror of a subtree it
# is skipped and the missing files listed — rsync them down first:
#   rsync -av <IFB>:<remote>/results/<sub>/ results/<sub>/
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

IFB_HOST="${IFB_HOST:-core.cluster.france-bioinformatique.fr}"
REMOTE_DIR="${REMOTE_DIR:-/shared/projects/npf_abinitio/ABCfold_ifb_AGO1_MiRs_complexes}"
LOCAL_DIR="${LOCAL_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# micromamba root holding the shared backend envs (envs/) + package cache (pkgs/)
MAMBA_ROOT="${MAMBA_ROOT:-/shared/projects/npf_abinitio/conda}"
CONDA_PKGS="${CONDA_PKGS:-$MAMBA_ROOT/pkgs}"

# subtrees (relative to the repo root) safe to purge from IFB once mirrored here
SUBTREES=(results/abcfold results/metadata)

PURGE=0
[ "${1:-}" = "--purge" ] && PURGE=1

cd "$LOCAL_DIR"
echo "local  : $LOCAL_DIR"
echo "remote : $IFB_HOST:$REMOTE_DIR"
echo "mode   : $([ $PURGE = 1 ] && echo 'PURGE (will reclaim on IFB)' || echo 'dry-run (report only)')"
echo

remote_count() { ssh "$IFB_HOST" "find '$REMOTE_DIR' -xdev 2>/dev/null | wc -l"; }
before=$(remote_count)
echo "IFB checkout file count (start): $before"
echo

# ── step 1: mirrored-output subtrees ────────────────────────────────────────
for sub in "${SUBTREES[@]}"; do
  echo "── $sub ──────────────────────────────────────────────"
  if ! ssh "$IFB_HOST" "test -d '$REMOTE_DIR/$sub'"; then
    echo "  not present on IFB — nothing to do"; echo; continue
  fi
  mkdir -p "$sub"

  # rsync dry-run, IFB -> local. Flag files whose CONTENT differs / is missing
  # locally: itemize flag chars 3-4 contain s (size), c (checksum) or + (new).
  bad=()
  while IFS= read -r line; do
    [ -n "$line" ] && bad+=("$line")
  done < <(
    rsync -ni -a --no-perms --no-owner --no-group \
      "$IFB_HOST:$REMOTE_DIR/$sub/" "$sub/" 2>/dev/null \
    | awk '/^>f/ { f=substr($0,3,2); if (f ~ /[sc+]/) print substr($0,12) }'
  )

  if [ "${#bad[@]}" -gt 0 ]; then
    echo "  ✗ NOT a complete local mirror — ${#bad[@]} file(s) missing/differing:"
    printf '      %s\n' "${bad[@]:0:15}"
    [ "${#bad[@]}" -gt 15 ] && echo "      … and $(( ${#bad[@]} - 15 )) more"
    echo "  → pull them first:  rsync -av $IFB_HOST:$REMOTE_DIR/$sub/ $sub/"
    echo "  → skipping (nothing deleted)"; echo; continue
  fi

  rc=$(ssh "$IFB_HOST" "find '$REMOTE_DIR/$sub' -xdev | wc -l")
  rsz=$(ssh "$IFB_HOST" "du -sh '$REMOTE_DIR/$sub' | cut -f1")
  echo "  ✓ fully mirrored locally ($rc files, $rsz on IFB)"
  if [ $PURGE = 1 ]; then
    ssh "$IFB_HOST" "rm -rf '$REMOTE_DIR/$sub'"
    echo "  ✓ deleted on IFB"
  else
    echo "  (dry-run) would delete: $IFB_HOST:$REMOTE_DIR/$sub"
  fi
  echo
done

# ── step 2: conda package cache (never the envs) ────────────────────────────
echo "── conda package cache ($CONDA_PKGS) ─────────────────"
pkgs_before=$(ssh "$IFB_HOST" "find '$CONDA_PKGS' -xdev 2>/dev/null | wc -l; du -sh '$CONDA_PKGS' 2>/dev/null | cut -f1" | paste -sd' ' -)
echo "  before: $pkgs_before  (files size)"
echo "  envs are NOT touched:"
ssh "$IFB_HOST" "ls '$(dirname "$CONDA_PKGS")/envs' 2>/dev/null | sed 's/^/      /'"
if [ $PURGE = 1 ]; then
  # login shell so ~/.bashrc puts micromamba + MAMBA_ROOT_PREFIX on PATH
  ssh "$IFB_HOST" "bash -lc '
    export MAMBA_ROOT_PREFIX=\"\${MAMBA_ROOT_PREFIX:-$MAMBA_ROOT}\";
    MM=\$(command -v micromamba || echo \$HOME/.local/bin/micromamba);
    if [ -x \"\$MM\" ]; then \"\$MM\" clean --all --yes;
    elif command -v conda >/dev/null; then conda clean --all --yes;
    else echo \"    !! micromamba/conda not found — clean $CONDA_PKGS manually\"; fi'"
  pkgs_after=$(ssh "$IFB_HOST" "find '$CONDA_PKGS' -xdev 2>/dev/null | wc -l; du -sh '$CONDA_PKGS' 2>/dev/null | cut -f1" | paste -sd' ' -)
  echo "  after:  $pkgs_after  (files size)"
else
  echo "  (dry-run) would run: micromamba clean --all --yes   (MAMBA_ROOT_PREFIX=$MAMBA_ROOT)"
fi
echo

if [ $PURGE = 1 ]; then
  after=$(remote_count)
  echo "IFB checkout file count: $before → $after  (freed $(( before - after )) inodes)"
fi
