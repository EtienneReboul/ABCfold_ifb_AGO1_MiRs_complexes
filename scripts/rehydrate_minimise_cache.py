#!/usr/bin/env python3
"""
rehydrate_minimise_cache.py — carry the minimise / PLIP cache across a re-cluster
================================================================================

Adding a backend (e.g. Boltz) to results/abcfold/<complex>/ and re-running the
postprocessing pipeline makes `pose_cluster` re-fit PCA + clustering on the
whole pooled set, so cluster labels shift and `select_top_n_per_cluster`
renames every staged CIF:

    old:  selected/cluster_7/rank_75_openfold3_seed5_samplenan.cif
    new:  selected/cluster_3/rank_12_openfold3_seed5_samplenan.cif   (same model)

`minimize_cif` / `run_plip` / `plip_to_csv` outputs are keyed on {cluster}/{fname},
so Snakemake would recompute ~every model even though the input structure is
byte-identical. But a minimised PDB (and its PLIP report) depend only on the
*model* — pinned by its ABCfold `source_cif`, which is stable across re-clusters —
never on which cluster / rank slot it lands in.

Two phases, because `select_top_n_per_cluster` overwrites selected_models.csv:

  1. BEFORE re-running postprocessing:
       python scripts/rehydrate_minimise_cache.py --snapshot --all
     records, per complex,  source_cif -> (cluster, fname)  from the current
     results/<complex>/selected_models.csv into results/<complex>/.rehydrate_snapshot.tsv

  2. Re-run postprocessing ONLY up to the selection checkpoint:
       <env wrapper> snakemake -s worflows/postprocessing/Snakefile --cores 4 --use-conda \
         --until select_top_n_per_cluster \
         results/<complex>/all_selected_summary.csv ...

  3. Rehydrate — hard-link each unchanged model's cached minimise / PLIP output
     to its new {cluster}/{fname} path and bump mtimes so Snakemake sees it as
     current:
       python scripts/rehydrate_minimise_cache.py --apply --all

  4. Run the rest with mtime triggers:
       <env wrapper> snakemake -s worflows/postprocessing/Snakefile --cores 5 --use-conda \
         --rerun-triggers mtime --rerun-incomplete --keep-going \
         results/<complex>/all_selected_summary.csv \
         results/<complex>/all_selected_summary_mir_ligands.csv ...

Only the genuinely new models (Boltz) are then left for Snakemake to compute.
`--dry-run` with `--apply` reports the counts without touching anything.
Idempotent.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import time
from pathlib import Path

SNAPSHOT = ".rehydrate_snapshot.tsv"

# per-model cached artefacts:  <tree>/<cluster>/<fname>/<rel formatted with {f}=fname>
# each entry is (rel-path template, dependency layer). The layer sets the
# relinked file's mtime = now + layer*LAYER_GAP so Snakemake's `--rerun-triggers
# mtime` sees the full chain as current: cif < minimized.pdb < _fixed.pdb <
# plip report.txt < plip summary.csv (and likewise for plip_mir_ligands). Using
# a plain enumerate() index here is the classic bug — plip/report.txt then lands
# at the same mtime as minimized/pdb and *older* than _fixed.pdb, so PLIP re-runs
# for every carried-forward model.
LAYER_GAP = 4  # seconds between dependency layers (> fs mtime granularity)
ARTEFACTS = {
    # {f}_energy.csv sits next to the minimised pdb — the domain notebooks'
    # energy filter reads it; without it a carried-forward model is treated as
    # "unassessed" and always kept (so e.g. the RF3 energy-blow-up models that
    # should be dropped slip back in).
    "minimized":        [("{f}.pdb", 0), ("{f}_energy.csv", 0), ("{f}_fixed.pdb", 1)],
    "plip":             [("{f}_report/{f}_report.txt", 2), ("{f}_report/csv/summary.csv", 3)],
    "plip_mir_ligands": [("{f}_report/{f}_report.txt", 2), ("{f}_report/csv/summary.csv", 3)],
}


def read_selected_models(complex_dir: Path):
    """yield (source_cif, cluster, fname) for each row of selected_models.csv."""
    csv_path = complex_dir / "selected_models.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)
    with csv_path.open() as fh:
        for row in csv.DictReader(fh):
            yield (row["source_cif"],
                   str(row["cluster"]),
                   Path(row["staged_cif"]).stem)


def do_snapshot(complex_dir: Path) -> int:
    rows = list(read_selected_models(complex_dir))
    out = complex_dir / SNAPSHOT
    with out.open("w") as fh:
        for src, cl, fn in rows:
            fh.write(f"{src}\t{cl}\t{fn}\n")
    print(f"  snapshot: {len(rows)} models -> {out}")
    return len(rows)


def load_snapshot(complex_dir: Path) -> dict[str, tuple[str, str]]:
    snap = complex_dir / SNAPSHOT
    if not snap.is_file():
        raise FileNotFoundError(
            f"{snap} missing — run `--snapshot` before the postprocessing re-run")
    m: dict[str, tuple[str, str]] = {}
    for line in snap.read_text().splitlines():
        if not line.strip():
            continue
        src, cl, fn = line.split("\t")
        m[src] = (cl, fn)
    return m


def link(src: Path, dst: Path, dry: bool) -> bool:
    """hard-link src -> dst (copy fallback). True if a new link was made."""
    if dst.exists():
        return False
    if dry:
        return True
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    return True


def do_apply(complex_dir: Path, dry: bool) -> dict[str, int]:
    old = load_snapshot(complex_dir)
    new_rows = list(read_selected_models(complex_dir))

    n_total = len(new_rows)
    n_new = n_unchanged = n_relinked = n_files = n_orphan_src = 0
    now = time.time()

    for src, new_cl, new_fn in new_rows:
        if src not in old:
            n_new += 1
            continue
        old_cl, old_fn = old[src]
        if (old_cl, old_fn) == (new_cl, new_fn):
            n_unchanged += 1
            continue

        moved = False
        for tree, rels in ARTEFACTS.items():
            old_dir = complex_dir / tree / old_cl / old_fn
            new_dir = complex_dir / tree / new_cl / new_fn
            if not old_dir.is_dir():
                continue
            for rel, layer in rels:
                s = old_dir / rel.format(f=old_fn)
                if not s.is_file():
                    continue
                d = new_dir / rel.format(f=new_fn)
                if link(s, d, dry):
                    n_files += 1
                    moved = True
                    if not dry:
                        t = now + layer * LAYER_GAP
                        os.utime(d, (t, t))
        if moved:
            n_relinked += 1
        else:
            # in snapshot, moved slot, but no cached artefact on disk to carry
            n_orphan_src += 1

    return dict(total=n_total, new=n_new, unchanged=n_unchanged,
                relinked=n_relinked, files=n_files, orphan=n_orphan_src)


def resolve_complexes(root: Path, names, want_all: bool):
    if names:
        return [root / n for n in names]
    if want_all:
        return sorted(p.parent for p in root.glob("*/selected_models.csv"))
    return []


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--postproc-root", default="results")
    ap.add_argument("--complex", action="append", dest="complexes",
                    help="complex name; repeatable")
    ap.add_argument("--all", action="store_true",
                    help="every <root>/*/selected_models.csv")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--snapshot", action="store_true",
                      help="phase 1: record the current selection")
    mode.add_argument("--apply", action="store_true",
                      help="phase 3: relink cached outputs to the new paths (default)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(args.postproc_root)
    complexes = resolve_complexes(root, args.complexes, args.all)
    if not complexes:
        ap.error("pass --complex <name> (repeatable) or --all")
    do_snap = args.snapshot and not args.apply

    grand: dict[str, int] = {}
    for cdir in complexes:
        print(f"\n=== {cdir.name} ===")
        if not cdir.is_dir():
            print("  !! not found"); continue
        try:
            if do_snap:
                do_snapshot(cdir)
                continue
            stats = do_apply(cdir, args.dry_run)
        except FileNotFoundError as e:
            print(f"  !! {e}"); continue
        print(f"  selected models        : {stats['total']}")
        print(f"  carried forward        : {stats['relinked'] + stats['unchanged']} "
              f"(relinked {stats['relinked']} / {stats['files']} files, "
              f"unchanged {stats['unchanged']})")
        print(f"  new -> compute          : {stats['new']}")
        if stats['orphan']:
            print(f"  moved but no cache     : {stats['orphan']}  (Snakemake will compute)")
        for k, v in stats.items():
            grand[k] = grand.get(k, 0) + v

    if not do_snap and grand:
        carried = grand.get('relinked', 0) + grand.get('unchanged', 0)
        print(f"\n{'DRY RUN — ' if args.dry_run else ''}TOTAL: "
              f"{carried}/{grand.get('total', 0)} carried forward "
              f"({grand.get('files', 0)} files linked), "
              f"{grand.get('new', 0) + grand.get('orphan', 0)} left for Snakemake.")
        if not args.dry_run:
            print("Next: snakemake … --rerun-triggers mtime --rerun-incomplete --keep-going")


if __name__ == "__main__":
    main()
