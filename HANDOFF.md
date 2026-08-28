# Handoff prompt — paste this into a new Claude Code session on this project

I'm continuing work on `ABCfold_ifb_AGO1_MiRs_complexes`. Read `README.md`
first for the full pipeline overview, then this note for what it doesn't cover.

## What this project is

A rebuild of `../../esther_project` (manual AF3-webserver AGO1/FBW2/miR
predictions + a PLIP-only Snakefile) using the exact protocol of the
sibling `../../drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes`:
ABCfold launches AlphaFold3 + Boltz-2 + Chai-1 + OpenFold3 + Protenix +
RosettaFold3 **together** per complex on the IFB cluster, then a local
postprocessing Snakefile does AGO1-anchored pose clustering → top-N
selection → minimize → fix_pdb → PLIP → aggregate.

Six complexes (README's table): `ago1_fbw2`, `ago1_fbw2_mir165a/168/393a`,
`ago1_fbw2_ask1_cul1`, `ago1_fbw2_ask1_cul1_mir168`. Anchor = AGO1 (chain
A) for all. RNA chains excluded from the pose-clustering feature vector
(they ride along in the minimized/PLIP structures).

## Current state: scaffold built, preprocessing test run

- All `configs/`, `worflows/`, `scripts/`, `envs/`, `notebooks/` exist.
  YAML parses, Python compiles, all 7 notebooks parse + syntax-check.
- `scripts/*.py` are byte-identical to the DRB2 pipeline's (a handful of
  docstring path examples updated). `scripts/symlink_overfolded_samples.py`
  was dropped (DRB2-specific).
- `worflows/processing/submit_abcfold.sh` = DRB2's, only `--job-name`
  changed (`abcfold_ago1_mirs_array`). Its IFB paths
  (`/shared/projects/npf_abinitio/conda/envs/…`) are unchanged and correct
  for this account — same `npf_abinitio` project, shared backend envs.
- `worflows/postprocessing/Snakefile` = DRB2's, but the third PLIP pass
  (`plip_drb2_drb4`) was dropped and `plip_rna_ligands` → `plip_mir_ligands`
  (config key + rule names), porting `esther_project/configs/mir_ligands.yaml`.
- Preprocessing (stage 1) was run locally as a validation test — see the
  session log / `logs/` for which complexes completed.

**No ABCfold job, no PLIP run, no real structure of any kind on the
cluster yet.**

## Key design decisions

- **Anchor = AGO1 (chain A)** for every complex; `partner_chains` = the
  other protein chains (`[B]` for AGO1/FBW2, `[B, C, D]` for the ASK1/CUL1
  ones). RNA is never in `partner_chains`.
- **Two PLIP passes**, mirroring `esther_project`'s two config views:
  main `plip:` = receptor (AGO1) vs. everything else in one call;
  `plip_mir_ligands:` (miR complexes only) = `--chains` omitted, no
  `--dnareceptor`, so each miR nucleotide is its own SMALLMOLECULE ligand.
- **`min_core_frac: 0.4`** — AGO1 is ~1050 aa with a Gly-rich disordered
  N-terminal ~190 aa; 0.4 keeps the folded PAZ/MID/PIWI lobe as the
  Kabsch reference.
- Notebook domain boundaries are **best-effort, not curated** — AGO1
  (N_ext/N/PAZ/MID/PIWI), FBW2 (F-box/WD40), CUL1 (coarse 2-way). Each
  notebook has `USE_WINDOWS = True` to fall back to fixed-width windows,
  and auto-detects per-chain residue offsets from the PLIP output (the
  DRB2 notebooks hardcoded `CHAIN_OFFSET`; here it's inferred).

## Open items — need a real IFB run to verify

1. **ASK1/CUL1 complex GPU/mem/time sizing not tuned.** `submit_abcfold.sh`
   defaults (`gpu:h200:1`, 250G, 2880min) are the DRB2 pipeline's
   RNA-complex profile — a reasonable start for ~2265-residue jobs but
   not measured. Ask before bumping.
2. **Chai-1 / Protenix / Boltz will likely fail on the two ASK1/CUL1
   complexes** (token caps / OOM — see README). This is an architectural
   limit, not a bug; the four smaller complexes should run all six.
3. **`mmseqs2msa` multimer behaviour** — believed correct (walks every
   `protein` entry), but inspect each `fold_input.resolved.json` before
   submitting: every protein chain should get its own MSA + ~20
   templates; RNA chains get an empty `unpairedMsa` patch (the DRB2
   pipeline patched `abcfold` on the shared IFB env for this — check
   whether that patch is still in place: `abcfold/scripts/abc_script_utils.py`
   `check_input_json()` should guard `sequence_type == "protein"` before
   injecting `templates`/`pairedMsa`).
4. **`~/.abcfold_config.ini` / backend envs** — user-scoped, already exist
   from sibling projects. Confirm (`du -sh /shared/projects/npf_abinitio/conda/envs/*`)
   before assuming `--prime` isn't needed.

## Where to pick up

Inspect the `fold_input.resolved.json` files, rsync `data/fold_inputs/` to
the IFB checkout, then `submit_abcfold.sh --test` (task 0) and check
`results/abcfold/<complex>/` per-backend CIF counts before the full array.
