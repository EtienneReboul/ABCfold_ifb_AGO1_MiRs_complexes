# HANDOFF — start a fresh Claude Code session in *this* folder

`ABCfold_ifb_AGO1_MiRs_complexes` — an ABCfold 6-backend rebuild of
`../../esther_project` (which did manual AF3-webserver AGO1/FBW2/miR
predictions + a PLIP-only Snakefile). Same protocol as the sibling
`../../drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes`. Read
`README.md` first for the pipeline overview.

## Current state (scaffold built 2026-08-28)

- All `configs/`, `worflows/`, `scripts/`, `envs/`, `notebooks/` exist.
  YAML parses, Python compiles, all notebooks syntax-check. Git-initialised.
- **Preprocessing validated locally**: stage 1a run for all 6 complexes
  (`data/fold_inputs/<complex>/fold_input.json` — gitignored); stage 1b
  (MMseqs2 MSA) run for `ago1_fbw2` only →
  `data/fold_inputs/ago1_fbw2/fold_input.resolved.json` (both protein
  chains got `unpairedMsa` + 20 templates, 20 seeds — correct shape).
- **Nothing has run on IFB.** No ABCfold job, no PLIP, no structures.
- The preprocessing test borrowed the DRB2 project's prebuilt abcfold
  conda env; this repo's own `envs/preprocessing.yaml` was never built.

Six complexes (`config.yaml` `complexes:`), anchor = AGO1 (chain A) for
all, RNA excluded from the pose-clustering feature vector:
`ago1_fbw2`, `ago1_fbw2_mir165a` / `_mir168` / `_mir393a`,
`ago1_fbw2_ask1_cul1`, `ago1_fbw2_ask1_cul1_mir168`.

---

## TASK A — replace the best-effort notebook domain boundaries with verified ones

`notebooks/<complex>_domain_analysis.ipynb` (6 of them) carry a
`CURATED_DOMAINS` dict with **approximate, unverified** boundaries (AGO1
N-ext/N/PAZ/MID/PIWI, FBW2 F-box/WD40, ASK1 single domain, CUL1 coarse
2-way). Each notebook also has a `USE_WINDOWS = True` fallback so the
heatmaps render regardless — but the axis labels are only meaningful once
`CURATED_DOMAINS` is right.

**Do not hand-edit the 6 notebooks.** They are generated. Edit
`CURATED_DOMAINS` (and `CHAIN_LENGTHS` if a sequence length ever changes)
at the top of `notebooks/generate_domain_notebooks.py`, then:

```bash
python notebooks/generate_domain_notebooks.py     # overwrites all 6 in place
```

### How to get the numbers

Scan the **exact** sequences in `configs/*.yaml` — not a UniProt
accession, because the esther_project constructs may be a different
isoform/truncation and you want coordinates 1-based against what is
actually folded.

Extract them:

```bash
python - <<'EOF'
import yaml, glob
seen = {}
for f in sorted(glob.glob("configs/*.yaml")):
    for s in yaml.safe_load(open(f))["sequences"]:
        if s["type"] == "protein" and s["name"] not in seen:
            seen[s["name"]] = s["sequence"]
for name, seq in seen.items():
    print(f">{name} len={len(seq)}")
    print(seq)
EOF
```

Gives AGO1 (1050 aa), FBW2 (317 aa), ASK1 (160 aa), CUL1 (738 aa).

1. **ScanProsite** — `https://prosite.expasy.org/scanprosite/` (has a REST
   API: `https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq=<SEQ>&output=json`).
   Returns every PROSITE pattern/profile hit with start–end. This is the
   direct "PROSITE-precision" source the DRB2 project used. Relevant
   profiles: PAZ `PS50821`, PIWI `PS50822` (AGO1); F-box `PS50181`, WD40
   `PS50082` / `PS50294` (FBW2). No PROSITE profile exists for the SKP1
   fold or the cullin repeats.
2. **InterProScan** — `https://www.ebi.ac.uk/interpro/` or its REST/CLI.
   Runs PROSITE **plus** Pfam + SMART + CDD and gives the consolidated
   InterPro domain call. Needed for AGO1 (PROSITE only covers PAZ/PIWI;
   the N / ArgoL1 / ArgoMID / ArgoL2 lobes are Pfam-only: `PF02170`,
   `PF08699`, `PF16486`, `PF16488`, `PF02171`) and for ASK1/CUL1.
3. **UniProt "Family & Domains"** as a cross-check —
   AGO1_ARATH `O04379`, CUL1_ARATH `Q94AH6`, ASK1/SKP1-like-1 `Q39255`.
   FBW2: look up the Arabidopsis locus on UniProt (accession not
   confirmed here).

A fresh session *may* be able to hit the ScanProsite / InterProScan REST
APIs directly (WebFetch) and populate `CURATED_DOMAINS` for you — ask it to.

---

## TASK B — run the pipeline on the IFB Core Cluster

### 0. One-time local setup

```bash
conda env create -f envs/pipeline.yaml && conda activate af3-ifb-ago1-mirs-pipeline
```

### 1. Preprocessing (local, needs internet) — finish the other 5 complexes

```bash
snakemake -s worflows/preprocessing/Snakefile --cores 2 --use-conda
```

Produces `data/fold_inputs/<complex>/fold_input.resolved.json` for all 6.
`ago1_fbw2` is already done. **Inspect each resolved JSON** before
submitting: every `protein` chain must have a non-empty `unpairedMsa` and
`templates` (~20); RNA chains get an empty `unpairedMsa` patch (see the
ABCfold patch note in step 3).

> **macOS `--use-conda` gotcha** (hit repeatedly on the DRB2 pipeline): if
> per-rule conda activation crashes with a `conda` `TypeError`
> (`_get_deactivate_scripts(None)`), launch snakemake as:
> ```bash
> env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV CONDA_SHLVL=0 \
>   PATH="$(dirname $(which snakemake)):/usr/local/bin:/usr/bin:/bin" \
>   snakemake -s worflows/... --cores N --use-conda
> ```
> `/usr/local/bin` must stay on PATH — the postprocessing `run_plip` rule
> calls bare `docker` (Docker Desktop symlink lives there).

### 2. Copy inputs to the IFB checkout

IFB Core Cluster: SSH host `core.cluster.france-bioinformatique.fr`
(`~/.ssh/config`, user `ereboul`). Account/QoS `npf_abinitio` / `normal`,
partitions `gpu` / `fast` / `long`. Shared storage
`/shared/projects/npf_abinitio/`.

Check out this repo on IFB (e.g. under
`/shared/projects/npf_abinitio/`), then from the local machine:

```bash
rsync -av data/fold_inputs/  <ifb>:/shared/projects/npf_abinitio/ABCfold_ifb_AGO1_MiRs_complexes/data/fold_inputs/
```

### 3. Submit ABCfold (on an IFB login node, from the repo root)

```bash
bash worflows/processing/submit_abcfold.sh --dry-run    # show the plan
bash worflows/processing/submit_abcfold.sh --test        # task 0 only (QoS-safe)
# inspect results/abcfold/<complex>/ , then:
bash worflows/processing/submit_abcfold.sh               # full array (skips done complexes)
```

- **`--prime` is NOT needed.** ABCfold's config (`~/.abcfold_config.ini`)
  is user-scoped on this account and the backend micromamba envs
  (`abcfold-boltz-py311`, `-chai-`, `-openfold-`, `-protenix-`,
  `-rosetta-py312`) + the `metadata-compress` env already exist under
  `/shared/projects/npf_abinitio/conda/envs/`. Sanity-check first:
  `du -sh /shared/projects/npf_abinitio/conda/envs/*` (multi-GB, not stubs)
  and `cat ~/.abcfold_config.ini`.
- The script sources `/etc/profile.d/modules.sh` itself and resolves
  `abcfold` / `micromamba` by absolute path — do NOT "simplify" those
  away (a non-interactive `ssh host 'bash ...'` / sbatch script gets
  neither `module` nor conda-env PATHs otherwise; cost real debugging
  time on the DRB2 pipeline).
- **ABCfold `check_input_json()` patch** — the DRB2 pipeline patched the
  SHARED env's `abcfold/scripts/abc_script_utils.py` `check_input_json()`
  to guard `sequence_type == "protein"` before injecting
  `templates`/`pairedMsa` (RNA chains otherwise fail schema validation
  after the empty-`unpairedMsa` patch). Check it's still in place
  (`.bak-*` backups sit alongside); if the env was rebuilt, re-apply.
- **GPU sizing** — `submit_abcfold.sh` defaults: `gpu:h200:1`, 250G,
  2880min (48h), excludes `gpu-node-7,9`. `gpu-node-4` = 4×H200 (141GB
  VRAM each, highest on the cluster). This profile is the DRB2 RNA
  complex's; fine for the 4 small AGO1/FBW2(+miR) complexes (~1367 aa,
  ~5h each expected). For the two ASK1/CUL1 complexes (~2265 aa) it's a
  reasonable start but **untuned — ask the user before bumping**.
- **Expect 3/6 backends on the two ASK1/CUL1 complexes.** At ~2265 aa,
  Chai-1 (hard 2048-token cap), Protenix (hard 2560-token cap) and Boltz
  (OOM even on H200) all fail while `abcfold` still exits 0. After those
  two runs, count per-backend CIFs in `results/abcfold/<complex>/` — do
  **not** trust the SLURM `COMPLETED` state alone. This is an
  architectural limit, not a bug; AF3 + OpenFold3 + RosettaFold3 should
  still produce structures.

### 4. Copy results back + postprocess (local)

```bash
rsync -av <ifb>:.../results/abcfold/  results/abcfold/
rsync -av <ifb>:.../results/metadata/ results/metadata/
snakemake -s worflows/postprocessing/Snakefile --cores 4 --use-conda   # + the macOS env wrapper from step 1
```

Postprocessing = compress metadata → AGO1-anchored pose clustering →
top-N/cluster selection → minimize (ChimeraX) → fix_pdb (miR complexes) →
PLIP → aggregate, **plus** a second `plip_mir_ligands` pass for miR
complexes (`--chains` omitted, no `--dnareceptor`, each miR nt its own
ligand — port of `esther_project/configs/mir_ligands.yaml`). Deliverables
per complex: `results/<complex>/all_selected_summary.csv`
(+ `..._mir_ligands.csv` for miR complexes), `pose_clusters.csv`,
`pose_clusters_pca.svg`.

`chimerax_bin` in `config.yaml` is `/Applications/ChimeraX-1.12.app/...`
— adjust for the local machine if needed.

### 5. Notebooks

```bash
conda env create -f envs/notebook.yaml
python -m ipykernel install --user --name abcfold-ago1-mirs-notebook
```

Then run each `notebooks/<complex>_domain_analysis.ipynb` and
`notebooks/pose_clustering.ipynb`. Do TASK A first if you want meaningful
domain axes.

---

## Autonomy note

The user has previously said, for the DRB2 pipeline: once a `--test` run
validates clean, they're fine with the session autonomously fixing
clear-root-cause failures (env/PATH issues confirmed by investigation, not
guesses), re-testing, and submitting the full multi-day array without
checking back — **but only for a conditional task they've explicitly
stated** ("if the test passes, launch the full run"). Ambiguous resource
tradeoffs or science/results judgement calls → stop and leave a status
trail. Don't assume this standing authorization without them stating it
for this pipeline.
