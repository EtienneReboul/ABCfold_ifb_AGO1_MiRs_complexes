# ABCfold IFB AGO1 / FBW2 / miR Complex Pipeline

Sibling of
[`../../drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes`](../../drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes)
and a rebuild of [`../../esther_project`](../../esther_project) (which
predicted AGO1 / FBW2 / miR complexes by hand-submitting AF3-webserver
replicas and running a PLIP-only Snakefile — no ABCfold, no pose
clustering). Same automated protocol as the DRB2 pipeline:
[ABCfold](https://github.com/rigdenlab/ABCFold) launches AlphaFold3,
Boltz-2, Chai-1, OpenFold3, Protenix and RosettaFold3 **together** per
complex on the IFB Core Cluster; the resulting ensemble is clustered by
rigid-anchor structural alignment (anchor = AGO1); the models per cluster
are minimized and run through PLIP to detect interactions.

## Why ABCfold instead of resampling AF3 alone

The sibling `ABCfold_NPF_pipeline` project found that AF3 alone, even with
30 seeds, does not recover every conformation a second,
architecturally-different model (Boltz-2) finds for the same sequences.
This pipeline applies the same fix: run six independent architectures per
complex instead of just resampling one.

## Complexes modelled

Six, all carried over from `esther_project`'s `data/` tree. Chains are
lettered proteins-first, then any miR RNA last.

| Complex | Chains | ~Residues | Anchor | miR-ligand PLIP pass |
|---|---|---|---|---|
| `ago1_fbw2` | AGO1 (A), FBW2 (B) | 1367 | AGO1 | — |
| `ago1_fbw2_mir165a` | AGO1 (A), FBW2 (B), miR165a (C) | 1367 + 21 nt | AGO1 | yes |
| `ago1_fbw2_mir168` | AGO1 (A), FBW2 (B), miR168 (C) | 1367 + 21 nt | AGO1 | yes |
| `ago1_fbw2_mir393a` | AGO1 (A), FBW2 (B), miR393a (C) | 1367 + 22 nt | AGO1 | yes |
| `ago1_fbw2_ask1_cul1` | AGO1 (A), FBW2 (B), ASK1 (C), CUL1 (D) | 2265 | AGO1 | — |
| `ago1_fbw2_ask1_cul1_mir168` | AGO1 (A), FBW2 (B), ASK1 (C), CUL1 (D), miR168 (E) | 2265 + 21 nt | AGO1 | yes |

Protein/RNA sequences are copied verbatim from `esther_project`'s
AF3-webserver `*_job_request.json` files — see `configs/ago1_fbw2.yaml`
etc. AGO1 = *Arabidopsis thaliana* ARGONAUTE 1 (1050 aa), FBW2 =
F-box/WD40 protein FBW2 (317 aa), ASK1 = SKP1-like (160 aa), CUL1 =
CULLIN 1 (738 aa).

> **Backend token limits** — the two ASK1/CUL1 complexes are ~2265–2286
> residues, close to the size where the DRB2 pipeline's RNA complex hit
> hard architectural caps: **Chai-1** (2048 tokens), **Protenix** (2560
> tokens) and **Boltz** (OOM even on an H200). `abcfold` swallows
> per-backend failures and still exits 0, so after those two runs check
> per-backend CIF counts, not just the SLURM state. The four
> AGO1/FBW2(+miR) complexes are ~1367–1389 residues and should clear
> every backend.

---

## Pipeline overview

```text
┌───────────────────────────────────────────────────────────────┐
│  PRE-PROCESSING (local, needs internet)                        │
│  worflows/preprocessing/Snakefile                                │
│                                                                 │
│  1a. Fold input   fold_input.json per complex                    │
│                   (AlphaFold3-dialect JSON, ABCfold's native       │
│                   input — every chain, every seed)                  │
│  1b. MMseqs2      MSA + top-hit templates from the ColabFold        │
│                   webserver via ABCfold's own `mmseqs2msa` CLI,      │
│                   embedded into fold_input.resolved.json             │
└──────────────────────────┬──────────────────────────────────────┘
                           │  rsync data/fold_inputs/
                           v
┌───────────────────────────────────────────────────────────────┐
│  PROCESSING (IFB cluster)  worflows/processing/submit_abcfold.sh   │
│                                                                     │
│  2. ABCfold run — one `abcfold -abcopr ...` per complex,             │
│     AlphaFold3 + Boltz-2 + Chai-1 + OpenFold3 + Protenix +            │
│     RosettaFold3 TOGETHER on the same fold_input.resolved.json.        │
│     Each array task compresses its own confidence sprawl to            │
│     model_metadata.parquet before rsync.                               │
└──────────────────────────┬──────────────────────────────────────┘
                           │  rsync results/abcfold/, results/metadata/
                           v
┌───────────────────────────────────────────────────────────────┐
│  POST-PROCESSING (local)  worflows/postprocessing/Snakefile        │
│                                                                     │
│  3a. compress_abcfold_metadata (local fallback)                      │
│  3b. pose_cluster    AGO1 rigid-core Kabsch + hierarchical RMSD        │
│                      clustering of the partner protein chain(s),       │
│                      pooled across all backends x seeds (RNA excluded  │
│                      from the feature vector)                          │
│  3c. select_top_n_per_cluster                                          │
│  3d-3h. minimize (ChimeraX) -> fix_pdb (miR complexes) -> PLIP          │
│         (Docker) -> pliparser -> aggregate                              │
│         + a second `plip_mir_ligands` pass for miR complexes            │
│           (--chains omitted, no --dnareceptor: every miR nt is its own  │
│            ligand — port of esther_project/configs/mir_ligands.yaml)    │
└───────────────────────────────────────────────────────────────┘
```

---

## Directory layout

```text
ABCfold_ifb_AGO1_MiRs_complexes/
├── config.yaml                        ← shared defaults (tracked by git)
├── config.local.yaml.example          ← optional personal overrides
├── envs/                              ← per-rule conda envs (pipeline, notebook, …)
├── configs/
│   ├── ago1_fbw2.yaml
│   ├── ago1_fbw2_mir165a.yaml  / _mir168.yaml / _mir393a.yaml
│   ├── ago1_fbw2_ask1_cul1.yaml
│   └── ago1_fbw2_ask1_cul1_mir168.yaml
├── worflows/
│   ├── preprocessing/Snakefile        ← stage 1 (local)
│   ├── processing/submit_abcfold.sh   ← stage 2 SLURM submission (cluster)
│   └── postprocessing/Snakefile       ← stage 3 (local)
├── scripts/                           ← copied verbatim from the DRB2 pipeline
│   ├── make_multimer_af3_input.py     ← stage 1a
│   ├── fetch_mmseqs2_msa.py           ← stage 1b
│   ├── abcfold_backends.py / compress_abcfold_metadata.py / parquet_utils.py
│   ├── pose_cluster_anchor.py         ← stage 3b
│   ├── select_top_n_per_cluster.py    ← stage 3c
│   └── sanitize_cif.py / minimize_cif.py / fix_pdb.py / aggregate_summaries.py / dssp_summary.py
├── notebooks/
│   ├── <complex>_domain_analysis.ipynb  ← one per complex — pooled / per-cluster /
│   │                                       per-backend domain-contact heatmaps +
│   │                                       residue-level interface map
│   └── pose_clustering.ipynb            ← interactive PCA / GMM / HDBSCAN on pose_clusters.csv
└── data/, results/, logs/             ← created automatically (gitignored)
```

Every `scripts/*.py` is byte-identical to the DRB2 pipeline's copy (only a
few docstring path examples were updated). That project confirmed each
ABCfold backend's output layout against real completed IFB runs — see its
module docstrings.

---

## Quick start

### 1. Install the controller environment (once)

```bash
conda env create -f envs/pipeline.yaml
conda activate af3-ifb-ago1-mirs-pipeline
```

### 2. Pre-processing (local, needs internet)

```bash
snakemake -s worflows/preprocessing/Snakefile --cores 2 --use-conda
```

Produces `data/fold_inputs/<complex>/fold_input.resolved.json` for all six
complexes. Inspect them by hand against `esther_project`'s
`*_job_request.json` sequences before submitting to the cluster.

### 3. Processing (IFB cluster)

```bash
# rsync data/fold_inputs/ to the IFB checkout first
bash worflows/processing/submit_abcfold.sh --dry-run    # show the plan
bash worflows/processing/submit_abcfold.sh --test        # single-task test
bash worflows/processing/submit_abcfold.sh               # full array
```

ABCfold's config (`~/.abcfold_config.ini`) is **user-scoped** on the IFB
`npf_abinitio` account and the backend micromamba envs already exist from
the sibling projects — so `--prime` is normally **not** needed here (check
`~/.abcfold_config.ini` and `du -sh /shared/projects/npf_abinitio/conda/envs/*`
first). Then `rsync` `results/abcfold/` and `results/metadata/` back.

### 4. Post-processing (local)

```bash
snakemake -s worflows/postprocessing/Snakefile --cores 4 --use-conda
```

> **macOS controller-env gotcha** (inherited from the DRB2 pipeline): if
> `--use-conda` per-rule activation crashes with a `conda` `TypeError`,
> launch as
> `env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV CONDA_SHLVL=0 PATH="<controller bin>:/usr/local/bin:/usr/bin:/bin" snakemake …`
> (`/usr/local/bin` must stay on PATH — the `run_plip` rule calls bare `docker`).

Produces, per complex, under `results/<complex>/`:
`pose_clusters.csv`, `rmsf_profile.svg`, `pose_clusters_pca.svg`,
`selected/cluster_<k>/…cif`, `minimized/`, `plip/`,
`all_selected_summary.csv` (+ `all_selected_summary_mir_ligands.csv` for
miR complexes).

### 5. Notebooks

```bash
conda env create -f envs/notebook.yaml
python -m ipykernel install --user --name abcfold-ago1-mirs-notebook
```

`notebooks/<complex>_domain_analysis.ipynb` — one per complex, structure
ported from the DRB2 pipeline's `rna_ds_drb2_drb4_domain_analysis.ipynb`.
**Domain boundaries in these notebooks are best-effort / approximate**
(AGO1 N-ext/N/PAZ/MID/PIWI, FBW2 F-box/WD40, coarse CUL1 split) — not a
curated PROSITE annotation like the DRB2 project's. Each notebook has a
`USE_WINDOWS = True` toggle that falls back to fixed-width residue windows
so every heatmap renders regardless; fill in `CURATED_DOMAINS` with
verified ranges when you have them.

---

## Relationship to `esther_project`

`esther_project` is to this repo what
`ab_initio_modelling_drbs_dcl4_ds_rna_complexes` was to the DRB2 ABCfold
pipeline: the manual-AF3-webserver predecessor. Its two PLIP config
"views" map onto this repo's two PLIP passes:

| esther_project config | this repo |
|---|---|
| `configs/AGO1_FBW2_interactions.yaml` (receptor = A) | main `plip:` block — receptor = anchor (AGO1) vs. every other chain |
| `configs/mir_ligands.yaml` (miR nt as ligands) | `plip_mir_ligands:` block — `--chains` omitted, no `--dnareceptor` |

The old webserver CIFs under `esther_project/data/` are **not** reused —
this repo regenerates every structure with the six ABCfold backends.
