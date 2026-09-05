#!/usr/bin/env python3
"""Regenerate one <complex>_domain_analysis.ipynb per AGO1/MiRs complex.

Structure mirrors the DRB2 pipeline's
notebooks/rna_ds_drb2_drb4_domain_analysis.ipynb (pooled / per-cluster /
per-backend domain-contact heatmaps + residue-level interface map).

Run this after editing CURATED_DOMAINS below (e.g. once you have
PROSITE/InterPro-verified domain boundaries — see HANDOFF.md task A):

    python notebooks/generate_domain_notebooks.py

Overwrites notebooks/<complex>_domain_analysis.ipynb in place. Only
regenerates notebooks/pose_clustering.ipynb if the DRB2 sibling repo is
present at the path below."""

import json
import pathlib

DST = pathlib.Path(__file__).resolve().parent.parent
NB_DIR = DST / "notebooks"
_SIBLING_POSE = pathlib.Path(
    "/Users/ereboul/projects/drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes/notebooks/pose_clustering.ipynb"
)
SRC_POSE = _SIBLING_POSE if _SIBLING_POSE.exists() else None


def _customize_pose_clustering(nb, complexes):
    """Post-process the DRB2-derived pose_clustering.ipynb:
      1. axis titles carry PCA explained variance (from
         results/<complex>/pose_clusters_pca_variance.json);
      2. one standalone example section per complex (the DRB2 source ships a
         single hard-coded complex).
    Assumes the DRB2 layout: cell 0 intro md, cell 1 helpers code, cell 2
    "## Load data" md, then a load cell + "## Examples" md + plot cells."""
    import secrets

    def _c(kind, text):
        base = {"cell_type": kind, "id": secrets.token_hex(4), "metadata": {},
                "source": text.splitlines(keepends=True)}
        if kind == "code":
            base.update(execution_count=None, outputs=[])
        return base

    OLD_LAYOUT = (
        '    fig.update_layout(\n'
        '        title=title, xaxis_title="PC1", yaxis_title="PC2",\n'
        '        template="plotly_white", height=620, width=820,\n'
        '        legend=dict(font=dict(size=9)),\n'
        '    )')
    NEW_LAYOUT = (
        '    var_path = (RESULTS_ROOT / complex_name / "pose_clusters_pca_variance.json") if complex_name else None\n'
        '    pca_var = json.loads(var_path.read_text()) if var_path and var_path.exists() else {}\n'
        '    x_title = (f"PC1 ({pca_var[\'pc1_explained_variance\']:.1%} explained variance)"\n'
        '               if pca_var else "PC1")\n'
        '    y_title = (f"PC2 ({pca_var[\'pc2_explained_variance\']:.1%} explained variance)"\n'
        '               if pca_var else "PC2")\n'
        '\n'
        '    fig.update_layout(\n'
        '        title=title, xaxis_title=x_title, yaxis_title=y_title,\n'
        '        template="plotly_white", height=620, width=820,\n'
        '        legend=dict(font=dict(size=9)),\n'
        '    )')
    patched = 0
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        if OLD_LAYOUT in src:
            cell["source"] = src.replace(OLD_LAYOUT, NEW_LAYOUT).splitlines(keepends=True)
            patched += 1
    if patched != 1:
        raise SystemExit(f"pose_clustering axis-variance patch: expected 1 site, hit {patched} "
                         "(DRB2 source layout changed -- update _customize_pose_clustering)")

    cells = nb["cells"][:2]
    cells.append(_c("markdown", "## Load data\n"))
    cells.append(_c("code",
        "COMPLEXES = [\n" + "".join(f'    "{c}",\n' for c in complexes) + "]\n"
        "\n"
        "# results/<complex>/pose_clusters.csv exists for all complexes regardless of\n"
        "# postprocessing progress (compress -> pose_cluster -> select runs early).\n"
        "# load_pose_clusters prints an n-models + PC1/PC2 explained-variance line each.\n"
        "dfs = {name: load_pose_clusters(name) for name in COMPLEXES}\n"))
    cells.append(_c("markdown",
        "## Examples\n\n"
        "One `plot_pca` call per cell: each opens its own figure, so you can\n"
        "zoom/pan/hover one plot without a later cell's output replacing it. Every\n"
        "section is standalone -- it rebinds its own `df_<complex>` from `dfs`.\n"))
    for c in complexes:
        cells.append(_c("markdown", f"### `{c}`\n"))
        cells.append(_c("code", f'df_{c} = dfs["{c}"]\ndf_{c}.head()\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}")  # pipeline\'s own hierarchical-RMSD clusters\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}", color_by="backend", cluster_method=None)  # coloured by backend\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}", cluster_method="gmm", n_components="auto")  # GMM auto (BIC-knee)\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}", cluster_method="gmm", n_components=4)  # GMM manual k -- adjust per complex\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}", cluster_method="hdbscan", n_components="auto")  # HDBSCAN auto (Optuna/DBCV)\n'))
        cells.append(_c("code",
            f'plot_pca(df_{c}, complex_name="{c}", models={{**{{b: True for b in df_{c}["backend"].unique()}}, "alphafold3": False}})  # ablation: AF3 excluded\n'))
    nb["cells"] = cells
    return nb

# ── per-complex parameters ───────────────────────────────────────────────────
# chain letter -> (display name, kind)   kind in {"protein", "rna"}
COMPLEXES = {
    "ago1_fbw2": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein")},
        "plip_chains": "[['A'], ['B']]",
        "fix_pdb": False,
        "has_mir": False,
    },
    "ago1_fbw2_mir165a": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein"), "C": ("miR165a", "rna")},
        "plip_chains": "[['A'], ['B', 'C']]",
        "fix_pdb": True,
        "has_mir": True,
    },
    "ago1_fbw2_mir168": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein"), "C": ("miR168", "rna")},
        "plip_chains": "[['A'], ['B', 'C']]",
        "fix_pdb": True,
        "has_mir": True,
    },
    "ago1_fbw2_mir393a": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein"), "C": ("miR393a", "rna")},
        "plip_chains": "[['A'], ['B', 'C']]",
        "fix_pdb": True,
        "has_mir": True,
    },
    "ago1_fbw2_ask1_cul1": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein"),
                   "C": ("ASK1", "protein"), "D": ("CUL1", "protein")},
        "plip_chains": "[['A'], ['B', 'C', 'D']]",
        "fix_pdb": False,
        "has_mir": False,
    },
    "ago1_fbw2_ask1_cul1_mir168": {
        "chains": {"A": ("AGO1", "protein"), "B": ("FBW2", "protein"),
                   "C": ("ASK1", "protein"), "D": ("CUL1", "protein"), "E": ("miR168", "rna")},
        "plip_chains": "[['A'], ['B', 'C', 'D', 'E']]",
        "fix_pdb": True,
        "has_mir": True,
    },
}

# chain lengths (from configs/*.yaml sequences) — used for the auto chain-offset
# search. All four confirmed 1:1 against the canonical UniProt entry
# (AGO1 O04379 folds 1050 vs 1048 canonical — a +2 N-terminal shift, so the
# coordinates below, taken from a scan of the *folded* sequence, are already
# construct-native; FBW2 Q9ZPE4 / ASK1 Q39255 / CUL1 Q94AH6 match exactly).
CHAIN_LENGTHS = {"AGO1": 1050, "FBW2": 317, "ASK1": 160, "CUL1": 738,
                 "miR165a": 21, "miR168": 21, "miR393a": 22}

# Domain tables — VERIFIED 2026-08-28 against the exact folded sequences via
# ScanProsite (REST, PSScan.cgi) + InterProScan 5 (EBI REST: Pfam / SMART /
# PROSITE profiles+patterns / CDD / SUPERFAMILY / Gene3D), cross-checked against
# UniProt "Family & Domains". Boundaries are 1-based inclusive, contiguous, and
# non-overlapping (linkers folded into the nearest domain). PROSITE profile
# calls are quoted verbatim where one exists (AGO1 PAZ/PIWI, CUL1 cullin
# homology); otherwise the consolidated InterPro/Pfam call is used. The notebook
# still has a fixed-width window fallback (USE_WINDOWS) so heatmaps render
# regardless.
#
# Provenance per domain (accession @ scanned span):
#   AGO1  N_ext_Grich  UniProt disorder 1-143/165-188 + Pfam PF12764 Gly-rich 75-172
#         ArgoN        Pfam PF16486 190-325
#         ArgoL1       Pfam PF08699 / SMART SM01163 336-388
#         PAZ          PROSITE profile PS50821 390-503 (Pfam PF02170 411-519)
#         ArgoL2       Pfam PF16488 529-575
#         MID          Pfam PF16487 ArgoMid 586-661
#         PIWI         PROSITE profile PS50822 678-999 (Pfam PF02171 679-998) + C-term tail
#   FBW2  F_box        UniProt 7-54 / Pfam PF12937 F-box-like 17-56
#         LRR_solenoid SUPERFAMILY SSF52047 RNI-like 36-227 / Gene3D LRR 54-244
#                      (NB: no curated WD40 despite the gene name — the C-lobe
#                       scans as an LRR / ribonuclease-inhibitor-like solenoid)
#         C_tail       245-317, acidic/aromatic-rich extension, no domain call
#   ASK1  SKP1_POZ     Pfam PF03931 Skp1_POZ 4-63 / SUPERFAMILY POZ 5-64
#         SKP1_dimer   Pfam PF01466 111-158 / SUPERFAMILY 83-158
#                      (UniProt "interaction with F-box" 102-160)
#   CUL1  cullin_repeats  Pfam PF00888 31-483 / SUPERFAMILY SSF74788 8-378 (N-term α-solenoid)
#         cullin_homology PROSITE profile PS50069 383-613 / SUPERFAMILY 382-651
#         cullin_CTD      Pfam PF10557 Cullin_Nedd8 668-730 + PROSITE pattern
#                         PS01256 711-738 + Gene3D winged-helix 654-738
CURATED_DOMAINS = {
    "AGO1": [
        ("N_ext_Grich", 1, 189),     # disordered Gly/Gln-rich N-terminal extension
        ("ArgoN", 190, 335),
        ("ArgoL1", 336, 389),
        ("PAZ", 390, 503),
        ("ArgoL2", 504, 585),
        ("MID", 586, 677),
        ("PIWI", 678, 1050),
    ],
    "FBW2": [
        ("F_box", 1, 54),
        ("LRR_solenoid", 55, 244),
        ("C_tail", 245, 317),
    ],
    "ASK1": [
        ("SKP1_POZ", 1, 82),
        ("SKP1_dimer", 83, 160),
    ],
    "CUL1": [
        ("cullin_repeats", 1, 382),   # N-terminal cullin-repeat α-solenoid
        ("cullin_homology", 383, 651),
        ("cullin_CTD", 652, 738),     # 4-helix bundle + winged-helix / Nedd8 / RBX1-binding
    ],
}


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(keepends=True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.rstrip("\n").splitlines(keepends=True)}


def build_notebook(complex_name, spec):
    chains = spec["chains"]
    protein_partners = [c for c, (_, k) in chains.items() if k == "protein" and c != "A"]
    rna_partners = [c for c, (_, k) in chains.items() if k == "rna"]
    couples = [c for c in chains if c != "A"]
    recv_name = chains["A"][0]

    chain_names_py = "{" + ", ".join(f'"{c}": "{nm}"' for c, (nm, _) in chains.items() if c != "A") + "}"
    couples_py = "[" + ", ".join(f'"{c}"' for c in couples) + "]"
    protein_partner_names = {c: chains[c][0] for c in protein_partners}
    rna_partner_names = {c: chains[c][0] for c in rna_partners}

    domain_lines = []
    for c in protein_partners + ["A"]:
        nm = chains[c][0]
        if nm in CURATED_DOMAINS:
            rows = ",\n        ".join(f'("{lbl}", {s}, {e})' for lbl, s, e in CURATED_DOMAINS[nm])
            domain_lines.append(f'    "{nm}": [\n        {rows},\n    ]')
    curated_py = "CURATED_DOMAINS = {\n" + ",\n".join(domain_lines) + ",\n}"

    chain_len_py = "{" + ", ".join(
        f'"{c}": {CHAIN_LENGTHS[chains[c][0]]}' for c in chains) + "}"
    receptor_len = CHAIN_LENGTHS[recv_name]

    mir_note = ""
    if spec["has_mir"]:
        mir_note = (
            "\nThis complex has a miR RNA chain. The **main** PLIP pass "
            f"(`plip.chains = {spec['plip_chains']}`, `dnareceptor: true`) folds the "
            "miR into the *receptor* group alongside AGO1, so it only surfaces "
            "AGO1-vs-protein-partner contacts. The per-nucleotide miR view comes "
            "from the separate `plip_mir_ligands` pass — load "
            "`all_selected_summary_mir_ligands.csv` (section at the end) for that.\n"
        )

    cells = []

    cells.append(md(
        f"# {recv_name} ("
        + ") / ".join(f"{nm} ({c})" for c, (nm, _) in chains.items())
        + ") Domain Contact Analysis\n\n"
        "**Kernel:** `abcfold-ago1-mirs-notebook` (`envs/notebook.yaml`) — install once:\n"
        "```\n"
        "conda env create -f envs/notebook.yaml\n"
        "conda activate abcfold-ago1-mirs-notebook\n"
        "python -m ipykernel install --user --name abcfold-ago1-mirs-notebook\n"
        "```\n\n"
        f"PLIP contacts from `results/{complex_name}/all_selected_summary.csv`, produced by\n"
        "`worflows/postprocessing/Snakefile` (stage 3f `run_plip` + 3h `aggregate`).\n"
        "Same structure as the sibling DRB2 pipeline's "
        "`notebooks/rna_ds_drb2_drb4_domain_analysis.ipynb`:\n"
        "pooled / per-pose-cluster / per-backend domain-contact heatmaps + a\n"
        "residue-level interface map.\n"
        f"{mir_note}\n"
        "**Receptor** is always the pose-cluster anchor, AGO1 (chain A). Every\n"
        "heatmap is one *couple*: AGO1 vs. one partner chain.\n\n"
        "---\n\n"
        "**Domain boundaries below are PROSITE / InterPro-verified** (ScanProsite +\n"
        "InterProScan 5, scanned against the exact folded sequences on 2026-08-28,\n"
        "cross-checked against UniProt — see the provenance block in\n"
        "`notebooks/generate_domain_notebooks.py`). PROSITE profile spans are quoted\n"
        "verbatim where one exists (AGO1 PAZ/PIWI, CUL1 cullin homology). Set\n"
        "`USE_WINDOWS = True` in the domain-definitions cell to fall back to plain\n"
        "fixed-width residue windows instead (every heatmap still renders, axis\n"
        "labels just become `1-60`, `61-120`, …).\n"
    ))

    cells.append(code(
        "import pandas as pd\n"
        "import numpy as np\n"
        "import plotly.express as px\n"
        "import plotly.graph_objects as go\n"
        "from plotly.subplots import make_subplots\n"
        "from pathlib import Path\n\n"
        "ROOT = Path(\"..\")\n"
        f"RESULTS_DIR = ROOT / \"results\" / \"{complex_name}\"\n"
        f"RECEPTOR_CHAIN, RECEPTOR_NAME = \"A\", \"{recv_name}\"\n"
        f"LIGAND_CHAIN_NAMES = {chain_names_py}\n"
        f"PROTEIN_PARTNERS = {json.dumps(protein_partner_names)}\n"
        f"RNA_PARTNERS = {json.dumps(rna_partner_names)}\n"
        f"CHAIN_LENGTHS = {chain_len_py}\n"
        f"RECEPTOR_LEN = {receptor_len}\n"
        "FIGURES_DIR = RESULTS_DIR / \"figures\" / \"domain_analysis\"\n\n"
        "TEMPLATE = \"plotly_white\"\n"
        "BACKEND_PALETTE = px.colors.qualitative.Set2\n"
        "ITYPE_PALETTE = px.colors.qualitative.Set1\n\n"
        "def out_path(subdir, filename):\n"
        "    out_dir = FIGURES_DIR / subdir\n"
        "    out_dir.mkdir(parents=True, exist_ok=True)\n"
        "    return out_dir / filename\n\n"
        "def save_fig(fig, filename, subdir=\"\"):\n"
        "    out = out_path(subdir, filename)\n"
        "    fig.write_html(out, include_plotlyjs=\"cdn\")\n"
        "    print(f\"Saved: {out}\")\n"
        "    fig.show()\n"
    ))

    cells.append(md("## Load data"))
    cells.append(code(
        "csv_path = RESULTS_DIR / \"all_selected_summary.csv\"\n"
        "df = pd.read_csv(csv_path)\n"
        "df = df.rename(columns={\"replica\": \"cluster\", \"model\": \"fname\"})\n"
        "df[\"cluster\"] = df[\"cluster\"].astype(int)\n\n"
        "sel = pd.read_csv(RESULTS_DIR / \"selected_models.csv\")\n"
        "sel[\"fname\"] = sel[\"staged_cif\"].apply(lambda p: Path(p).stem)\n"
        "sel = sel[[\"fname\", \"cluster\", \"backend\", \"seed\", \"sample_index\", \"ranking_score\", \"ptm\", \"iptm\"]]\n\n"
        "df = df.merge(sel, on=[\"fname\", \"cluster\"], how=\"left\", validate=\"many_to_one\")\n"
        "n_missing_meta = df[\"backend\"].isna().sum()\n"
        "if n_missing_meta:\n"
        "    print(f\"WARNING: {n_missing_meta} contact rows have no matching selected_models.csv entry\")\n\n"
        "n_models = df.groupby([\"cluster\", \"fname\"]).ngroups\n"
        "print(f\"{csv_path}: {len(df)} contact rows, {df['cluster'].nunique()} pose cluster(s), \"\n"
        "      f\"{n_models} model(s), before the energy filter below\")\n"
        "print(\"\\ncontact rows per backend:\")\n"
        "print(df.groupby(\"backend\").size().rename(\"rows\").to_frame())\n"
    ))

    cells.append(md(
        "## Filter out numerically-unconverged structures\n\n"
        "Robust (median / MAD) modified z-score on each minimised structure's final\n"
        "ChimeraX energy — same approach as the DRB2 pipeline's notebooks.\n\n"
        "1. **Whole-backend cut:** if a backend has >50 % of its energy-assessed models\n"
        "   flagged as outliers, the *entire* backend is dropped.\n"
        "2. **Missing `_energy.csv` ⇒ keep:** minimisations that produced a `*_fixed.pdb`\n"
        "   but no parseable energy trajectory are kept (\"not assessed\", not \"bad\") for\n"
        "   a non-dropped backend.\n"
    ))
    cells.append(code(
        "def read_final_energy(energy_csv_path):\n"
        "    if not energy_csv_path.exists():\n"
        "        return np.nan\n"
        "    last_energy = np.nan\n"
        "    with open(energy_csv_path) as fh:\n"
        "        next(fh, None)\n"
        "        for line in fh:\n"
        "            parts = line.strip().split(\",\")\n"
        "            if len(parts) < 2:\n"
        "                continue\n"
        "            try:\n"
        "                last_energy = float(parts[1])\n"
        "            except ValueError:\n"
        "                pass\n"
        "    return last_energy\n\n"
        "energy_rows = []\n"
        "for pdb_path in sorted(RESULTS_DIR.glob(\"minimized/*/*/*.pdb\")):\n"
        "    if pdb_path.stem.endswith((\"_fixed\", \"_amber\", \"_nonprot\")):\n"
        "        continue\n"
        "    cluster = int(pdb_path.parent.parent.name)\n"
        "    energy_csv = pdb_path.with_name(pdb_path.stem + \"_energy.csv\")\n"
        "    energy_rows.append({\"fname\": pdb_path.stem, \"cluster\": cluster,\n"
        "                        \"final_energy\": read_final_energy(energy_csv)})\n\n"
        "energy_all = pd.DataFrame(energy_rows).merge(\n"
        "    sel[[\"fname\", \"cluster\", \"backend\"]], on=[\"fname\", \"cluster\"], how=\"left\")\n"
        "energy_df = energy_all.dropna(subset=[\"final_energy\"]).copy()\n\n"
        "MOD_Z_THRESHOLD = 3.5\n"
        "pooled_median = energy_df[\"final_energy\"].median()\n"
        "pooled_mad = (energy_df[\"final_energy\"] - pooled_median).abs().median()\n"
        "energy_df[\"energy_mod_z\"] = 0.6745 * (energy_df[\"final_energy\"] - pooled_median) / pooled_mad\n"
        "energy_df[\"energy_ok\"] = energy_df[\"energy_mod_z\"].abs() <= MOD_Z_THRESHOLD\n\n"
        "print(f\"Pooled final-energy median={pooled_median:,.0f} kJ/mol, MAD={pooled_mad:,.0f}\")\n"
        "print(energy_df.groupby(\"backend\")[\"final_energy\"].agg([\"count\", \"min\", \"median\", \"max\"]).round(1).to_string())\n\n"
        "flag_rate = (1 - energy_df.groupby(\"backend\")[\"energy_ok\"].mean()).rename(\"flagged_frac\")\n"
        "DROP_BACKENDS = sorted(flag_rate[flag_rate > 0.5].index)\n"
        "print(f\"\\nflagged fraction by backend (|mod-z| > {MOD_Z_THRESHOLD}):\")\n"
        "print(flag_rate.round(3).to_string())\n"
        "print(f\"\\n=> whole-backend drop (>50% flagged): {DROP_BACKENDS or 'none'}\")\n"
    ))
    cells.append(code(
        "assessed_ok = set(zip(energy_df.loc[energy_df[\"energy_ok\"], \"cluster\"],\n"
        "                      energy_df.loc[energy_df[\"energy_ok\"], \"fname\"]))\n"
        "unassessed = set(zip(energy_all.loc[energy_all[\"final_energy\"].isna(), \"cluster\"],\n"
        "                     energy_all.loc[energy_all[\"final_energy\"].isna(), \"fname\"]))\n"
        "keep_pairs = assessed_ok | unassessed\n\n"
        "drop_fnames = set(sel.loc[sel[\"backend\"].isin(DROP_BACKENDS), \"fname\"])\n"
        "df_keys = pd.MultiIndex.from_arrays([df[\"cluster\"], df[\"fname\"]])\n"
        "n_models_before = df.groupby([\"cluster\", \"fname\"]).ngroups\n\n"
        "df = df[df_keys.isin(keep_pairs) & ~df[\"fname\"].isin(drop_fnames)].copy()\n\n"
        "n_models_after = df.groupby([\"cluster\", \"fname\"]).ngroups\n"
        "print(f\"{n_models_after} / {n_models_before} model(s) kept after the energy filter\")\n"
        "print(df.groupby(\"backend\")[\"fname\"].nunique().rename(\"n_models\").to_frame())\n"
    ))

    cells.append(md(
        "## Domain definitions\n\n"
        "`CURATED_DOMAINS` below is the PROSITE / InterPro-verified table (see the\n"
        "header note above). `USE_WINDOWS = True` ignores it and bins every chain\n"
        "into fixed-width residue windows instead."
    ))
    cells.append(code(
        "USE_WINDOWS = False        # True -> fixed-width residue windows instead of CURATED_DOMAINS\n"
        "WINDOW = 60               # residues per window when USE_WINDOWS\n\n"
        "# PROSITE / InterPro-verified boundaries (1-based inclusive), scanned\n"
        "# against the exact folded sequences 2026-08-28. See generate_domain_notebooks.py.\n"
        f"{curated_py}\n\n"
        "def windows_for(chain_len, w=None):\n"
        "    w = w or WINDOW\n"
        "    return [(f\"{s}-{min(s+w-1, chain_len)}\", s, min(s + w - 1, chain_len))\n"
        "            for s in range(1, chain_len + 1, w)]\n\n"
        "def domains_for(name, chain_len):\n"
        "    if USE_WINDOWS or name not in CURATED_DOMAINS:\n"
        "        return windows_for(chain_len)\n"
        "    return CURATED_DOMAINS[name]\n\n"
        "def make_domain_mapper(domain_ranges):\n"
        "    intervals = pd.IntervalIndex.from_tuples(\n"
        "        [(s, e) for _, s, e in domain_ranges], closed=\"both\")\n"
        "    labels = [l for l, _, _ in domain_ranges]\n"
        "    def mapper(resnr_series):\n"
        "        idx = intervals.get_indexer(resnr_series.astype(float))\n"
        "        return pd.Series([labels[i] if i != -1 else pd.NA for i in idx],\n"
        "                         index=resnr_series.index, dtype=\"object\")\n"
        "    return mapper\n\n"
        "RECEPTOR_DOMAINS = domains_for(RECEPTOR_NAME, RECEPTOR_LEN)\n"
        "RECEPTOR_LABELS = [l for l, _, _ in RECEPTOR_DOMAINS]\n"
        "receptor_mapper = make_domain_mapper(RECEPTOR_DOMAINS)\n\n"
        "LIGAND_DOMAINS = {c: domains_for(nm, CHAIN_LENGTHS[c]) for c, nm in PROTEIN_PARTNERS.items()}\n"
        "LIGAND_LABELS = {c: [l for l, _, _ in d] for c, d in LIGAND_DOMAINS.items()}\n"
        "LIGAND_MAPPERS = {c: make_domain_mapper(d) for c, d in LIGAND_DOMAINS.items()}\n"
    ))

    cells.append(md(
        "## Chain-residue offset detection\n\n"
        "After `fix_pdb` (pdb4amber + PDBFixer) the whole complex is usually renumbered\n"
        "continuously across chains rather than restarting each chain at 1. This cell\n"
        "auto-detects the per-chain offset: for each ligand chain it finds the offset\n"
        "that lands every `resnr_lig` inside `[1, chain_length]`, so domain mapping is\n"
        "on that chain's own numbering. (If the receptor `resnr` isn't already 1-based\n"
        "it's offset too.)"
    ))
    cells.append(code(
        "def best_offset(values, chain_len):\n"
        "    v = values.dropna().astype(int)\n"
        "    if v.empty:\n"
        "        return 0\n"
        "    cand = int(v.min()) - 1\n"
        "    for off in (0, cand):\n"
        "        if v.sub(off).between(1, chain_len).all():\n"
        "            return off\n"
        "    return cand\n\n"
        "rec_off = best_offset(df[\"resnr\"], RECEPTOR_LEN)\n"
        "df[\"resnr_raw\"] = df[\"resnr\"]\n"
        "df[\"resnr\"] = df[\"resnr\"] - rec_off\n"
        "df[\"receptor_domain\"] = receptor_mapper(df[\"resnr\"])\n\n"
        "df[\"resnr_lig_raw\"] = df[\"resnr_lig\"]\n"
        "CHAIN_OFFSET = {}\n"
        "for chain, clen in CHAIN_LENGTHS.items():\n"
        "    if chain == RECEPTOR_CHAIN:\n"
        "        continue\n"
        "    mask = df[\"reschain_lig\"] == chain\n"
        "    if not mask.any():\n"
        "        continue\n"
        "    off = best_offset(df.loc[mask, \"resnr_lig\"], clen)\n"
        "    CHAIN_OFFSET[chain] = off\n"
        "    df.loc[mask, \"resnr_lig\"] = df.loc[mask, \"resnr_lig\"] - off\n"
        "    bad = df.loc[mask & ~df[\"resnr_lig\"].between(1, clen)]\n"
        "    tag = \"OK\" if bad.empty else f\"WARNING {len(bad)} rows outside [1,{clen}]\"\n"
        "    print(f\"chain {chain} ({LIGAND_CHAIN_NAMES.get(chain, '?')}): offset {off} -> {tag}\")\n\n"
        "df[\"ligand_domain\"] = pd.NA\n"
        "for chain, mapper in LIGAND_MAPPERS.items():\n"
        "    m = df[\"reschain_lig\"] == chain\n"
        "    df.loc[m, \"ligand_domain\"] = mapper(df.loc[m, \"resnr_lig\"])\n\n"
        "n_unmapped = df[\"receptor_domain\"].isna().sum()\n"
        "print(f\"\\nUnmapped {RECEPTOR_NAME} residues: {n_unmapped} ({100*n_unmapped/max(len(df),1):.1f}%)\")\n"
    ))

    cells.append(md(
        "## Interactor couples in this dataset\n\n"
        "One heatmap per couple: AGO1 (fixed receptor) vs. each partner. Couples with\n"
        "zero rows are expected for chains folded into the PLIP receptor group — the\n"
        "check below makes that explicit."
    ))
    cells.append(code(
        f"COUPLES = {couples_py}\n\n"
        "print(\"Contact rows per couple (AGO1 vs. partner):\")\n"
        "for c in COUPLES:\n"
        "    n = (df[\"reschain_lig\"] == c).sum()\n"
        "    n_models_c = df.loc[df[\"reschain_lig\"] == c, [\"cluster\", \"fname\"]].drop_duplicates().shape[0]\n"
        "    status = \"POPULATED\" if n else \"EMPTY -- folded into the PLIP receptor group or no contact\"\n"
        "    print(f\"  {RECEPTOR_NAME} x {LIGAND_CHAIN_NAMES[c]:8s}: {n:6d} rows / {n_models_c:3d} model(s) -- {status}\")\n"
    ))

    cells.append(code(
        "def domain_pair_matrix(data, ligand_chain, n_models_norm):\n"
        "    \"\"\"(x_labels, y_labels, rate matrix) for a protein couple; None if empty.\"\"\"\n"
        "    labels = LIGAND_LABELS[ligand_chain]\n"
        "    sub = data[data[\"reschain_lig\"] == ligand_chain].dropna(subset=[\"receptor_domain\", \"ligand_domain\"])\n"
        "    if sub.empty:\n"
        "        return None\n"
        "    ct = (sub.groupby([\"receptor_domain\", \"ligand_domain\"], observed=True).size()\n"
        "          .unstack(fill_value=0).reindex(index=RECEPTOR_LABELS, columns=labels, fill_value=0))\n"
        "    return labels, RECEPTOR_LABELS, (ct / n_models_norm if n_models_norm else ct)\n\n"
        "def rna_pos_matrix(data, ligand_chain, n_models_norm):\n"
        "    sub = data[data[\"reschain_lig\"] == ligand_chain].dropna(subset=[\"receptor_domain\"])\n"
        "    if sub.empty:\n"
        "        return None\n"
        "    nt = sorted(sub[\"resnr_lig\"].dropna().unique())\n"
        "    ct = (sub.groupby([\"resnr_lig\", \"receptor_domain\"], observed=True).size()\n"
        "          .unstack(fill_value=0).reindex(index=nt, columns=RECEPTOR_LABELS, fill_value=0))\n"
        "    return RECEPTOR_LABELS, [str(int(p)) for p in nt], (ct / n_models_norm if n_models_norm else ct)\n\n"
        "def couple_matrix(data, ligand_chain, n_models_norm):\n"
        "    if ligand_chain in RNA_PARTNERS:\n"
        "        return rna_pos_matrix(data, ligand_chain, n_models_norm)\n"
        "    return domain_pair_matrix(data, ligand_chain, n_models_norm)\n\n"
        "def heatmap_from_matrix(res, ligand_chain, title, filename, subdir=\"domain_contacts\"):\n"
        "    if res is None:\n"
        "        print(f\"No contacts -- skipping '{title}'.\")\n"
        "        return\n"
        "    x_labels, y_labels, rate = res\n"
        "    ligand_name = LIGAND_CHAIN_NAMES[ligand_chain]\n"
        "    is_rna = ligand_chain in RNA_PARTNERS\n"
        "    fig = go.Figure(go.Heatmap(\n"
        "        z=rate.values, x=x_labels, y=y_labels, colorscale=\"YlOrRd\",\n"
        "        text=[[f\"{v:.2f}\" if v > 0 else \"\" for v in row] for row in rate.values],\n"
        "        texttemplate=\"%{text}\", textfont=dict(size=8),\n"
        "        hovertemplate=(f\"{RECEPTOR_NAME}: %{{{'x' if is_rna else 'y'}}}<br>{ligand_name}: \"\n"
        "                       f\"%{{{'y' if is_rna else 'x'}}}<br>%{{z:.3f}} contacts/model<extra></extra>\"),\n"
        "        colorbar=dict(title=\"Mean<br>contacts/<br>model\")))\n"
        "    fig.update_layout(\n"
        "        title=title,\n"
        "        xaxis_title=f\"{RECEPTOR_NAME} domain\" if is_rna else f\"{ligand_name} domain\",\n"
        "        yaxis_title=f\"{ligand_name} nt position\" if is_rna else f\"{RECEPTOR_NAME} domain\",\n"
        "        yaxis=dict(autorange=\"reversed\"), template=TEMPLATE,\n"
        "        width=max(560, len(x_labels) * 90), height=max(420, len(y_labels) * (16 if is_rna else 64)))\n"
        "    save_fig(fig, filename, subdir)\n"
    ))

    cells.append(md("## One heatmap per couple (all backends, all clusters pooled)"))
    cells.append(code(
        "n_models_total = df.groupby([\"cluster\", \"fname\"]).ngroups\n"
        "for c in COUPLES:\n"
        "    heatmap_from_matrix(\n"
        "        couple_matrix(df, c, n_models_total), c,\n"
        "        f\"{RECEPTOR_NAME} x {LIGAND_CHAIN_NAMES[c]} -- pooled (n={n_models_total} models)\",\n"
        "        f\"{RECEPTOR_NAME.lower()}_{LIGAND_CHAIN_NAMES[c].lower()}_heatmap_pooled.html\")\n"
    ))

    cells.append(md(
        "## Per-pose-cluster breakdown\n\n"
        "One small-multiple panel per couple: one heatmap per pose cluster "
        "(`scripts/pose_cluster_anchor.py`'s `cluster` column), shared colour scale\n"
        "within the panel. If there's only one cluster this collapses to a single tile."
    ))
    cells.append(code(
        "clusters = sorted(df[\"cluster\"].unique())\n"
        "cluster_n_models = df.groupby(\"cluster\")[\"fname\"].nunique()\n"
        "print(f\"{len(clusters)} pose cluster(s): \" + \", \".join(f'{k} (n={cluster_n_models[k]})' for k in clusters))\n\n"
        "def panel_by(group_col, group_values, group_n_models, tag, subdir):\n"
        "    for c in COUPLES:\n"
        "        ligand_name = LIGAND_CHAIN_NAMES[c]\n"
        "        mats = {}\n"
        "        for g in group_values:\n"
        "            res = couple_matrix(df[df[group_col] == g], c, group_n_models[g])\n"
        "            if res is not None:\n"
        "                mats[g] = res\n"
        "        if not mats:\n"
        "            print(f\"No {RECEPTOR_NAME}-{ligand_name} contacts for any {tag} -- skipping.\")\n"
        "            continue\n"
        "        present = list(mats)\n"
        "        zmax = max(r[2].values.max() for r in mats.values()) or 1\n"
        "        ncols = min(3, len(present)); nrows = -(-len(present) // ncols)\n"
        "        fig = make_subplots(rows=nrows, cols=ncols,\n"
        "                            subplot_titles=[f\"{tag} {g} (n={group_n_models[g]})\" for g in present])\n"
        "        for i, g in enumerate(present):\n"
        "            x_labels, y_labels, rate = mats[g]\n"
        "            r, cc = i // ncols + 1, i % ncols + 1\n"
        "            fig.add_trace(go.Heatmap(z=rate.values, x=x_labels, y=y_labels, colorscale=\"YlOrRd\",\n"
        "                zmin=0, zmax=float(zmax), showscale=bool(g == present[-1]),\n"
        "                text=[[f\"{v:.2f}\" if v > 0 else \"\" for v in row] for row in rate.values],\n"
        "                texttemplate=\"%{text}\", textfont=dict(size=7),\n"
        "                hovertemplate=f\"{tag} {g}<br>%{{y}} / %{{x}}<br>%{{z:.3f}}/model<extra></extra>\"),\n"
        "                row=r, col=cc)\n"
        "        fig.update_yaxes(autorange=\"reversed\")\n"
        "        fig.update_layout(title=f\"{RECEPTOR_NAME} x {ligand_name} by {tag}\", template=TEMPLATE,\n"
        "                          width=max(760, 360 * ncols), height=max(430, 150 * nrows))\n"
        "        save_fig(fig, f\"{RECEPTOR_NAME.lower()}_{ligand_name.lower()}_heatmap_by_{tag}.html\", subdir)\n\n"
        "panel_by(\"cluster\", clusters, cluster_n_models, \"cluster\", \"per_cluster\")\n"
    ))

    cells.append(md(
        "## Per-backend breakdown\n\n"
        "Cross-architecture agreement check: do the six ABCfold backends place the\n"
        "same domain–domain contacts? One panel per couple, one heatmap per surviving\n"
        "backend, shared colour scale."
    ))
    cells.append(code(
        "backends = sorted(df[\"backend\"].dropna().unique())\n"
        "backend_n_models = df.groupby(\"backend\")[\"fname\"].nunique()\n"
        "print(f\"{len(backends)} backend(s): \" + \", \".join(f'{b} (n={backend_n_models[b]})' for b in backends))\n"
        "panel_by(\"backend\", backends, backend_n_models, \"backend\", \"per_backend\")\n"
    ))

    cells.append(md("## Interaction-type breakdown"))
    cells.append(code(
        "print(\"Interaction type counts:\")\n"
        "display(df[\"interaction_type\"].value_counts().rename(\"count\").to_frame())\n\n"
        "print(\"\\nInteraction type by receptor-ligand chain pair:\")\n"
        "display(df.groupby([\"reschain\", \"reschain_lig\", \"interaction_type\"], observed=True)\n"
        "        .size().rename(\"count\").reset_index().sort_values(\"count\", ascending=False).head(20))\n\n"
        "per_backend_itype = (df.groupby([\"backend\", \"interaction_type\"], observed=True).size()\n"
        "    .div(df.groupby(\"backend\")[\"fname\"].nunique(), level=\"backend\")\n"
        "    .rename(\"contacts_per_model\").reset_index())\n"
        "fig = px.bar(per_backend_itype, x=\"backend\", y=\"contacts_per_model\", color=\"interaction_type\",\n"
        "             color_discrete_sequence=ITYPE_PALETTE, template=TEMPLATE,\n"
        "             title=f\"{RECEPTOR_NAME} interface contacts per model, by interaction type and backend\")\n"
        "fig.update_layout(width=780, height=460, xaxis_title=\"\", yaxis_title=\"contacts / model\")\n"
        "save_fig(fig, \"interaction_types_by_backend.html\")\n"
    ))

    cells.append(md(
        "## Residue-level interface map\n\n"
        "Single-residue resolution for the busiest AGO1 and partner residues, split by\n"
        "interaction type, plus the 2-D AGO1-residue × partner-residue contact map —\n"
        "the predicted interface footprint, pooled across surviving backends. Runs for\n"
        "the first populated protein couple; edit `COUPLE_FOR_RESIDUE_MAP` to switch."
    ))
    first_prot = protein_partners[0] if protein_partners else "None"
    cells.append(code(
        f"COUPLE_FOR_RESIDUE_MAP = {first_prot!r}   # a protein partner chain letter\n"
        "TOP_N = 25\n"
        "bc = df[df[\"reschain_lig\"] == COUPLE_FOR_RESIDUE_MAP].copy()\n"
        "n_norm = df.groupby([\"cluster\", \"fname\"]).ngroups\n"
        "partner_name = LIGAND_CHAIN_NAMES.get(COUPLE_FOR_RESIDUE_MAP, \"partner\")\n\n"
        "def residue_rate_bar(data, resnr_col, restype_col, chain_label, filename):\n"
        "    if data.empty:\n"
        "        print(f\"{chain_label}: no rows -- skipping.\"); return pd.DataFrame({resnr_col: [], restype_col: []})\n"
        "    per_res_total = (data.groupby([resnr_col, restype_col], observed=True).size()\n"
        "                     .div(n_norm).rename(\"rate\").reset_index()\n"
        "                     .sort_values(\"rate\", ascending=False).head(TOP_N))\n"
        "    order = per_res_total[resnr_col].tolist()\n"
        "    per_res_itype = (data.groupby([resnr_col, restype_col, \"interaction_type\"], observed=True).size()\n"
        "                     .div(n_norm).rename(\"rate\").reset_index())\n"
        "    per_res_itype = per_res_itype[per_res_itype[resnr_col].isin(order)].copy()\n"
        "    per_res_itype[\"label\"] = (per_res_itype[restype_col].astype(str)\n"
        "                              + per_res_itype[resnr_col].astype(int).astype(str))\n"
        "    label_order = [f\"{per_res_total.loc[per_res_total[resnr_col] == r, restype_col].iloc[0]}{int(r)}\"\n"
        "                   for r in order]\n"
        "    fig = px.bar(per_res_itype, x=\"label\", y=\"rate\", color=\"interaction_type\",\n"
        "                 category_orders={\"label\": label_order}, color_discrete_sequence=ITYPE_PALETTE,\n"
        "                 template=TEMPLATE, title=f\"{chain_label} interface residues -- top {TOP_N} by contact rate (n={n_norm})\")\n"
        "    fig.update_layout(width=950, height=460, xaxis_title=f\"{chain_label} residue\",\n"
        "                      yaxis_title=\"contacts / model\", xaxis_tickangle=-45)\n"
        "    save_fig(fig, filename, \"residue_interface\")\n"
        "    return per_res_total\n\n"
        "top_rec = residue_rate_bar(bc, \"resnr\", \"restype\", RECEPTOR_NAME, f\"{RECEPTOR_NAME.lower()}_top_residues.html\")\n"
        "top_lig = residue_rate_bar(bc, \"resnr_lig\", \"restype_lig\", partner_name, f\"{partner_name.lower()}_top_residues.html\")\n"
        "display(top_rec.reset_index(drop=True)); display(top_lig.reset_index(drop=True))\n"
    ))
    cells.append(code(
        "if not bc.empty and len(top_rec) and len(top_lig):\n"
        "    tr = top_rec[\"resnr\"].tolist(); tl = top_lig[\"resnr_lig\"].tolist()\n"
        "    pair = bc[bc[\"resnr\"].isin(tr) & bc[\"resnr_lig\"].isin(tl)]\n"
        "    mat = (pair.groupby([\"resnr\", \"resnr_lig\"], observed=True).size()\n"
        "           .div(n_norm).unstack(fill_value=0)\n"
        "           .reindex(index=sorted(tr), columns=sorted(tl), fill_value=0))\n"
        "    rtick = {int(r): f\"{bc.loc[bc['resnr']==r,'restype'].iloc[0]}{int(r)}\" for r in mat.index}\n"
        "    ltick = {int(c): f\"{bc.loc[bc['resnr_lig']==c,'restype_lig'].iloc[0]}{int(c)}\" for c in mat.columns}\n"
        "    fig = go.Figure(go.Heatmap(z=mat.values, x=[ltick[c] for c in mat.columns], y=[rtick[r] for r in mat.index],\n"
        "        colorscale=\"YlOrRd\", hovertemplate=f\"{RECEPTOR_NAME} %{{y}}<br>{partner_name} %{{x}}<br>%{{z:.3f}}/model<extra></extra>\",\n"
        "        colorbar=dict(title=\"contacts<br>/ model\")))\n"
        "    fig.update_layout(title=f\"{RECEPTOR_NAME} x {partner_name} residue-residue contact map -- top {TOP_N} each (n={n_norm})\",\n"
        "        xaxis_title=f\"{partner_name} residue\", yaxis_title=f\"{RECEPTOR_NAME} residue\",\n"
        "        yaxis=dict(autorange=\"reversed\"), template=TEMPLATE, width=900, height=760)\n"
        "    save_fig(fig, f\"{RECEPTOR_NAME.lower()}_{partner_name.lower()}_residue_contact_map.html\", \"residue_interface\")\n"
    ))

    if spec["has_mir"]:
        cells.append(md(
            "## miR-ligand PLIP pass (per-nucleotide)\n\n"
            f"`results/{complex_name}/all_selected_summary_mir_ligands.csv` — the second\n"
            "PLIP pass with `--chains` omitted and no `--dnareceptor`, so every miR\n"
            "nucleotide is reported as its own ligand. Here we filter to AGO1 (chain A)\n"
            "contacts and map each nucleotide position × AGO1 domain."
        ))
        cells.append(code(
            "mir_csv = RESULTS_DIR / \"all_selected_summary_mir_ligands.csv\"\n"
            "if mir_csv.exists():\n"
            "    mdf = pd.read_csv(mir_csv).rename(columns={\"replica\": \"cluster\", \"model\": \"fname\"})\n"
            "    mdf = mdf.merge(sel, on=[\"fname\"], how=\"left\", suffixes=(\"\", \"_sel\"))\n"
            "    # keep rows where AGO1 (chain A) is on either side of the contact\n"
            "    a_side = mdf[mdf[\"reschain\"] == \"A\"].copy()\n"
            "    a_side[\"ago1_domain\"] = receptor_mapper(a_side[\"resnr\"] - rec_off)\n"
            "    nt = sorted(a_side[\"resnr_lig\"].dropna().unique())\n"
            "    n_norm_m = a_side.groupby([\"cluster\", \"fname\"]).ngroups or 1\n"
            "    ct = (a_side.dropna(subset=[\"ago1_domain\"]).groupby([\"resnr_lig\", \"ago1_domain\"], observed=True)\n"
            "          .size().div(n_norm_m).unstack(fill_value=0)\n"
            "          .reindex(index=nt, columns=RECEPTOR_LABELS, fill_value=0))\n"
            "    fig = go.Figure(go.Heatmap(z=ct.values, x=RECEPTOR_LABELS, y=[str(int(p)) for p in nt],\n"
            "        colorscale=\"YlOrRd\", colorbar=dict(title=\"contacts<br>/ model\"),\n"
            "        hovertemplate=\"AGO1 domain %{x}<br>miR nt %{y}<br>%{z:.3f}/model<extra></extra>\"))\n"
            "    fig.update_layout(title=f\"miR nucleotide x AGO1 domain contact rate (n={n_norm_m} models)\",\n"
            "        xaxis_title=\"AGO1 domain\", yaxis_title=\"miR nt position\",\n"
            "        yaxis=dict(autorange=\"reversed\"), template=TEMPLATE, width=760, height=max(400, len(nt) * 18))\n"
            "    save_fig(fig, \"mir_nt_x_ago1_domain_heatmap.html\", \"mir_ligands\")\n"
            "else:\n"
            "    print(f\"{mir_csv} not found -- run the plip_mir_ligands pass first.\")\n"
        ))

    cells.append(md("## Export contact tables"))
    cells.append(code(
        "export_dir = out_path(\"tables\", \"\")\n"
        "for c in COUPLES:\n"
        "    for tag, data in [(\"pooled\", df)] + [(b, df[df[\"backend\"] == b]) for b in backends]:\n"
        "        res = couple_matrix(data, c, data.groupby([\"cluster\", \"fname\"]).ngroups)\n"
        "        if res is None:\n"
        "            continue\n"
        "        _, _, rate = res\n"
        "        p = export_dir / f\"{RECEPTOR_NAME.lower()}_{LIGAND_CHAIN_NAMES[c].lower()}_rate_{tag}.csv\"\n"
        "        rate.to_csv(p)\n"
        "        print(f\"Saved: {p}\")\n"
    ))

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "abcfold-ago1-mirs-notebook",
                           "language": "python", "name": "abcfold-ago1-mirs-notebook"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


for name, spec in COMPLEXES.items():
    nb = build_notebook(name, spec)
    p = NB_DIR / f"{name}_domain_analysis.ipynb"
    p.write_text(json.dumps(nb, indent=1))
    print("wrote", p)

# ── pose_clustering.ipynb: copy + retarget from the DRB2 sibling (if present) ─
if SRC_POSE is None:
    print("skipped pose_clustering.ipynb (DRB2 sibling repo not found) — "
          "the existing notebooks/pose_clustering.ipynb is left untouched")
else:
    rep = SRC_POSE.read_text().replace("rna_ds_dcl4_drb2_drb4", "ago1_fbw2_mir168")
    rep = rep.replace("abcfold-drbs-notebook", "abcfold-ago1-mirs-notebook")
    rep = rep.replace(
        "../../NPF-ab-initio-modelling/ABCfold_NPF_pipeline",
        "../../drb2_modelling/ABCfold_ifb_drbs_dcl4_ds_rna_complexes (sibling)")
    pose_nb = _customize_pose_clustering(json.loads(rep), list(COMPLEXES))
    (NB_DIR / "pose_clustering.ipynb").write_text(json.dumps(pose_nb, indent=1) + "\n")
    print("wrote", NB_DIR / "pose_clustering.ipynb", f"({len(pose_nb['cells'])} cells, {len(COMPLEXES)} complexes)")
