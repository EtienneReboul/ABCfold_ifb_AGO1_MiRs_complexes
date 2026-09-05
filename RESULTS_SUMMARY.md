# AGO1 / FBW2 / miR complex modelling — results summary

*Generated 2026-09-04 from the completed ABCfold 6-backend rebuild + postprocessing
(pose clustering → top-N selection → ChimeraX minimise → PLIP → domain analysis).
Numbers below are pooled over all retained models per complex; see the executed
notebooks and `results/<complex>/figures/domain_analysis/**` for the full plots.*

---

## 1. Headline findings

1. **FBW2 grips AGO1 through its C-terminal half, not its F-box.** The FBW2
   `LRR_solenoid` + `C_tail` (res ~55–317) carry essentially the whole
   interface. The FBW2 `F_box` (res 1–54) is nearly free: only a low-level
   AGO1 `N_ext_Grich`×`F_box` contact (~2.8/model in the bare AGO1–FBW2 run)
   that **collapses to ~0.5/model as soon as ASK1/SKP1 is present** — i.e. the
   F-box is available for its real partner. The models are geometrically
   compatible with FBW2 acting as an SCF substrate receptor that presents AGO1.

2. **Two distinguishable AGO1–FBW2 contact modes:**
   - a **focused aromatic/hydrophobic patch**: AGO1 **PAZ** (Y466, Y471, F473,
     R474, R430, R438, K424/436, N487, P390) against the **FBW2 LRR→C_tail
     junction** (Y140, L169, T198?, V248, F251, Y252, D253, D260, D264). These
     are the highest-rate *individual* residue–residue pairs (0.3–0.55/model) and
     recur across backends → the most likely genuine anchor.
   - a **diffuse electrostatic sheet**: AGO1's basic **PIWI** face (R889, R907,
     R942, R973, R976, K888, R506, R222 …) against the acidic **FBW2 C_tail**
     (D280, E279, D264, E285, Y252, W275, E288 …). Very contact-rich in bulk
     (dominates the domain-level rate table) but spread over many weak pairs and
     variable geometry — partly the mutual attraction of two highly charged
     surfaces.

3. **In the SCF context (ASK1/CUL1 complexes)** the picture is internally
   consistent: **CUL1's N-terminal arm** (`cullin_repeats` + `cullin_homology`)
   packs against **AGO1's N-terminal lobe** (`N_ext_Grich` + `ArgoN`;
   N_ext×cullin_homology ≈ 18–21 contacts/model — the single largest domain pair
   in those runs), **ASK1 barely touches AGO1** (≤1.8/model, as expected for the
   SKP1 adaptor that bridges F-box↔CUL1), and the FBW2 C_tail→AGO1 PIWI grip
   persists (≈57/model). So the ensemble places AGO1 as the substrate, clamped
   between the FBW2 receptor and the cullin scaffold arm.

4. **The miRNA sits where canonical Argonaute biology predicts.** In all four
   miR complexes the miR-nucleotide PLIP pass puts ~75–85 % of miR contacts on
   **PIWI** (~37–40/model) and **PAZ** (~10–19/model), with the **5′ / seed
   nucleotides (positions ~1–10)** carrying the most contact (top nt ≈ 4–7
   contacts/model), tapering toward the 3′ end. `MID` scores near-zero, but that
   is likely the domain-boundary definition (the 5′-nucleotide pocket residues
   fall inside the `PIWI` 678–1050 block as drawn) — read PIWI+MID together here.
   The miRNA makes **no direct contact with FBW2** in any model.

5. **Backend consensus is good for the small complexes, poor for the SCF ones.**
   `ago1_fbw2`, `ago1_fbw2_mir168`, `ago1_fbw2_mir393a` collapse to 2–3 pose
   clusters holding 70–95 % of models; `ago1_fbw2_ask1_cul1(_mir168)` fragment
   into 7–8 clusters. Clusters also tend to **segregate by backend family** (e.g.
   `ago1_fbw2` cluster 1 = Chai-1/Protenix/AF3, cluster 2 = RoseTTAFold3/OpenFold3)
   — the diffusion backends agree on the interface chemistry but not the exact
   rigid-body pose.

---

## 2. Pipeline status & caveats

| complex | selected | usable after energy filter | backends present | pose clusters | deliverables |
|---|---|---|---|---|---|
| ago1_fbw2 | 500 | 469 | AF3, Chai-1, OpenFold3, Protenix, RF3 | 3 | ✅ summary |
| ago1_fbw2_mir165a | 500 | 475 | 5 | 6 | ✅ summary + mir_ligands |
| ago1_fbw2_mir168 | 500 | 478 | 5 | 3 | ✅ summary + mir_ligands |
| ago1_fbw2_mir393a | 500 | 465 | 5 | 2 | ✅ summary + mir_ligands |
| ago1_fbw2_ask1_cul1 | 400 (399) | 367 | AF3, OpenFold3, Protenix, RF3 | 7 | ✅ summary (399/400) |
| ago1_fbw2_ask1_cul1_mir168 | 400 | 371 | 4 | 8 | ✅ summary + mir_ligands (399/400) |

**Caveats to keep in mind when reading the numbers:**

- **Boltz produced zero structures for all six complexes** (known backend bug in
  `run_boltz.py`'s OOM string-match — see HANDOFF.md). The ensembles are 4–5
  backends, not 6.
- **Chai-1 is absent from the two ASK1/CUL1 complexes** (2265 aa > Chai's
  2048-token cap) — those runs are 4 backends (AF3, OpenFold3, Protenix, RF3).
- **RoseTTAFold3 minimisation instability:** 27–40 % of RF3 models per complex
  blew up during ChimeraX minimisation (final energy positive / ~1e18–1e24
  kJ/mol) and are removed by the MAD energy filter (|mod-z| > 3.5). ~66–78 RF3
  models survive per complex (of 100). No whole backend was dropped.
- **1 model permanently lost:** `ago1_fbw2_ask1_cul1` /
  `rank_75_rosettafold3_seed11_sample1.0` — minimisation diverged to NaN in the
  first 2000 steps (a known RF3-input failure mode `minimize_cif.py` refuses to
  save). Its complex was aggregated over the remaining 399/400.
- **Contact *rates* are not domain-length-normalised** — big domains (PIWI 373
  aa, cullin_homology ~270 aa) accumulate more raw contacts. Compare the
  residue-level tables for specificity.
- **AGO1×miR shows "0 rows" in the main (anchor) PLIP pass** for every miR
  complex — that pass groups the miRNA with the AGO1 receptor. All miRNA contact
  data comes from the separate per-nucleotide `plip_mir_ligands` pass
  (`all_selected_summary_mir_ligands.csv`).

---

## 3. Per-complex detail

### 3.1 `ago1_fbw2` (baseline AGO1–FBW2)

- 82 114 contact rows / 500 models (469 kept). Backends: AF3 100, Chai-1 100,
  OpenFold3 99, Protenix 99, RF3 70.
- **3 pose clusters:** 1 (n=231, 46 % — Chai-1/Protenix/AF3-led), 2 (n=154, 31 %
  — RF3/OpenFold3-led), 3 (n=115, 23 % — mixed). PC1+PC2 = 73 % of variance.
- **Interaction mix:** H-bonds 33 261 · hydrophobic 21 283 · salt bridges
  14 486 · π-cation 1 501 · π-stacking 587. All AGO1(A)–FBW2(B).
- **Domain hotspots (contacts/model):** PIWI×C_tail **55.3**, N_ext×C_tail 14.9,
  PAZ×C_tail 14.7, ArgoN×C_tail 14.3, ArgoL2×C_tail 12.7, PAZ×LRR 8.4,
  N_ext×LRR 7.2. (F_box column ≈ 0.)
- **Top AGO1 interface residues:** R942, R222, R889, R973, R506, K728, K888,
  R438, R907, R340 (PIWI + PAZ + ArgoN, all basic).
- **Top FBW2 interface residues:** D280, E279, Y252, D264, E285, Y265, E292,
  W275, E288, E276 (C_tail, acidic + aromatic).
- **Most reproducible residue pairs:** AGO1 E469–FBW2 R222, AGO1 R474–FBW2 Y140,
  AGO1 K424/K436–FBW2 D264, AGO1 P390–FBW2 L169, AGO1 Y466–FBW2 F251/Y252
  (PAZ ↔ FBW2 LRR/C_tail junction).

### 3.2 `ago1_fbw2_mir165a` / `_mir168` / `_mir393a` (AGO1–FBW2 + one miRNA)

The AGO1–FBW2 interface is **the same interface as the baseline but weaker /
more variable** once the miRNA is present:

| | mir165a | mir168 | mir393a |
|---|---|---|---|
| models kept | 475 | 478 | 465 |
| pose clusters | 6 (2 dominant, 74 %) | 3 (all sizeable) | 2 (74 % in one) |
| PIWI×C_tail rate | 20.3 | 21.1 | 19.6 |
| PAZ×LRR rate | 8.0 | 5.3 | 7.8 |
| top AGO1 residues | R474, R907, K716, K418, R421 | R907, R474, K716, K418, Q720 | R474, R907, R421, K418, R283 |
| top FBW2 residues | W295, Y252, W313, E292, E288 | W295, Y252, W313, E292, E288 | Y252, W295, W313, Y265, D264 |

- FBW2 interface residues are **highly consistent across all three miRNAs and the
  baseline** (Y252, W295, W313, E288/292, the 250–300 acidic/aromatic stretch) —
  the FBW2 side of the interface is robust to which miRNA is loaded.
- AGO1 R474 / R907 / K716 / K418 recur as the top AGO1 anchors in every miR
  complex (PAZ R474/K418; PIWI R907; MID K716).
- **miRNA placement** (`plip_mir_ligands`, contacts/model): PIWI ≈ 39–40,
  N_ext_Grich ≈ 14–17, PAZ ≈ 10–19, ArgoN ≈ 10–13, ArgoL1/L2 ≈ 4–6, MID ≈ 0.
  Seed nucleotides (pos ~1–10) dominate; miR168 also shows a 3′-end (pos 21)
  contact. miR–FBW2 contacts: none.
- H-bond-dominated miRNA interface (miR: ~32 k H-bonds vs ~13 k salt bridges vs
  ~1 k hydrophobic per complex) — i.e. mostly backbone-phosphate / 2′-OH
  recognition, as expected for RNA in an Argonaute.

### 3.3 `ago1_fbw2_ask1_cul1` (+ SCF^FBW2 partial: ASK1/SKP1 + CUL1)

- 89 122 rows / 399 models (367 kept). Backends: AF3 100, OpenFold3 99,
  Protenix 98, RF3 67. **7 pose clusters**, the largest (n=163, 41 %) RF3/OpenFold3;
  next (n=104, 26 %) almost pure Protenix; next (n=74, 18 %) almost pure AF3 —
  **strong backend segregation, low consensus**. PC1+PC2 = 67 %.
- **AGO1–FBW2** (54 824 rows / 367 models): PIWI×C_tail **57.2**, ArgoN×C_tail
  15.5, ArgoL2×C_tail 13.7, PAZ×C_tail 12.4, N_ext×C_tail 10.3, PIWI×LRR 5.8.
  Top residue pairs move toward PIWI (AGO1 E852, Q857, K728, R907, R942 vs FBW2
  Y252, W225, F278, D296/297, D307) — the FBW2 grip spreads over more of PIWI in
  the bigger assembly.
- **AGO1–CUL1** (20 756 rows / 340 models): N_ext×cullin_homology **20.6**,
  N_ext×cullin_repeats 7.4, N_ext×cullin_CTD 7.4, PIWI×cullin_homology 3.9,
  PIWI×cullin_CTD 4.8, ArgoN×cullin_homology 2.2. → CUL1's N-terminal repeat arm
  lies along AGO1's Gly-rich N-extension and ArgoN.
- **AGO1–ASK1** (2 260 rows / 270 models, weak): ArgoN×SKP1_POZ 1.8,
  MID×SKP1_POZ 1.1, N_ext×SKP1_POZ 0.8 — glancing, non-specific.
- Interaction mix (all partners): H-bonds 39 103 · hydrophobic 21 756 ·
  salt bridges 15 043 · π-cation 1 492 · π-stacking 446.

### 3.4 `ago1_fbw2_ask1_cul1_mir168` (full assembly + miR168)

- 56 394 rows / 400 models (371 kept). Backends: AF3 100, OpenFold3 99,
  Protenix 97, RF3 71. **8 pose clusters**, largest n=157 (39 %, RF3/OpenFold3);
  n=89 (22 %) and n=73 (18 %) next. **PC1+PC2 only 49 %** — the most
  conformationally heterogeneous ensemble of the six (6 chains, most DOF).
- **AGO1–FBW2** (30 787 rows / 368 models): PIWI×C_tail 30.3, MID×C_tail 7.2,
  N_ext×C_tail 7.9, ArgoL2×C_tail 6.5, ArgoN×C_tail 4.7 — same interface,
  roughly half the contact density of the miR-free SCF run (the miRNA competes
  for AGO1 surface).
- **AGO1–CUL1** (17 149 rows / 361 models): N_ext×cullin_homology **17.4**,
  N_ext×cullin_repeats 4.8, N_ext×cullin_CTD 5.6, PIWI×cullin_homology 3.9 —
  CUL1 arm ↔ AGO1 N-lobe, same as 3.3.
- **AGO1–ASK1** (1 775 rows / 241 models): ArgoN/N_ext×SKP1_POZ ≈ 1.0, weak.
- **miR168** (`plip_mir_ligands`, 399 models): PIWI 37.2, PAZ 17.8, ArgoN 11.2,
  N_ext 9.2, ArgoL2 6.0; seed nt (pos ~1–10) dominant; no miR–FBW2/ASK1/CUL1
  contact.
- Interaction mix: H-bonds 25 695 · hydrophobic 14 114 · salt bridges 8 754 ·
  π-cation 956 · π-stacking 192.

---

## 4. Cross-complex synthesis

- **The AGO1–FBW2 interface is one interface**, present and chemically identical
  in all six complexes: FBW2's acidic/aromatic C-terminal solenoid+tail (res
  ~140–315; key residues Y252, W295, W313, F251, E279/288/292, D264/280) against
  AGO1's PAZ aromatic patch (Y466/471, F473, R474) plus a broad basic PIWI/N
  surface (R907, R942, R889, K716, K418, R474). The F-box is never engaged.
- **Adding a miRNA** roughly halves AGO1–FBW2 contact density and shifts weight
  from the focused PAZ patch toward the diffuse PIWI sheet, but does not move the
  interface. The miRNA occupies the canonical PIWI/PAZ nucleic-acid channel
  (seed-first) and is fully sequestered from FBW2.
- **Adding ASK1/CUL1** introduces a second, independent interface — CUL1
  N-terminal arm on AGO1's N-extension/ArgoN — and sharply increases pose
  heterogeneity. ASK1 does not contact AGO1. This is the geometry expected if
  SCF^FBW2 ubiquitinates AGO1 with FBW2 as the substrate-recognition subunit.
- **Confidence ranking of the poses:** `ago1_fbw2` ≈ `ago1_fbw2_mir393a` >
  `ago1_fbw2_mir168` ≈ `ago1_fbw2_mir165a` > `ago1_fbw2_ask1_cul1` >
  `ago1_fbw2_ask1_cul1_mir168` (by cluster concentration, backend agreement, and
  PCA variance captured).

---

## 5. Where things are

- **Per-model contact tables:** `results/<complex>/all_selected_summary.csv`
  (AGO1-anchored PLIP) and `..._mir_ligands.csv` (per-nucleotide miR PLIP, 4 miR
  complexes). 12 columns: `replica,model,resnr,restype,reschain,resnr_lig,
  restype_lig,reschain_lig,dist,ligcoo,protcoo,interaction_type`.
- **Pose clustering:** `results/<complex>/pose_clusters.csv` (pc1/pc2 + cluster +
  backend/seed/scores), `pose_clusters_pca.svg`, `pose_clusters_pca_variance.json`.
- **Domain-analysis figures:**
  `results/<complex>/figures/domain_analysis/` —
  `domain_contacts/*_heatmap_pooled.html`, `per_cluster/`, `per_backend/`,
  `residue_interface/*_top_residues.html` + `*_residue_contact_map.html`,
  `mir_ligands/mir_nt_x_<domain>_heatmap.html`, and `tables/*_rate_*.csv`
  (domain×domain contacts/model, pooled and per backend).
- **Executed notebooks (with outputs):**
  `notebooks/<complex>_domain_analysis.ipynb` (6) and
  `notebooks/pose_clustering.ipynb` (interactive GMM/HDBSCAN re-clustering, all 6
  complexes). Kernel used: `python3` from the `abcfold-drbs-notebook` micromamba
  env (the repo's own `abcfold-ago1-mirs-notebook` was never built; identical
  deps).
- **Cluster→structure symlinks:** `results/<complex>/reannotated/<method>/cluster_*/`
  (≤20 CIFs/cluster, for loading a whole cluster into ChimeraX/PyMOL).

---

## 6. Open items / suggested next steps

1. **`ago1_fbw2_ask1_cul1` is complete but its `all_selected_summary.csv` covers
   399/400 models** (1 RF3 model dropped). Cosmetic; note it in any figure
   captions. To fully close it, exclude
   `selected/cluster_7/rank_75_rosettafold3_seed11_sample1.0.cif` from the
   selection set and let snakemake finish cleanly.
2. **Recover Boltz** (optional, multi-day) — relax the OOM check in
   `abcfold/boltz/run_boltz.py` and re-run `--models b`. Would add a 6th backend
   to the consensus, most valuable for the two low-consensus SCF complexes.
3. **Domain boundaries for the miR analysis:** the MID/PIWI split (585/586,
   677/678) hides the 5′-nucleotide pocket inside PIWI. If MID-vs-PIWI
   partitioning of the miRNA matters, refine those two boundaries (or add a
   `MID_5p_pocket` sub-range) in `notebooks/generate_domain_notebooks.py`.
4. **RoseTTAFold3 minimisation blow-ups** (27–40 % of RF3 models) — worth a look
   at whether it's a systematic clash in RF3 outputs (e.g. RNA geometry) that a
   pre-minimisation fix could rescue, rather than discarding ~30 models/complex.
5. **Validate the two AGO1–FBW2 contact modes** against any experimental
   handle from `esther_project` (crosslinking, mutagenesis) — the PAZ↔FBW2-LRR
   aromatic patch is the falsifiable prediction.

---

*Working tree note: `notebooks/generate_domain_notebooks.py` + all 7 notebooks
are modified (generator bug-fixes for plotly ≥6 and the multi-complex
`pose_clustering.ipynb`; notebooks now carry executed outputs). Not committed —
review and commit at your discretion. Automation scripts for the weekend run are
in the session scratchpad (`weekend_driver.sh`, `weekend_finalize.sh`).*
