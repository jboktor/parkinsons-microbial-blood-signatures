# Analysis Plan: Validating PD-Associated CDR3β Signatures Against External Single-Cell and TCR Datasets

*Draft — March 30, 2026*

---

## Inputs from Unpublished Study

**Signature A — Meta-clonotype signatures:** ClustTCR-derived clusters of CDR3β sequences that are statistically enriched or depleted in PD vs. HC repertoires. Each cluster can be represented as the constituent exact sequences, and optionally converted to HMM motif models for flexible matching.

**Signature B — CDR3-risk-score (HLA-dosage motif):** A CDR3β motif whose prevalence correlates with cumulative dosage of protective HLA alleles in PD. A continuous score is assigned per TCR by existing scoring code.

**Scoring infrastructure:** The user has working code that accepts a CDR3β (+ V/J gene) and returns a continuous score for Signature B, and cluster membership for Signature A.

---

## Question 1: Which T-cell subtypes are enriched for each signature?

### Rationale

If the CDR3β signatures reflect a biologically meaningful immune process in PD, they should not be uniformly distributed across all T-cell states. Understanding whether the signatures concentrate in naive vs. antigen-experienced cells, CD4 vs. CD8 compartments, or regulatory T cells constrains the mechanistic interpretation (e.g., autoimmune effector response vs. failed tolerance).

### Required data type

Paired scRNA-seq + TCR-seq at the single-cell level, with T-cell subtype annotations. This is the only way to link a CDR3β sequence to a specific cellular phenotype.

### Dataset selection

| Priority | Study | Why | Cells with paired GEX+TCR | Subtypes | Tissue | Accession |
|----------|-------|-----|---------------------------|----------|--------|-----------|
| **PRIMARY** | Wang et al. 2021 (Study 1) | Largest paired 5' GEX + V(D)J TCR dataset focused on T cells in PD; 21 T-cell subsets annotated; includes naive, CM, EM, TEMRA, CTL, Treg; HLA data available for cross-referencing with Signature B | ~84,000 | 21 T-cell subsets | Blood + CSF | GSE134578 (CSF); new blood data in paper |
| **SECONDARY** | Moquin-Beaudry et al. 2025 (Study 2) | Independent paired scRNA + TCR-seq; 38 immune subtypes; different cohort (Montreal) provides replication | ~78,876 (fraction with TCR) | 38 subtypes | Blood | Check paper |
| **TERTIARY (brain)** | Ran/Ma et al. 2025 (Study 4) | Only dataset with paired snRNA-seq + TCR in brain tissue; small T-cell numbers but unique tissue | ~535 T cells (PD SN) | T-cell states in tissue context | Brain (SN) | GSE303705 |

**Not suitable:** Studies 3, 6, 8, 9 (no TCR-seq); Studies 7, 11 (reanalyses without TCR); Studies 12a, 12b (brain snRNA without TCR, T cells are minor populations).

### Analysis workflow

**Step 1 — Data acquisition and TCR extraction.**
Obtain the processed Seurat/AnnData objects and TCR contig files from Wang et al. and Moquin-Beaudry et al. For Wang et al., the blood data may require contacting the authors or accessing the GSA/CNCB archive; CSF TCR data is at GSE134578. For Ran/Ma et al., download from GSE303705 and the Google Drive link.

**Step 2 — Score every TCR-bearing cell.**
For each cell that has both a gene-expression profile and a productive CDR3β chain:

- Compute the **CDR3-risk-score** (Signature B) using existing scoring code.
- Determine **meta-clonotype membership** (Signature A) by checking whether the CDR3β falls within any ClustTCR cluster (exact match or HMM-based scoring).

Attach these scores as metadata columns to the single-cell object.

**Step 3 — Subtype-level enrichment analysis.**

For **Signature B (continuous score)**:

- Compare CDR3-risk-score distributions across annotated T-cell subtypes using a Kruskal-Wallis test (non-parametric; score distributions will not be normal).
- Follow up with pairwise Dunn's tests (BH-corrected) to identify which subtypes differ.
- Visualize with violin plots overlaid on the UMAP, stratified by subtype.

For **Signature A (meta-clonotype membership, binary)**:

- For each meta-clonotype cluster (or aggregated enriched/depleted sets), compute the proportion of cells belonging to each T-cell subtype that match.
- Compare against the expected proportion under the null (proportion of all TCR-bearing cells in that subtype). Use a Fisher's exact test per subtype or a chi-squared test across all subtypes.
- Compute odds ratios with 95% CI for each subtype.

**Step 4 — Address sub-questions directly.**

*Q1a: Naive vs. antigen-experienced.*

- Dichotomize cells into "naive" (CD4 naive, CD8 naive, recent thymic emigrants) vs. "antigen-experienced" (all memory, effector, TEMRA, CTL subtypes) based on the study's own annotations.
- Wilcoxon rank-sum test comparing CDR3-risk-score between the two groups.
- Fisher's exact test for meta-clonotype enrichment (proportion matching in naive vs. experienced).
- **Expected outcome:** If the signature is disease-relevant, antigen-experienced cells should be enriched (they have undergone clonal selection). A null result (equal distribution) would suggest the motif is germline-encoded rather than selection-driven.

*Q1b: CD4 vs. CD8 vs. Treg.*

- Group cells into CD4 conventional (Th1, Th2, Th17, Tfh, CM, EM), CD8 (naive, CM, EM, TEMRA, CTL), and Treg.
- Kruskal-Wallis across these three broad compartments, then pairwise Dunn's tests.
- Of particular interest: CD4 CTLs were identified as expanded in PD by Wang et al. — test whether these carry disproportionately high signature scores.
- For Tregs specifically, test whether signature-positive Tregs differ transcriptionally from signature-negative Tregs (differential expression within Treg cluster, conditioned on signature status).

**Step 5 — Replication.**
Repeat Steps 2–4 using Moquin-Beaudry et al. as an independent validation cohort. Compare effect sizes and direction of enrichment. A meta-analytic combination (fixed-effects or random-effects on the odds ratios) can be computed if both datasets yield the same subtype annotations.

### Power considerations

Wang et al. provides ~84,000 cells with both modalities. Even if only 50% carry a productive CDR3β, ~42,000 cells across 21 subtypes gives excellent power to detect moderate enrichment (OR > 1.5) in any subtype containing >500 cells. The Ran/Ma brain dataset (535 T cells) will have limited power for subtype stratification but can detect large effects.

---

## Question 2: Do T cells with known reactivity to PD-associated antigens harbor enrichment for the signatures?

### Rationale

If the CDR3β signatures capture TCRs that are functionally relevant to PD pathogenesis, then TCRs with experimentally verified reactivity to PD-relevant antigens should be enriched for the signature. Recent work has expanded the target antigen landscape well beyond α-synuclein: Williams et al. (JCI, 2024) demonstrated that PINK1 is a major T cell target in PD, with additional responses to GBA, SOD1, Parkin, OGDH, and LRRK2. The mechanistic basis was established by Matheoud et al. (Cell, 2016), who showed that PINK1/Parkin deficiency leads to mitochondrial antigen presentation (MitAP) on MHC-I, directly linking mitochondrial dysfunction to autoimmune T cell activation. Testing across this full antigen panel is critical — a signature that tracks with multiple PD antigens (not just α-syn) provides much stronger evidence of disease-relevant immune recognition.

### Required data type

TCR sequences from antigen-specific T cells with verified reactivity to PD-relevant peptides, plus background TCRs. Alternatively, T cell reactivity measurements (FluoroSPOT, ICS) with donor HLA typing permit indirect analysis of HLA-restricted epitope–signature relationships.

### Dataset selection

| Priority | Study | Antigens tested | What it provides | Limitation |
|----------|-------|-----------------|------------------|------------|
| **PRIMARY (TCR sequences)** | Singhania et al. 2021 (Study 15) + Lindestam Arlehamn et al. 2020 (Study 5) | α-synuclein | CDR3β sequences of α-syn-reactive T cells from 6 HLA-typed PD donors; AIM-sorted single-cell TCR-seq; pertussis-specific TCRs as negative control | 6 donors; sequences may need author contact |
| **PRIMARY (epitope–HLA mapping)** | Williams et al. 2024, JCI (Study 14) | PINK1, GBA, SOD1, Parkin, OGDH, LRRK2 | FluoroSPOT reactivity for 6 antigen pools; 34 PINK1 epitopes mapped to HLA class II restrictions (DRB1*15:01, DRB4*01:01, DRB1*04, DPB1*04:01, DQB1*03:02); 39 PD + 39 HC; HLA-typed donors | No TCR sequences (FluoroSPOT only); but HLA restriction patterns enable indirect signature connection |
| **SECONDARY (discovery)** | Sulzer et al. 2017 (Study 13) | α-synuclein | Original discovery: α-syn T cell reactivity in 67 PD + 36 HC; HLA-DRB1*15:01/DRB5*01:01 associations | ICS data, no TCR sequences |
| **SECONDARY (prodromal)** | Johansson et al. 2025 (Study 16) | PINK1, α-synuclein | PINK1 + α-syn reactivity in prodromal PD (iRBD); temporal staging data | FluoroSPOT only; HLA typing status unconfirmed |
| **SECONDARY (brain homology)** | Ran/Ma et al. 2025 (Study 4) | α-synuclein (homology) | Brain-infiltrating TCR clonotypes annotated for α-syn-reactive sequence homology | Homology-based, not direct experimental verification |
| **SUPPORTING (mechanistic)** | Garretti et al. 2023 (Study 17) | α-syn (HLA-DRB1*15:01-restricted) | In vivo causal demonstration: α-syn32-46 + HLA-DRB1*15:01 triggers enteric PD features; CD4+ Th1/Th17 | Mouse model — provides epitope-level mechanistic context, not human TCR data |
| **SUPPORTING (mechanistic)** | Matheoud et al. 2016 (Study 18) | Mitochondrial antigens (PINK1/Parkin pathway) | PINK1/Parkin deficiency → MitAP → MHC-I presentation → CD8+ T cell activation | Mechanistic framework; no human TCR sequences |
| **EXTERNAL DB** | VDJdb / McPAS-TCR / IEDB | Various PD-relevant | Any deposited α-syn-reactive or PD-associated TCR sequences | Sparse coverage for PD antigens |

### Analysis workflow

**Step 1 — Assemble antigen-specific TCR sequences.**

*1a — α-synuclein-reactive CDR3β sequences (primary set):*
- Obtain CDR3β sequences from Singhania et al. 2021 (Study 15) supplementary tables — these are the actual single-cell TCR-seq data from AIM-sorted α-syn-reactive T cells across 6 PD donors. This is the Sette lab's follow-up to Lindestam Arlehamn 2020 with the actual TCR sequences.
- Obtain the pertussis-specific TCRs from the same donors (Lindestam Arlehamn 2020, Study 5) — critical negative control (foreign-antigen-reactive TCRs from the same individuals).
- Query VDJdb (vdjdb.cdr3.net) and McPAS-TCR for additional α-syn-reactive or PD-associated TCR entries.
- From Ran/Ma et al. (Study 4), extract brain TCR clonotypes reported to have homology to known α-syn-reactive sequences.

*1b — PINK1 epitope–HLA restriction map (for indirect analysis):*
- From Williams et al. 2024 (Study 14) supplementary data, extract the 34 mapped PINK1 epitope–HLA restriction pairs. Key restrictions: DRB1*15:01, DRB4*01:01, DRB1*04:01/04:04, DPB1*04:01, DQB1*03:02.
- Note the overlap: DRB1*15:01 restricts both α-syn epitopes (Sulzer 2017, Garretti 2023) and PINK1 epitopes (Williams 2024). DRB1*04:01/04:04 are protective alleles in PD and also restrict PINK1 epitopes. This directly connects the HLA-dosage CDR3β signature (Signature B) to the antigen-specificity question.

**Step 2 — Score antigen-specific TCRs.**

- Apply CDR3-risk-score (Signature B) to each α-syn-specific TCR from Singhania/Lindestam Arlehamn.
- Determine meta-clonotype membership (Signature A) for each.

**Step 3 — Construct the background distribution.**

Three complementary approaches:

*Approach A — Within-donor background (preferred):*
For each donor in the Lindestam Arlehamn/Singhania studies, use the full unsorted repertoire (if available) as background. Score all background TCRs and compare to the antigen-specific subset from the same donor. This controls for HLA-driven repertoire shaping.

*Approach B — Permutation test against matched repertoires:*
From paired scRNA+TCR datasets (Wang et al. or Moquin-Beaudry et al.), randomly sample N TCRs (N = number of α-syn-specific TCRs) from PD patients. Compute mean CDR3-risk-score. Repeat 10,000 times for null distribution. Report empirical p-value.

*Approach C — Foreign-antigen control:*
Compare α-syn-specific TCR scores to pertussis-specific TCRs from the same donors. Paired Wilcoxon rank-sum test. Effect size: Cohen's d or rank-biserial correlation.

**Step 4 — Meta-clonotype overlap analysis.**

- Compute the overlap: what fraction of α-syn-specific TCRs fall within any PD-enriched meta-clonotype? Any PD-depleted meta-clonotype?
- Compare to the expected overlap rate under the null (from Step 3).
- One-sided Fisher's exact test: are α-syn-specific TCRs over-represented in PD-enriched meta-clonotypes?
- Separately test under-representation in PD-depleted meta-clonotypes.

**Step 5 — HLA restriction–signature convergence analysis (NEW).**

This analysis leverages the remarkable convergence between the HLA-dosage CDR3β signature and the known HLA restrictions of PD antigen epitopes:

*5a — HLA-stratified α-syn enrichment:*
- Among α-syn-specific TCRs from donors carrying protective HLA alleles (DRB1*04:01/04:04), test whether CDR3-risk-score is higher or lower than from donors without. This directly connects Signature B to the HLA-dosage model.
- Among donors carrying risk alleles (DRB1*15:01/DRB5*01:01), test the converse.

*5b — Cross-antigen HLA restriction analysis (Williams 2024 data):*
- The Williams et al. epitope–HLA restriction map reveals that:
  - DRB1*15:01 restricts both α-syn AND PINK1 epitopes → risk allele
  - DRB4*01:01 restricts PINK1 epitopes → linked to DRB1*04 (protective haplotype)
  - DRB1*04:01/04:04 restrict PINK1 epitopes → protective alleles
- Test whether the CDR3-risk-score differentially tracks with TCRs restricted by risk vs. protective HLA alleles. Specifically: if TCRs restricted by protective HLA alleles (DRB1*04-restricted PINK1 epitopes) show different signature scores than TCRs restricted by risk alleles (DRB1*15:01-restricted epitopes), this would directly validate the HLA-dosage mechanism.
- This analysis may be indirect (matching HLA types of donors to their reactivity profiles) since Williams et al. used FluoroSPOT rather than TCR-seq. However, the HLA typing + reactivity data allows donor-level correlation.

*5c — Epitope-stratified analysis:*
- Stratify by epitope restriction class: MHC-I-restricted (CD8-presented) vs. MHC-II-restricted (CD4-presented) epitopes.
- Matheoud et al. (2016) showed PINK1/Parkin deficiency leads to MHC-I presentation of mitochondrial antigens → CD8+ T cells. Williams et al. (2024) found predominantly MHC-II-restricted PINK1 epitopes → CD4+ T cells. Test whether the CDR3β signatures differentially enrich for CD4- vs. CD8-restricted antigen-specific TCRs (complementing Q1b).

**Step 6 — Cross-antigen comparison (NEW).**

Using donor-level data from Williams et al. 2024 (39 PD + 39 HC with FluoroSPOT across 6 antigen pools):

- For donors who also appear in paired scRNA+TCR datasets (if any overlap exists), test whether donors with strong PINK1/GBA/LRRK2 T cell responses have repertoires enriched for the CDR3β signatures compared to non-responders.
- Even without donor overlap, the Williams et al. HLA-typed cohort allows testing whether donors carrying protective HLA alleles show different antigen reactivity profiles — connecting HLA-dosage to antigen specificity at the cohort level.
- Test the hierarchy: does the signature correlate more strongly with PINK1 reactivity (dominant response, 42.5% of PD) than with less frequent responses (GBA, LRRK2, SOD1)?

### Power considerations

The Singhania et al. 2021 study identified TCRs from 6 PD donors with "surprisingly diverse" repertoires — if ~50–200 unique α-syn-specific CDR3β sequences are available, permutation approaches will have adequate power to detect a 0.5 SD shift in mean CDR3-risk-score. The pertussis control analysis is well-powered (same donors, paired design). The Williams et al. 2024 cohort (39+39) provides adequate power for donor-level HLA-stratified analyses of reactivity profiles. The cross-antigen HLA restriction analysis is hypothesis-generating but builds on established restriction patterns from multiple independent studies.

---

## Question 3: Do TCRs from PD-relevant tissues (substantia nigra, CSF, brain) harbor enrichment for the signatures?

### Rationale

If the CDR3β signatures mark T cells that participate in the neuroinflammatory process, then TCRs recovered from the actual sites of neurodegeneration (substantia nigra, CSF) should be enriched compared to peripheral blood. This tests whether the peripherally-defined signature has tissue-infiltrating relevance.

### Required data type

TCR sequences from brain tissue (SN, striatum) or CSF, ideally with matched blood from the same donors for paired comparison.

### Dataset selection

| Priority | Study | Why | Tissue TCRs available | Paired blood? | Accession |
|----------|-------|-----|----------------------|---------------|-----------|
| **PRIMARY** | Ran/Ma et al. 2025 (Study 4) | Only dataset with TCR-seq from post-mortem substantia nigra; clonally expanded CD8+ T cells; TCRs with α-syn homology | TCR clonotypes from SN (535 T cells in PD); also cingulate cortex | No (post-mortem brain only) | GSE303705 + GitHub |
| **PRIMARY** | Wang et al. 2021 (Study 1) | Paired blood + CSF TCR-seq from same PD patients; direct within-patient tissue comparison | ~84K cells with TCR across blood and CSF | Yes (blood + CSF from overlapping patients) | GSE134578 (CSF) |
| **SUPPORTING** | Gate/Cantoni et al. 2025 (Study 10) | Large CSF scRNA-seq compendium; no TCR-seq, but TRUST4 computational TCR extraction from 5' scRNA-seq reads may be feasible | Potential for computationally extracted CDR3β from CSF 5' scRNA-seq | Yes (CSF + blood) | Synapse syn51730532 |
| **SUPPORTING** | Guan et al. 2022 (Study 11) / GSE141578 | CSF scRNA-seq (10x 5'); source data for computational TCR extraction if 5' chemistry was used | 18,553 CSF cells; potential CDR3β extraction | No | GSE141578 |

**Brain snRNA-seq without TCR (Studies 12a, 12b):** These capture infiltrating T cells as minor populations but used 3' chemistry snRNA-seq, making computational TCR extraction infeasible. Not suitable for this question.

### Analysis workflow

**Step 1 — Assemble tissue-derived TCR sequences.**

*From Ran/Ma et al. (brain):*
- Download processed data from GSE303705 and/or the Google Drive link.
- Extract all TCR clonotypes recovered from substantia nigra and cingulate cortex. These include clonally expanded CD8+ T cells.
- Separately flag TCRs the authors annotated as having homology to α-syn-reactive clones.

*From Wang et al. (CSF):*
- Obtain CSF TCR data from GSE134578.
- Obtain paired blood TCR data from the same patients (contact authors or check GSA/CNCB).
- This is uniquely valuable because it enables within-patient tissue comparisons.

*From Gate/Cantoni et al. (CSF, computational):*
- If the original data used 10x 5' chemistry (confirmed for new samples), apply TRUST4 or scRepertoire to computationally extract CDR3β sequences from the BAM/FASTQ files at Synapse syn51730532.
- This is a higher-effort step but yields CSF TCRs from a large multi-disease cohort (6 PD patients).

**Step 2 — Score all tissue-derived TCRs.**

Apply both the CDR3-risk-score and meta-clonotype membership scoring to every tissue-derived TCR.

**Step 3 — Statistical tests.**

*Analysis 3A — Within-patient tissue comparison (Wang et al.):*
This is the most powerful design because it controls for donor-level confounders (genetics, disease stage, HLA).

- For each PD patient with both blood and CSF TCR data, compute the mean CDR3-risk-score in CSF TCRs vs. blood TCRs.
- Paired Wilcoxon signed-rank test across patients (H₀: no difference in score between compartments).
- Compute the proportion of CSF TCRs falling in PD-enriched meta-clonotypes vs. the proportion of blood TCRs from the same patient. Paired McNemar-style test or paired proportion test.
- Visualize with paired dot plots (blood vs. CSF per patient).

*Analysis 3B — Brain vs. external blood reference (Ran/Ma et al.):*
Since Ran/Ma et al. lacks paired blood, use blood TCRs from Wang et al. or Moquin-Beaudry et al. as the reference.

- Compare CDR3-risk-score distribution of SN-derived TCRs vs. blood-derived TCRs using a Wilcoxon rank-sum test.
- Bootstrap confidence interval on the difference in means to quantify effect size.
- For meta-clonotype membership: Fisher's exact test comparing the proportion of SN TCRs in enriched meta-clonotypes vs. blood TCRs.
- **Caveat:** This comparison is cross-cohort and cannot control for donor-level effects. Report effect sizes alongside p-values.

*Analysis 3C — Clonal expansion interaction:*
Within the Ran/Ma brain data and the Wang CSF data, clonally expanded TCRs (clone size > 1) are already identified.

- Test whether clonally expanded tissue TCRs have higher CDR3-risk-scores than singleton tissue TCRs (Wilcoxon rank-sum).
- Test whether the most expanded clones (top 10% by clone size) are more likely to fall within PD-enriched meta-clonotypes (Fisher's exact).
- This addresses whether the signature marks the specific clones undergoing active selection in PD-relevant tissue.

*Analysis 3D — α-syn-homologous brain TCRs (Ran/Ma et al.):*
Ran/Ma et al. specifically flagged brain TCRs with sequence homology to known α-syn-reactive clones.

- Compare CDR3-risk-scores of α-syn-homologous brain TCRs vs. non-homologous brain TCRs (Wilcoxon).
- This connects Q2 (antigen specificity) with Q3 (tissue localization) in the most disease-relevant context.

### Power considerations

The Ran/Ma brain dataset has ~535 T cells from PD substantia nigra with TCR data. This is small but sufficient to detect large effects (OR > 2.0 for meta-clonotype enrichment). The Wang et al. CSF dataset is larger and the paired design greatly increases power. The computational TCR extraction from Gate/Cantoni adds volume but introduces noise from the extraction algorithm.

---

## Execution Order and Dependencies

The analyses above have logical dependencies. The recommended execution order:

```
Phase 0: Data acquisition (parallel)
  ├─ Obtain Wang et al. blood + CSF paired data (contact authors / GSA)
  ├─ Download Ran/Ma et al. from GSE303705 + Google Drive
  ├─ Download Moquin-Beaudry et al. (check Brain article for accession)
  ├─ Obtain α-syn-specific CDR3β sequences from Singhania et al. 2021 supplementary data
  ├─ Contact Lindestam Arlehamn / Sette lab for pertussis-specific TCR controls
  ├─ Extract PINK1 epitope–HLA restriction map from Williams et al. 2024 (JCI) supplementary data
  ├─ Extract α-syn epitope–HLA associations from Sulzer et al. 2017 and Garretti et al. 2023
  ├─ Query VDJdb + McPAS-TCR for PD-relevant TCRs
  └─ Optionally: download Gate/Cantoni from Synapse syn51730532

Phase 1: Scoring (depends on Phase 0)
  └─ Score ALL external TCRs with both Signature A and Signature B

Phase 2: Q1 — Subtype enrichment (depends on Phase 1)
  ├─ Primary analysis: Wang et al. paired GEX+TCR
  ├─ Replication: Moquin-Beaudry et al.
  └─ Brain context: Ran/Ma et al. (small N)

Phase 3: Q2 — Antigen specificity (depends on Phase 1; parallel with Phase 2)
  ├─ Score α-syn-specific TCRs (Singhania 2021 + Lindestam Arlehamn 2020)
  ├─ Permutation test against background
  ├─ Pertussis negative control comparison
  ├─ HLA restriction–signature convergence (risk vs. protective allele restrictions)
  ├─ Cross-antigen HLA restriction analysis (Williams 2024 PINK1 epitope map)
  ├─ Epitope-stratified analysis: MHC-I (PINK1/Parkin MitAP) vs. MHC-II (α-syn, PINK1)
  └─ Cross-antigen reactivity hierarchy (PINK1 > GBA > LRRK2 > SOD1)

Phase 4: Q3 — Tissue enrichment (depends on Phase 1; parallel with Phases 2–3)
  ├─ Within-patient CSF vs. blood (Wang et al.)
  ├─ Brain vs. external blood (Ran/Ma et al.)
  ├─ Clonal expansion interaction
  └─ α-syn-homologous brain TCR sub-analysis (connects Q2 + Q3)

Phase 5: Integration and multiple testing correction
  ├─ Apply BH-FDR correction across all tests within each question
  ├─ Summarize effect sizes and confidence intervals
  ├─ Build integrated figure panel
  └─ Compile HLA restriction convergence narrative (connecting Signature B
     to antigen-specific HLA restrictions across α-syn, PINK1, and MitAP pathways)
```

---

## Multiple Testing Correction Strategy

Within each question, apply Benjamini-Hochberg FDR correction at q < 0.05. Across questions, treat Q1/Q2/Q3 as separate hypothesis families (no cross-question correction needed since they address distinct biological questions). Within Q1, the subtype comparisons (21 subtypes × 2 signature types) constitute ~42 tests; within Q2, the expanded antigen panel increases tests to ~12–18 (core antigen enrichment tests + HLA restriction stratification + cross-antigen comparisons); within Q3, the tissue comparisons constitute ~6–8 tests.

For permutation-based tests (Q2 background comparisons), report empirical p-values directly — these are already multiple-comparison-aware within the permutation framework.

Note on the Q2 expansion: the addition of HLA restriction–signature convergence analyses (Step 5) and cross-antigen comparisons (Step 6) are partially hypothesis-generating. Pre-register the core antigen enrichment tests (Steps 2–4) as confirmatory, and clearly label the HLA restriction and cross-antigen analyses as exploratory/hypothesis-generating.

---

## Summary: Dataset-to-Question Mapping

| Dataset | Q1 (Subtypes) | Q2 (Antigen specificity) | Q3 (Tissue) | Priority |
|---------|:-:|:-:|:-:|----------|
| **Wang et al. 2021** (Study 1) | **PRIMARY** (paired GEX+TCR, 21 subtypes) | Background repertoire | **PRIMARY** (paired blood vs. CSF) | Essential |
| **Moquin-Beaudry 2025** (Study 2) | **REPLICATION** (paired GEX+TCR) | Background repertoire | — | High |
| **Ran/Ma 2025** (Study 4) | Tertiary (brain, small N) | α-syn-homologous TCRs | **PRIMARY** (brain TCRs) | Essential |
| **Singhania 2021** (Study 15) | — | **PRIMARY** (α-syn CDR3β sequences) | — | Essential |
| **Lindestam Arlehamn 2020** (Study 5) | — | **PRIMARY** (HLA-typed donors + pertussis control) | — | Essential |
| **Williams et al. 2024, JCI** (Study 14) | — | **PRIMARY** (PINK1/GBA/SOD1/Parkin/OGDH/LRRK2 + HLA restrictions) | — | Essential |
| **Sulzer et al. 2017** (Study 13) | — | Supporting (α-syn discovery, HLA-DRB1*15:01) | — | High |
| **Johansson et al. 2025** (Study 16) | — | Supporting (prodromal PINK1 + α-syn) | — | Recommended |
| **Garretti et al. 2023** (Study 17) | — | Supporting (mechanistic: α-syn–HLA-DRB1*15:01 causal) | — | Context |
| **Matheoud et al. 2016** (Study 18) | — | Supporting (mechanistic: PINK1/Parkin MitAP → MHC-I) | — | Context |
| **Gate/Cantoni 2025** (Study 10) | — | — | Supporting (CSF, computational TCR extraction) | Optional |
| **Guan 2022 / GSE141578** (Study 11) | — | — | Supporting (CSF, computational TCR extraction) | Optional |
| VDJdb / McPAS-TCR / IEDB | — | Supporting (curated antigen-specific TCRs) | — | Recommended |

### Essential datasets (6 total)

The six essential datasets are **Wang et al. 2021**, **Moquin-Beaudry et al. 2025**, **Ran/Ma et al. 2025**, **Singhania et al. 2021**, **Lindestam Arlehamn et al. 2020**, and **Williams et al. 2024 (JCI)**. Together they provide: paired GEX+TCR from blood and CSF (Q1, Q3), brain-infiltrating TCRs (Q3), α-syn-specific CDR3β sequences with HLA context (Q2), and the broadest PD antigen–HLA restriction map available (Q2), covering all three questions with statistical rigor.

### Key HLA convergence points for Signature B validation

The HLA-dosage CDR3β signature (Signature B) is associated with protective HLA alleles. The antigen-specific literature reveals a striking convergence:

- **DRB1*15:01 (risk allele):** Restricts α-syn epitopes Y39 and S129-phospho (Sulzer 2017), restricts α-syn32-46 causing enteric PD features (Garretti 2023), AND restricts PINK1 epitopes (Williams 2024)
- **DRB1*04:01/04:04 (protective alleles):** Restrict PINK1 epitopes (Williams 2024)
- **DRB4*01:01:** Restricts PINK1 epitopes; in LD with DRB1*04 haplotype (Williams 2024)
- **DRB5*01:01 (risk allele):** Associated with α-syn T cell reactivity (Sulzer 2017); in LD with DRB1*15:01

This means the same HLA alleles that define the CDR3β signature's association with PD risk/protection are the alleles that restrict presentation of α-syn and PINK1 epitopes to T cells. Testing whether signature-positive TCRs are enriched among antigen-specific T cells restricted by these specific alleles is the most direct validation of the HLA-dosage model.