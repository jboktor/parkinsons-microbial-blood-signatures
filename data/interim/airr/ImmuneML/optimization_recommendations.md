# ImmuneML Transfer-Study Optimization Recommendations (Dataset-Aware)

**Date:** 2026-04-14
**Scope:** 95/96 completed jobs from the S2S + LOSO matrix

---

## 0. Dataset Structure — Hard Constraints for Tuning

Every recommendation below is grounded in the actual shape of this dataset:

| Study | N samples | Clonal volume (min / median / max) |
|---|---|---|
| **BioFIND** | 174 | 501 / 1077 / 4661 |
| **PDBP** | 1234 | 500 / 1916 / 16000 |
| **PPMI** | 1421 | 500 / 3088 / 14896 |
| **HBS** | 746 | 936 / 3903 / 8464 |

- Every repertoire was **pre-filtered to clonal_volume ≥ 500** (see `AIRR_ml_prep.qmd` line 304).
- **Median depth 1,000-4,000 sequences/repertoire.** This is ~100× shallower than the datasets DeepRC / CompAIRR were designed for (Emerson 2017 had ~10⁵ seqs/repertoire; Rawat 2020 used ~10⁶).
- Metadata includes binned covariates: `sex`, `age_at_baseline_binned` (1-5), `clonal_volume_binned` (1-5), `ethnicity`, `race` — use these for `PerformancePerLabel` to detect confounders.

**Implications for hyperparameters:**
1. DeepRC `sample_n_sequences` must be ≤ minimum repertoire size (~500). My earlier suggestion of 10,000 was wrong — that exceeds every repertoire.
2. Nested 5×5 CV on BioFIND means ~78 train / 20 val / 20 test per fold. DeepRC will overfit badly at this scale — restrict deep methods to PDBP/PPMI/HBS or LOSO folds that combine them.
3. CompAIRR `p_value_threshold=0.001` is statistically unreachable when each repertoire only contains ~1k unique clones — effectively no feature can pass Bonferroni-corrected significance.

---

## 1. CRITICAL: Current Transfer-Study Design Bug

Every current S2S config defines both `train_dataset` (e.g., PDBP) and `test_dataset` (e.g., PPMI), but the `TrainMLModel` instruction only references:

```yaml
instructions:
  S2S_pdbp_to_ppmi_kmer_logreg:
    type: TrainMLModel
    dataset: train_dataset    # ← only train_dataset is used!
    assessment:
      split_strategy: random
      training_percentage: 0.7  # ← 70/30 random split WITHIN PDBP
```

**What this means:** The "PDBP → PPMI" run is actually *PDBP train / 30% of PDBP test* — **it never touches PPMI at all.** The `test_dataset` block is imported but ignored. So our entire "transfer-study matrix" is currently measuring within-cohort generalization, not cross-cohort transfer.

### Fix: use `MLApplication` as a second instruction, chained after `TrainMLModel`

```yaml
instructions:
  train_on_pdbp:
    type: TrainMLModel
    dataset: pdbp_dataset
    assessment:
      split_strategy: stratified_k_fold
      split_count: 5
    selection:
      split_strategy: stratified_k_fold
      split_count: 5
    labels: [{case_control_other_latest: {positive_class: Case}}]
    # …

  apply_to_ppmi:
    type: MLApplication
    dataset: ppmi_dataset
    config_path: train_on_pdbp/optimal_case_control_other_latest/  # path to trained model
    label: case_control_other_latest
    metrics: [auc, balanced_accuracy, precision, recall, f1_micro]
```

This runs CV on the train study to pick the best model, then **actually applies it** to the held-out test study. Without this fix, **the transfer-study results so far are not transfer results** — they're within-PDBP, within-PPMI, etc., labelled with the wrong destination.

**Priority: fix this before spending more GPU time on DeepRC.**

---

## 2. Per-Model Recommendations (Dataset-Grounded)

### 2A. DeepRC — needs data-scaled retraining

| Param | Current | Recommended | Rationale |
|---|---|---|---|
| `n_updates` | 10,000 | **100,000** | Paper used 10⁵ as minimum; 3×10⁵ is overkill for small repertoires |
| `evaluate_at` | 2,000 | 5,000 | Paper default |
| `sample_n_sequences` | 2,000 | **500** | Must be ≤ smallest repertoire (clonal_volume floor is 500) |
| `kernel_size` | 5 | `[5, 7]` sweep | Start shorter — our CDR3β ~8-11 AAs |
| `n_kernels` | 16 | 32 | Paper optimum |
| `n_additional_convs` | 1 | 1 | Keep shallow to prevent overfit at ~1000-sample train sets |
| `l2_weight_decay` | 0.001 | **0.001 (keep)** | Strong regularization IS warranted at this depth |
| `learning_rate` | 5e-5 | 5e-5 | Keep |

**Do NOT run DeepRC on BioFIND-only training** (174 samples × 0.7 = 122; nested CV gets ~78 train). It will overfit. Restrict DeepRC to:
- Single-study training: PDBP, PPMI, HBS (all ≥ 746 samples)
- All 4 LOSO configs (train sets 2154-3401 samples)

That's 3 S2S-source studies × 3 test targets + 4 LOSO = **13 DeepRC runs instead of 16**.

### 2B. K-mer Models — Keep & Extend

K-mer models are the only ones showing real signal (0.43-0.60 range). Extend the encoding space:

```yaml
encodings:
  kmer3_all:         # current, keep
    KmerFrequency: {k: 3, sequence_encoding: continuous_kmer, normalization_type: relative_frequency, reads: all, scale_to_unit_variance: true, scale_to_zero_mean: true, sequence_type: amino_acid, region_type: IMGT_CDR3}
  kmer4_all:         # add — longer motifs
    KmerFrequency: {k: 4, sequence_encoding: continuous_kmer, normalization_type: relative_frequency, reads: all, scale_to_unit_variance: true, scale_to_zero_mean: true, sequence_type: amino_acid, region_type: IMGT_CDR3}
  gapped_2_2:        # add — non-contiguous motifs
    KmerFrequency: {k_left: 2, k_right: 2, min_gap: 0, max_gap: 3, sequence_encoding: gapped_kmer, normalization_type: relative_frequency, sequence_type: amino_acid, region_type: IMGT_CDR3}
```

For ML methods, add elastic net:

```yaml
log_reg_elastic:
  LogisticRegression:
    penalty: elasticnet
    l1_ratio: [0.1, 0.5, 0.9]
    C: [0.001, 0.01, 0.1, 1.0, 10.0]
    max_iter: 10000
    solver: saga
  model_selection_cv: true
  model_selection_n_folds: 5
```

### 2C. CompAIRR — Reduce memory AND relax threshold

```yaml
compairr_1mm:
  CompAIRRSequenceAbundance:
    compairr_path: /central/groups/MazmanianLab/jboktor/software/compairr/src/compairr
    p_value_threshold: [0.05, 0.01, 0.001]   # sweep; shallow depth needs relaxed
    ignore_genes: true                        # huge memory + time saving
    threads: 8
    sequence_batch_size: 5000                 # default 10000; reduce for LOSO large datasets
    keep_temporary_files: false
```

For LOSO CompAIRR (2800-3400 train samples), even this may OOM. Options:
- Raise SLURM `--mem=900GB`
- Or subsample the LOSO train set to ≤ 2000 samples (stratified by study and class) first

### 2D. Pubclone — Broaden hyperparameter grid

```yaml
seq_abundance:
  SequenceAbundance:
    p_value_threshold: [0.5, 0.1, 0.01]      # current 0.1 only; sweep
prob_binary:
  ProbabilisticBinaryClassifier:
    max_iterations: 1000                      # current 200
    update_rate: 0.01
```

Add `comparison_attributes: [amino_acid_sequence, v_gene, j_gene]` if V/J genes are in the AIRR files — this increases specificity.

### 2E. Assessment/Selection — Fix Statistical Power

**All configs currently use:**
```yaml
assessment: {split_strategy: random, split_count: 1, training_percentage: 0.7}
```
This gives a single point estimate, no CI, no way to detect overfitting.

**Recommended for within-study CV (before MLApplication to other studies):**
```yaml
assessment: {split_strategy: stratified_k_fold, split_count: 5}
selection:  {split_strategy: stratified_k_fold, split_count: 5}
```

**For transfer-study (PDBP→PPMI):** use the `MLApplication` pattern above. The `TrainMLModel` step does 5×5 nested CV on PDBP only; `MLApplication` applies the refit model to PPMI. Reports on PPMI become out-of-distribution performance.

---

## 3. Diagnostic Reports — Add Everywhere

| Report | What it tells us | Requirement |
|---|---|---|
| `ConfusionMatrix` | Is the model biased toward one class? | All models |
| `ROCCurve` / `ROCCurveSummary` | Per-fold ROC curves with mean ± SD | All models |
| `PrecisionRecallCurveSummary` | Better than ROC under class imbalance | All models |
| `TrainingPerformance` | Train vs test metrics — detects overfit | All models |
| `MLSettingsPerformance` | Compare settings side-by-side in one plot | All models |
| `Coefficients` (`n_largest: 100`) | Top features for interpretability | LogReg, SVM-linear, RF |
| `PerformancePerLabel` | Stratified by covariate — detect confounders | All models |
| `DeepRCMotifDiscovery` | Integrated gradients on PD-positive samples | DeepRC only |
| `KernelSequenceLogo` | Visualize CNN kernel motifs | DeepRC only |
| `DesignMatrixExporter` | Dump encoded feature matrix for external analysis | K-mer models |
| `SignificantKmerPositions` | Which CDR3 positions contribute | K-mer models |

**`PerformancePerLabel` alternative_labels to use** (all available in our metadata):

```yaml
perf_by_sex:
  PerformancePerLabel:
    alternative_label: sex
    metric: balanced_accuracy
perf_by_age_bin:
  PerformancePerLabel:
    alternative_label: age_at_baseline_binned
    metric: balanced_accuracy
perf_by_depth_bin:
  PerformancePerLabel:
    alternative_label: clonal_volume_binned   # ← critical: reveals depth-as-signal confound
    metric: balanced_accuracy
```

**Do NOT include `alternative_label: study`** in S2S configs — S2S trains within one study so there's no variation. Only use it in LOSO configs.

**One-time data characterization reports** (run once via `ExploratoryAnalysisInstruction`, not in every model config):

```yaml
label_dist: LabelDist
seq_count_dist: SequenceCountDistribution
shannon_diversity: ShannonDiversityOverview
clonotype_summary:
  RepertoireClonotypeSummary:
    color_label: case_control_other_latest
aa_freq_by_disease:
  AminoAcidFrequencyDistribution:
    label: case_control_other_latest
    split_by_label: true
    alignment: IMGT
    region_type: IMGT_CDR3
```

---

## 4. SLURM Strategy — CPU and GPU Stay SEPARATE

**The existing split in `shell_scripts/` is correct and must be preserved:**
- `immuneml_slurm_full_study_cpu.sh` → expansion partition, no GPU
- `immuneml_slurm_full_study_deeprc.sh` → gpu partition, 1 GPU

**Do not combine these into one array.** Each partition has different limits, each GPU tier (P100 vs H200) has different queue behavior, and mixing them in an array breaks scheduling.

### Proposed structure: Two independent arrays

**Array A — CPU models (expansion partition)**
- 80 jobs (5 CPU models × 16 study pairs)
- Generate from `immuneml_generate_S2S_yamls.py`
- Submit as `--array=0-79%30` (max 30 concurrent)
- Single SLURM script, reads yaml name from a manifest file indexed by `$SLURM_ARRAY_TASK_ID`

**Array B — DeepRC GPU (gpu partition, H200)**
- 13 jobs (DeepRC on 3 S2S-source studies × 3 targets + 4 LOSO; skip BioFIND-as-source)
- Submit as `--array=0-12%8` (max 8 concurrent H200s)
- Separate SLURM script with H200 `--gres` and 72h walltime

### CPU resource tiers (within Array A, use `--mem` per-job based on model)

Different CPU models have different memory profiles. Submit them as separate arrays per tier, not one monolithic array, to avoid over-allocating memory for fast jobs:

| Tier | Models | Partition | RAM | Time | CPUs | Array size |
|---|---|---|---|---|---|---|
| **Fast-CPU** | kmer_svm, kmer_rf, kmer_logreg | expansion | 96GB | 8h | 4 | 48 jobs (3 × 16 pairs) |
| **Heavy-CPU** | pubclone | expansion | 256GB | 24h | 8 | 16 jobs |
| **Heavy-CPU** | compairr | expansion | 900GB | 24h | 8 | 16 jobs |

### GPU tier (Array B)

| Tier | Model | Partition | GPU | RAM | Time | CPUs | Array size |
|---|---|---|---|---|---|---|---|
| **GPU** | DeepRC (spec-compliant) | gpu | h200 | 64GB | 48h | 8 | 13 jobs |

### Array job pattern

```bash
#!/bin/bash
#SBATCH --array=0-47%20
#SBATCH --partition=expansion
#SBATCH --mem=96GB
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4

MANIFEST=/path/to/cpu_jobs_manifest.txt
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $MANIFEST)
DESIGN=$(echo $LINE | cut -d' ' -f1)
CONFIG=$(echo $LINE | cut -d' ' -f2)

mamba activate immuneml
immune-ml /path/to/configs/$DESIGN/$CONFIG.yaml /path/to/results/$DESIGN/$CONFIG
```

Manifest file = plain text, one job per line: `S2S S2S_pdbp_to_ppmi_kmer_logreg`.

**Benefits over individual sbatch:**
- One cancel for the whole batch: `scancel <array_job_id>`
- One output dir pattern per array
- Automatic retry on failure with `--requeue`

---

## 5. Priority Ranked Actions

### P0 — Blocks further GPU spend
1. **Fix transfer-study design with `MLApplication`**. Current "S2S" results are within-study, not cross-study. Without this, no DeepRC retraining will produce meaningful transfer results.

### P1 — Statistical robustness
2. Switch all assessment splits to `stratified_k_fold` `split_count=5`.
3. Add `ConfusionMatrix`, `TrainingPerformance`, `PerformancePerLabel(clonal_volume_binned)` to every config.

### P2 — DeepRC retraining (after P0/P1)
4. Scale DeepRC params to **actual** repertoire depth: `sample_n_sequences: 500`, `n_updates: 100000`.
5. Skip DeepRC where N_train < 200 (i.e., skip BioFIND-as-source).

### P3 — Extend signal space
6. Add 4-mer, gapped k-mer encodings.
7. Add elastic net LogReg.
8. CompAIRR p-value sweep with `ignore_genes: true`.

### P4 — Infrastructure
9. Convert to two separate SLURM arrays (CPU array, GPU array).
10. Write centralized results aggregator Python script.
11. Run one-time `ExploratoryAnalysisInstruction` for diversity / motif overviews.

---

## 6. Interpretation Caveats Updated for This Dataset

- **Shallow depth ceiling:** at ~1000-4000 seqs/repertoire, no existing method is expected to achieve >0.65 bacc. Papers with 0.70+ used 10⁵-10⁶ seqs.
- **BioFIND-as-train wins may be a within-study overfit artifact** — recall the current bug treats "PDBP → PPMI" as "PDBP random 70/30 split." The strong BioFIND numbers may just reflect the model fitting to BioFIND's idiosyncratic clonal-volume distribution (median 1077, the lowest of any study). Once `MLApplication` is used for actual transfer, expect this to drop.
- **Expected true transfer bacc:** ~0.50-0.55 for k-mer models, 0.50 for DeepRC until retrained, 0.50 for CompAIRR/pubclone.

---

## Sources

- [DeepRC: NeurIPS 2020](https://proceedings.neurips.cc/paper/2020/hash/da4902cb0bc38210839714ebdcf0efc3-Abstract.html)
- [DeepRC GitHub](https://github.com/ml-jku/DeepRC)
- [CompAIRR paper](https://pubmed.ncbi.nlm.nih.gov/35852318/) — O(N²) distance matrix memory scaling
- [CompAIRR GitHub](https://github.com/uio-bmi/compairr) — `--no-matrix`, `--threads`
- [immuneML docs](https://docs.immuneml.uio.no/latest/) — reports, MLApplication, troubleshooting
- [Emerson et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28369033/) — original public-clone benchmark depth (~10⁵ seqs)
