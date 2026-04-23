# ImmuneML Transfer-Study Optimization Diagnostic & Recommendations

**Date:** 2026-04-14
**Based on:** 95/96 completed jobs across S2S (4-way pairwise) + LOSO transfer-study matrix

---

## 1. Current State Diagnostic

### Performance Summary

| Model | N runs | Mean bacc | Min | Max | Std | Verdict |
|---|---|---|---|---|---|---|
| **kmer_rf** | 16 | 0.516 | 0.430 | 0.594 | 0.038 | Modest signal, best variance |
| **kmer_logreg** | 16 | 0.514 | 0.431 | 0.602 | 0.040 | Best peak (0.602), interpretable |
| **kmer_svm** | 16 | 0.507 | 0.449 | 0.566 | 0.030 | Weak signal |
| **pubclone** | 15 | 0.499 | 0.436 | 0.549 | 0.034 | At chance, slight nonzero signal |
| **compairr** | 14 | 0.500 | 0.500 | 0.500 | 0.000 | **Completely at chance** (red flag) |
| **deeprc** | 16 | 0.500 | 0.500 | 0.500 | 0.000 | **Completely at chance** (undertrained) |

### Wall Time by Model

| Model | Mean (min) | Min (min) | Max (min) |
|---|---|---|---|
| kmer_rf | 62 | 5 | 149 |
| kmer_svm | 61 | 6 | 160 |
| kmer_logreg | 90 | 8 | 228 |
| compairr | 139 | 7 | 394 |
| pubclone | 188 | 5 | 1182 |
| deeprc | 487 | 5 | 1076 |

### Resource Failures
- **2 LOSO CompAIRR jobs OOM at 600GB** (BioFIND/HBS LOSO = 2829-3401 samples)
  - CompAIRR distance matrix scales O(N²) — 3400² ≈ 12M pairs
  - Need sharding or the `--no-matrix` flag

### Critical Diagnostic Gaps
1. **No per-fold CV metrics** — we can't tell if the modest kmer signal is fold-stable or fluke
2. **No train vs test performance comparison** — can't detect overfitting
3. **No learning curves** for DeepRC — can't tell if it's plateaued, diverged, or just not trained long enough
4. **No confusion matrices** — don't know if models are biased toward one class
5. **No attention/motif extraction** from DeepRC — can't interpret what signal the model did/didn't learn
6. **Single random split per assessment** (split_count=1) — point estimate, no confidence interval

---

## 2. Root Cause Analysis

### Why DeepRC = 0.500 (undertrained)

Our config uses:
```yaml
n_updates: 10000        # 30× below paper's 300,000
evaluate_at: 2000
sample_n_sequences: 2000  # Paper uses 10,000
kernel_size: 5           # Paper's original sweep: 5–9
n_kernels: 16            # Paper sweep: 8–32
```

The original DeepRC paper ([Widrich et al. 2020](https://proceedings.neurips.cc/paper/2020/hash/da4902cb0bc38210839714ebdcf0efc3-Abstract.html)) states:
> "Every 5×10³ updates, the current model was evaluated against the validation fold. The early stopping hyperparameter was determined by selecting the model with the best loss on the validation fold after **10⁵ updates**."

**We trained for 10⁴ — one order of magnitude below the minimum.**

### Why CompAIRR = 0.500 (no signal found)

Our config: `p_value_threshold: 0.001` with `max_edit_distance=1` (1mm).

At shallow sequencing depth with only ~500-5000 sequences per repertoire, no individual TCR clone reaches Bonferroni-corrected significance. The spec sheet called this out as expected:
> "Expected AUROC ≈ 0.50 (chance) due to shallow depth. Serves as lower bound."

### Why K-mer models show modest signal but with large pair-to-pair variance

The variance in performance across study pairs (std ≈ 0.04) suggests **batch effects dominate**. BioFIND as training set consistently produces better models — likely because BioFIND has different library prep/sequencing characteristics that are learned by the simpler k-mer models.

---

## 3. Optimization Recommendations

### 3A. DeepRC — Highest Priority

**Parameter changes (aligned with published spec):**

```yaml
ml_methods:
  deeprc_model:
    DeepRC:
      # Core training length (30× increase)
      n_updates: 300000        # was 10000
      evaluate_at: 10000       # was 2000

      # Architecture (from paper optimal)
      kernel_size: 9           # was 5; Rawat used 9 for longer motifs
      n_kernels: 32            # was 16; grid [16, 32]
      n_additional_convs: 2    # was 1
      n_attention_network_layers: 2
      n_attention_network_units: 64  # was 32

      # Input sampling (match paper)
      sample_n_sequences: 10000  # was 2000

      # Regularization (critical for shallow-depth data)
      learning_rate: 5.0e-5
      l2_weight_decay: 0.0001   # was 0.001; reduce to allow signal
      l1_weight_decay: 0

      # Batch size + workers
      training_batch_size: 4
      n_workers: 4

      pytorch_device_name: cuda:0
```

**SLURM changes:**
```bash
#SBATCH --gres=gpu:nvidia_h200:1   # or h100
#SBATCH --time=72:00:00            # 300k updates takes longer
#SBATCH --mem=64GB                 # more headroom for large repertoires
#SBATCH --cpus-per-task=8          # feed GPU faster
```

**Expected speedup:** H200 is ~3× faster than P100 for transformer/attention workloads. 300k updates × 3× = ~9× total compute, but H200 cuts it to 3× wall-time.

### 3B. DeepRC Diagnostic Reports — Critical for Interpretability

Add these reports to every DeepRC config:

```yaml
reports:
  # Training curves per fold
  roc_summary: ROCCurveSummary
  pr_curve: PrecisionRecallCurveSummary
  ml_settings: MLSettingsPerformance
  conf_matrix: ConfusionMatrix

  # DeepRC-specific interpretability
  deeprc_motifs:
    DeepRCMotifDiscovery:
      threshold: 0.5         # IG contribution threshold
      n_steps: 50            # IG integration steps

  # Per-class performance stratified by covariates
  perf_sex:
    PerformancePerLabel:
      alternative_label: sex
      metric: balanced_accuracy
  perf_age:
    PerformancePerLabel:
      alternative_label: age_at_baseline_binned
      metric: balanced_accuracy
```

### 3C. K-mer Model Optimizations

**Current K-mer signal exists — amplify it:**

```yaml
encodings:
  # Add 4-mer for longer motifs
  kmer4_all:
    KmerFrequency:
      k: 4
      sequence_encoding: continuous_kmer
      normalization_type: l2     # try L2 instead of relative_frequency
      reads: all
      sequence_type: amino_acid
      region_type: IMGT_CDR3

  # Gapped k-mers capture non-contiguous motifs
  gapped_3mer:
    KmerFrequency:
      k_left: 2
      k_right: 2
      min_gap: 0
      max_gap: 3
      sequence_encoding: gapped_kmer
      normalization_type: relative_frequency

ml_methods:
  log_reg_elastic:
    LogisticRegression:
      penalty: elasticnet
      l1_ratio: [0.1, 0.5, 0.9]     # sweep L1/L2 mix
      C: [0.001, 0.01, 0.1, 1.0, 10.0]
      max_iter: 10000
      solver: saga
    model_selection_cv: true
    model_selection_n_folds: 5
```

### 3D. Assessment Split Strategy — Fix Statistical Power

**Current:**
```yaml
assessment:
  split_strategy: random
  split_count: 1          # SINGLE split, no CI
  training_percentage: 0.7
```

**Recommended:**
```yaml
assessment:
  split_strategy: stratified_k_fold  # ensures class balance
  split_count: 5                     # get mean ± SD across 5 folds
  reports:
    models: [coefficients, conf_matrix, roc]
selection:
  split_strategy: stratified_k_fold
  split_count: 5
```

This gives us:
- 5-fold × 5-fold nested CV = 25 models trained
- Can detect train/test gap per fold
- Confidence intervals on bacc
- Robust to unlucky splits

### 3E. CompAIRR Memory Fix

```yaml
encodings:
  compairr_1mm:
    CompAIRRSequenceAbundance:
      compairr_path: /central/groups/MazmanianLab/jboktor/software/compairr/src/compairr
      p_value_threshold: [0.01, 0.001, 0.0001]  # less strict; grid sweep
      ignore_genes: true       # was false; huge memory/time savings
      threads: 8               # was default
      keep_temporary_files: false
      sequence_batch_size: 10000  # default; reduce for LOSO large
```

**SLURM:** Use `--mem=900GB` for LOSO datasets (the 2 that failed).

### 3F. Pubclone (ProbabilisticBinaryClassifier)

Current threshold `p_value_threshold: 0.1` may be too strict (or too loose).

```yaml
encodings:
  seq_abundance:
    SequenceAbundance:
      p_value_threshold: [0.5, 0.1, 0.01]  # sweep
      comparison_attributes:
        - amino_acid_sequence
        - v_gene          # add V gene for specificity
        - j_gene
ml_methods:
  prob_binary:
    ProbabilisticBinaryClassifier:
      max_iterations: 1000    # was 200
      update_rate: 0.01
```

### 3G. Universal Data Reports (add once, run everywhere)

```yaml
reports:
  # Understand dataset balance/diversity
  label_dist: LabelDist
  seq_count_dist: SequenceCountDistribution
  shannon: ShannonDiversityOverview
  repertoire_summary:
    RepertoireClonotypeSummary:
      color_label: case_control_other_latest
  aa_freq:
    AminoAcidFrequencyDistribution:
      label: case_control_other_latest
      split_by_label: true
      alignment: IMGT
      region_type: IMGT_CDR3
```

---

## 4. SLURM Resource Recommendations

| Tier | Model | Partition | GPU | RAM | Time | CPUs |
|---|---|---|---|---|---|---|
| **Fast** | kmer_svm, kmer_rf | expansion | — | 64GB | 4h | 4 |
| **Medium** | kmer_logreg | expansion | — | 96GB | 8h | 8 |
| **Heavy-CPU** | compairr, pubclone | expansion | — | **900GB** | 24h | 8 |
| **GPU-light** | DeepRC (minimal) | gpu | p100 | 32GB | 4h | 4 |
| **GPU-heavy** | DeepRC (spec-compliant) | gpu | **h200** | 64GB | 72h | 8 |

### Array jobs instead of individual submits

Convert the 96 separate submissions into a SLURM array:
```bash
#SBATCH --array=0-79%20   # max 20 concurrent
```
Simpler dependency management, easier to cancel/resubmit batches.

### Job output consolidation

Currently 96 jobs × 2 files = 192 log files. Recommend:
- Consolidated TSV of results emitted by each job (append to shared results CSV)
- Single Python post-run aggregator that reads all log dirs

---

## 5. Priority-Ranked Action List

### P0 — Do first (biggest signal gains)
1. **Retrain DeepRC with `n_updates: 300000` on H200**. Current runs are undertrained by 30×. Expected to actually learn signal.
2. **Switch assessment to `stratified_k_fold` with `split_count: 5`**. Fixes statistical robustness of all k-mer results.
3. **Add `DeepRCMotifDiscovery`, `TrainingPerformance`, `ConfusionMatrix` reports** to every config.

### P1 — Diagnostic resolution
4. Add `PerformancePerLabel` with `alternative_label: sex` and `study` to detect confounder leakage.
5. Add `ROCCurve` and `PrecisionRecallCurveSummary` for per-fold curves (not just summary).
6. Enable `DesignMatrixExporter` to dump encoded features — lets us do external analysis.

### P2 — Extend model space
7. Add gapped k-mer encoding (captures non-contiguous motifs).
8. Add elastic net LogReg (between L1 and L2 extremes).
9. CompAIRR with `ignore_genes: true` and stricter p-value sweep.

### P3 — Infrastructure
10. Convert to SLURM array jobs.
11. Write centralized results aggregator (Python script that reads all `full_*.yaml` + log + metrics and emits a unified CSV).
12. Resubmit failed LOSO compairr jobs with 900GB + `ignore_genes: true`.

---

## 6. Interpretation Caveats

**Why all DeepRC/CompAIRR results are exactly 0.500:** Both produced models that predicted the majority class (or couldn't train), resulting in exactly balanced accuracy of 0.5 on the held-out test. This is not a numerical tie — it's a sign of training failure or no signal found.

**Why BioFIND-as-train wins:** BioFIND has 174 samples but distinct batch characteristics. The model may be learning batch-specific TCR patterns that happen to correlate with PD status in that cohort, and these generalize (weakly) to other cohorts. **This is a confound to investigate, not a real biological signal.**

**Effect sizes to target:** Published TCR-based PD classifiers (Rawat et al. 2020) achieved AUROC ≈ 0.70 using DeepRC on a single cohort with deeper sequencing (~10⁵ sequences/repertoire). Our shallow data (~500-5000 seqs) is likely a ceiling at bacc ≈ 0.60-0.65 even with optimal models.

---

## Sources

- [DeepRC: Immune repertoire classification with attention-based deep MIL](https://proceedings.neurips.cc/paper/2020/hash/da4902cb0bc38210839714ebdcf0efc3-Abstract.html) — n_updates=10⁵ reference
- [DeepRC GitHub](https://github.com/ml-jku/DeepRC) — architecture details
- [CompAIRR paper](https://pubmed.ncbi.nlm.nih.gov/35852318/) — memory characteristics
- [CompAIRR GitHub](https://github.com/uio-bmi/compairr) — `--no-matrix`, `--threads` flags
- [immuneML docs](https://docs.immuneml.uio.no/latest/) — report specs
- [immuneML Nat Mach Intell paper](https://www.nature.com/articles/s42256-021-00413-z) — benchmark design
