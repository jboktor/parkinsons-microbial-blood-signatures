# DeepRC Hard Judgment Calls — Setting Up for Success

**Date:** 2026-04-14
**Goal:** Maximize probability DeepRC finds real signal given our shallow-depth, multi-cohort dataset.

---

## The Fundamental Tradeoff

DeepRC's attention-based multiple instance learning (MIL) needs:
- **Depth per bag:** enough sequences per repertoire for attention to isolate signal from noise
- **Number of bags:** enough training samples for the CNN+attention weights to generalize

Our constraint: we have neither. Current `clonal_volume ≥ 500` filter + small cohorts = shallow bags × few samples.

Four levers to pull (in order of impact):

1. **Depth threshold elevation** — trade sample count for per-bag resolution
2. **Dataset pooling** — trade study purity for training size
3. **Architecture + training regime** — right-size the model for the data
4. **Ensembling** — reduce MIL training noise

---

## Judgment Call 1: Raise depth threshold to 1000 for DeepRC only

**Current:** filter floor 500 sequences/repertoire.

**Quantified tradeoff:**

| Study | N@≥500 | N@≥1000 | N@≥1500 | N@≥2000 |
|---|---|---|---|---|
| BioFIND | 174 | **99 (-43%)** | 36 (-79%) | 14 (-92%) |
| PDBP | 1234 | 947 (-23%) | 712 (-42%) | 603 (-51%) |
| PPMI | 1421 | 1271 (-11%) | 1094 (-23%) | 960 (-32%) |
| HBS | 746 | 745 (-0%) | 739 (-1%) | 718 (-4%) |

**Judgment:** Raise the DeepRC-only threshold to **1000**. Keep k-mer models at 500.

Why 1000 not 1500:
- At 1000, we lose 43% of BioFIND (now 99 samples) but only 11-23% of PDBP/PPMI.
- At 1500, BioFIND drops to 36 — too small to use at all.
- 1000 sequences/bag × attention pooling gives roughly the signal resolution DeepRC needs.

Apply via metadata filter or via ImmuneML preprocessing step. Simplest: write a filtered metadata CSV for DeepRC only.

**Do NOT use BioFIND as a DeepRC training source at any threshold.** 99 samples × 70/30 split × nested CV = ~55 train / 14 val / 14 test per fold. Overfitting is mathematically guaranteed.

---

## Judgment Call 2: Skip single-study S2S training for DeepRC entirely

S2S training on one cohort gives at most ~1400 samples (PPMI) after threshold filter. LOSO gives 1800-3000. For MIL, more bags >> bigger bags once you're past the attention-resolvability floor.

**Recommended DeepRC training regimes (in priority order):**

| Regime | Train set | N@≥1000 | Viability |
|---|---|---|---|
| **LOSO held-out BioFIND** | PDBP+PPMI+HBS | 2963 | **Best** |
| **LOSO held-out PDBP** | PPMI+BioFIND+HBS | 2115 | Good |
| **LOSO held-out HBS** | PDBP+PPMI+BioFIND | 2317 | Good |
| **LOSO held-out PPMI** | PDBP+BioFIND+HBS | 1791 | Good |
| S2S PPMI-train (alt) | PPMI only | 1271 | Marginal |
| S2S PDBP-train (alt) | PDBP only | 947 | Marginal |
| ~~S2S HBS-train~~ | HBS only | 745 | Skip |
| ~~S2S BioFIND-train~~ | BioFIND only | 99 | Skip |

**Judgment:** Run **4 LOSO configs + 2 S2S-from-PPMI / S2S-from-PDBP as sanity checks = 6 total DeepRC runs** instead of 16. Drop BioFIND-source and HBS-source. The compute savings (~60%) goes into longer training per run.

---

## Judgment Call 3: Match DeepRC architecture to data scale

The Widrich et al. paper tested on 100-1000 samples × 10⁵ seqs/bag. We have 1000-3000 samples × ~1000 seqs/bag. The bag/sample balance is flipped — we have *relatively* more bags and shallower bags.

**Architecture for shallow bags + moderate samples:**

```yaml
deeprc_model:
  DeepRC:
    # Architecture: smaller CNN, stronger attention
    kernel_size: 7                  # was 5 — CDR3β medians 14 AAs, kernel=7 captures midrange motifs
    n_kernels: 32                   # was 16 — with 2000+ training bags we can afford this
    n_additional_convs: 1           # keep shallow — we need translation invariance, not depth
    n_attention_network_layers: 2
    n_attention_network_units: 64   # was 32 — bigger attention helps pool shallow bags

    # Sequence handling
    sample_n_sequences: 1000        # matches the new depth floor
    add_positional_information: true # CRITICAL — tells CNN which position in CDR3 each residue is
    consider_seq_counts: true        # CRITICAL — weight by clone count, not uniform
    sequence_counts_scaling_fn: log  # log-scales counts to prevent mega-clones dominating
    sequence_reduction_fraction: 0.1 # attention pre-selects top 10% before gradient step
    reduction_mb_size: 100

    # Training
    n_updates: 200000               # 2× the paper minimum; with early stopping on val loss
    evaluate_at: 5000
    training_batch_size: 16         # was 4 — H200 can handle it, stabler gradients
    learning_rate: 5.0e-5           # paper default
    n_torch_threads: 4
    n_workers: 4

    # Regularization — strong, because shallow bags = noisy attention
    l2_weight_decay: 0.001          # keep
    l1_weight_decay: 1.0e-5         # add small L1 for sparsity

    # Compute
    keep_dataset_in_ram: true        # LOSO datasets fit in 64GB; cuts I/O to zero
    pytorch_device_name: cuda:0
```

**Key additions not in our current configs:**

1. **`add_positional_information: true`** — The paper Table A4 shows this improves AUROC by ~0.03. Free win.

2. **`consider_seq_counts: true` + `sequence_counts_scaling_fn: log`** — Currently we sample uniformly across sequences. With log-weighting, expanded clones get more gradient signal (they're more likely antigen-specific) without mega-clones drowning out rare disease signatures.

3. **`sequence_reduction_fraction: 0.1`** — During each training update, attention pre-selects the top 10% of the 1000 sampled sequences for the gradient step. This is the *same* mechanism the attention layer uses at inference, applied during training. It should speed convergence and reduce noise.

4. **`training_batch_size: 16`** — Paper used 4 because their GPUs were smaller. H200 has 141GB VRAM; we can easily fit 16-32 bags. Larger batch = less noisy gradients.

5. **`keep_dataset_in_ram: true`** — LOSO datasets are ~2000-3000 bags × 1000 seqs × (few bytes per seq token) < 2GB. Easy to hold in RAM, eliminates I/O bottleneck.

---

## Judgment Call 4: Train longer with smart early stopping

**Current:** `n_updates: 10000, evaluate_at: 2000` → 5 checkpoints, no time to converge.

**Paper:** `n_updates: 300000, evaluate_at: 5000` → 60 checkpoints with best-val selection.

**Our proposal:** `n_updates: 200000, evaluate_at: 5000` → 40 checkpoints.

The paper's `300000` was for repertoire sizes 100× ours. With shallower bags, convergence is faster (less signal means gradient saturates sooner). 200k is a defensible middle ground.

**DeepRC's built-in training loop** saves best-val-loss checkpoint every `evaluate_at` updates. Effective early stopping is automatic — we just need to budget enough compute for the validation curve to plateau.

---

## Judgment Call 5: Ensemble 5 seeds → 1 prediction

MIL training is noisy. Different random seeds produce models with ~0.05 AUC variance on the same data. Ensemble averaging cuts that variance by √5.

**Implementation:** Run each DeepRC config 5× with different `--seed` and average probabilities at prediction time. ImmuneML doesn't natively support seed ensembling — we'd do it in a small post-hoc Python script that reads the 5 prediction CSVs and averages.

**Cost:** 5× compute per config. With the reduced run count (6 configs instead of 16), total DeepRC compute goes from 16 × 1 = 16 runs → 6 × 5 = 30 runs. Net ~2× compute, probably 2-3× better SNR.

**If compute is tight:** ensemble only the 4 LOSO configs (most important for transfer claims). Skip ensembling the S2S sanity checks.

---

## Judgment Call 6: Feed DeepRC the best-quality sequences, not uniform samples

Right now `sample_n_sequences: N` picks N sequences uniformly from each bag. But not all sequences are equal:
- **Expanded clones** (high `duplicate_count`) are more likely antigen-specific
- **Rare singleton clones** are mostly noise
- **Sequences flagged `productive: true`** are already our filter, but we could further filter by clone count

**Two ways to bias toward expanded clones:**

A. **`consider_seq_counts: true` + `sequence_counts_scaling_fn: log`** (built into DeepRC, mentioned above). Weights uniform sampling by log(count).

B. **Preprocessing: keep top 1000 by `duplicate_count`**, not a random 1000. This is a one-time filter we'd apply before ImmuneML imports the data. Harder to wire up but gives a cleaner prior.

**Recommendation:** Use (A) first — it's a config flag. Only add (B) if (A) doesn't produce signal after the threshold + architecture + ensemble changes.

---

## Judgment Call 7: Pre-train on healthy-vs-generated-background (moonshot)

DeepRC's bottleneck is "signal is rare in noisy bags." Traditional fix: more labeled data. Alternative fix: **pre-train on a different, easier task and transfer the learned kernels.**

Train DeepRC to distinguish **real human TCR repertoires** from **generated-from-OLGA-model background sequences** (which lack disease selection pressure). This task has essentially unlimited training data, the model learns "what a real TCR looks like" in its kernels, and those kernels transfer to downstream disease classification.

**Feasibility:** ImmuneML supports OLGA/SONIA via `sonnia` dependency (already installed). We'd need to write a custom pipeline — not a simple YAML change. High payoff if our 4-cohort data can't overcome the MIL small-sample problem.

**Recommendation:** Defer unless our P0-P2 recommendations still give 0.50 bacc.

---

## Summary — Combined DeepRC Recommended Regime

**Pre-run:**
1. Create `metadata_deeprc_filtered.csv` with `clonal_volume ≥ 1000` (loses 43% BioFIND, 23% PDBP, 11% PPMI, 0% HBS).
2. Generate LOSO metadata files from the filtered metadata.
3. Skip S2S-from-BioFIND and S2S-from-HBS entirely for DeepRC.

**Run config:** `n_updates=200000`, `evaluate_at=5000`, `sample_n_sequences=1000`, `kernel_size=7`, `n_kernels=32`, `n_attention_network_units=64`, `add_positional_information=true`, `consider_seq_counts=true`, `sequence_counts_scaling_fn=log`, `sequence_reduction_fraction=0.1`, `training_batch_size=16`, `keep_dataset_in_ram=true`, `l2_weight_decay=0.001`, `l1_weight_decay=1e-5`.

**SLURM:** H200, 48-72h walltime, 64GB RAM, 8 CPUs.

**Run plan:** 6 configs × 5 seeds = 30 runs, submitted as a **separate GPU-only SLURM array**.

**Expected outcome:** If there's a real TCR-based PD signal accessible at this depth, this regime will find it. If we still see 0.50 across the board after this, it's strong evidence that **~1000 seqs/repertoire is below the detection floor for TCR-based PD classification** — a legitimate negative result with publishable value.

---

## What I'm NOT Recommending (and why)

- **Do not** increase model size further (n_kernels > 32, n_attention_units > 64). Our sample count doesn't support it; overfit risk.
- **Do not** mix BioFIND into training pools. Its characteristics differ enough to introduce batch bias without adding much signal.
- **Do not** run all 16 original S2S pairs. 10 of them are starvation-doomed.
- **Do not** raise threshold to 1500+. Losing 79% of BioFIND makes it useless as a test target, and PDBP loses 42% without proportional signal gain.
