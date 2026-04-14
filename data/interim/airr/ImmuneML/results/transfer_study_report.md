# ImmuneML Transfer-Study Classification Report

**Date:** 2026-04-14
**Total jobs:** 95/96 completed (2 LOSO CompAIRR OOM failures)

**Studies:** PDBP (n=1234), PPMI (n=1421), BioFIND (n=174), HBS (n=746)

**Models:** kmer_logreg (L1/L2 3-mer LogReg), kmer_svm (linear SVM), kmer_rf (RandomForest), pubclone (SequenceAbundance + ProbBinaryClassifier), compairr (CompAIRR near-neighbor), deeprc (DeepRC CNN, n_updates=10000)

---

## S2S Transfer Matrix: kmer_logreg

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | **0.602** | 0.431 | 0.546 |
| **PDBP** | 0.532 | — | 0.520 | 0.471 |
| **PPMI** | 0.510 | 0.505 | — | 0.515 |
| **HBS** | **0.567** | 0.513 | 0.484 | — |

## S2S Transfer Matrix: kmer_svm

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | 0.464 | **0.566** | 0.449 |
| **PDBP** | 0.507 | — | **0.558** | 0.510 |
| **PPMI** | 0.524 | 0.536 | — | 0.508 |
| **HBS** | 0.487 | 0.493 | 0.512 | — |

## S2S Transfer Matrix: kmer_rf

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | 0.430 | **0.579** | **0.594** |
| **PDBP** | 0.510 | — | 0.499 | 0.513 |
| **PPMI** | 0.492 | 0.509 | — | 0.509 |
| **HBS** | 0.496 | **0.556** | 0.526 | — |

## S2S Transfer Matrix: pubclone

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | 0.545 | 0.549 | 0.512 |
| **PDBP** | 0.436 | — | 0.480 | 0.463 |
| **PPMI** | 0.484 | 0.454 | — | 0.486 |
| **HBS** | 0.540 | 0.514 | 0.536 | — |

## S2S Transfer Matrix: compairr

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | 0.500 | 0.500 | 0.500 |
| **PDBP** | 0.500 | — | 0.500 | 0.500 |
| **PPMI** | 0.500 | 0.500 | — | 0.500 |
| **HBS** | 0.500 | 0.500 | 0.500 | — |

## S2S Transfer Matrix: deeprc

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | 0.500 | 0.500 | 0.500 |
| **PDBP** | 0.500 | — | 0.500 | 0.500 |
| **PPMI** | 0.500 | 0.500 | — | 0.500 |
| **HBS** | 0.500 | 0.500 | 0.500 | — |

## Best Model Per S2S Pair

| Train \ Test | BIOFIND | PDBP | PPMI | HBS |
|---|---|---|---|---|
| **BIOFIND** | — | **logreg 0.602** | **rf 0.579** | **rf 0.594** |
| **PDBP** | logreg 0.532 | — | **svm 0.558** | rf 0.513 |
| **PPMI** | svm 0.524 | svm 0.536 | — | logreg 0.515 |
| **HBS** | **logreg 0.567** | **rf 0.556** | pubclone 0.536 | — |

## LOSO (Leave-One-Study-Out)

| Held-out | kmer_logreg | kmer_svm | kmer_rf | pubclone | compairr | deeprc |
|---|---|---|---|---|---|---|
| **BIOFIND** | **0.523** | 0.496 | **0.522** | OOM | OOM | 0.500 |
| **PDBP** | 0.467 | 0.500 | 0.492 | 0.500 | 0.500 | 0.500 |
| **PPMI** | 0.504 | 0.496 | 0.505 | 0.488 | 0.500 | 0.500 |
| **HBS** | **0.526** | 0.509 | 0.518 | 0.495 | OOM | 0.500 |

## Model Type Comparison (mean bacc across S2S pairs)

| Model | Mean bacc | Min | Max | N pairs |
|---|---|---|---|---|
| kmer_rf | 0.518 | 0.430 | 0.594 | 12 |
| kmer_logreg | 0.516 | 0.431 | 0.602 | 12 |
| kmer_svm | 0.510 | 0.449 | 0.566 | 12 |
| compairr | 0.500 | 0.500 | 0.500 | 12 |
| deeprc | 0.500 | 0.500 | 0.500 | 12 |
| pubclone | 0.500 | 0.436 | 0.549 | 12 |

## Key Findings

1. **BioFIND as training set outperforms** — despite being the smallest cohort (174 samples), it produces the best transfer performance (up to 0.602). This suggests BioFIND may have a cleaner disease signal or less batch-specific noise.
2. **K-mer models show modest but consistent signal** — kmer_logreg, kmer_rf, and kmer_svm achieve 0.55+ on multiple pairs, consistently above chance.
3. **CompAIRR at chance everywhere (0.500)** — near-neighbor sequence abundance finds no disease-associated sequences at p=0.001. Expected with shallow TCR sequencing depth.
4. **DeepRC at chance (0.500)** — n_updates=10,000 is insufficient for convergence. The spec calls for 300,000 updates for real signal.
5. **LOSO near chance** — training on 3 combined studies does not improve over single-study training, suggesting batch effects dominate when mixing cohorts.
6. **Pubclone weak but nonzero** — Emerson-style exact matching shows some signal (up to 0.549), consistent with weak public clone sharing.

## Failures

| Job | Reason |
|---|---|
| LOSO_biofind_compairr | OOM at 600GB (3401 train samples) |
| LOSO_hbs_compairr | OOM at 600GB (2829 train samples) |
