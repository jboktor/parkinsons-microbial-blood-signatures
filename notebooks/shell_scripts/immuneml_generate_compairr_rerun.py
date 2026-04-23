#!/usr/bin/env python3
"""Regenerate CompAIRR configs as combined TrainMLModel+MLApplication (single-run)."""

import os
from pathlib import Path
from itertools import permutations

BASE = Path("/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML")
CONFIGS = BASE / "configs"
META = BASE / "metadata"
REPO = BASE / "data_input" / "participant_repertoires"

STUDIES = ["pdbp", "ppmi", "biofind", "hbs"]

TEMPLATE = """definitions:
  datasets:
    train_dataset:
      format: AIRR
      params:
        path: {repo}
        metadata_file: {train_meta}
        is_repertoire: true
        import_productive: true
        import_with_stop_codon: false
        import_out_of_frame: false
        import_illegal_characters: false
        import_empty_nt_sequences: true
        import_empty_aa_sequences: false
        region_type: IMGT_CDR3
    test_dataset:
      format: AIRR
      params:
        path: {repo}
        metadata_file: {test_meta}
        is_repertoire: true
        import_productive: true
        import_with_stop_codon: false
        import_out_of_frame: false
        import_illegal_characters: false
        import_empty_nt_sequences: true
        import_empty_aa_sequences: false
        region_type: IMGT_CDR3

  encodings:
    compairr_1mm:
      CompAIRRSequenceAbundance:
        compairr_path: /central/groups/MazmanianLab/jboktor/software/compairr/src/compairr
        p_value_threshold: 0.001
        ignore_genes: true
        threads: 8

  ml_methods:
    log_reg_compairr:
      LogisticRegression:
        penalty: l1
        C: [0.1, 1.0]
      model_selection_cv: true
      model_selection_n_folds: 3

  reports:
    roc_summary: ROCCurveSummary

instructions:
  {train_name}:
    type: TrainMLModel
    dataset: train_dataset
    labels:
      - case_control_other_latest:
          positive_class: Case
    settings:
      - encoding: compairr_1mm
        ml_method: log_reg_compairr
    assessment:
      split_strategy: random
      split_count: 1
      training_percentage: 0.7
    selection:
      split_strategy: k_fold
      split_count: 3
    strategy: GridSearch
    optimization_metric: balanced_accuracy
    metrics: [auc, balanced_accuracy]
    reports:
      - roc_summary
    refit_optimal_model: true
    number_of_processes: 8

  {apply_name}:
    type: MLApplication
    dataset: test_dataset
    config_path: {zip_ref}
    number_of_processes: 4
    metrics:
      - accuracy
      - balanced_accuracy
      - auc
      - precision
      - recall
      - f1_micro
      - confusion_matrix
"""


def make(train_study, test_study, outdir, out_name, train_meta=None, test_meta=None):
    train_meta = train_meta or (META / f"metadata_pid_{train_study}.csv")
    test_meta = test_meta or (META / f"metadata_pid_{test_study}.csv")
    train_name = f"train_{train_study}_compairr"
    apply_name = f"apply_{train_study}_to_{test_study}_compairr"
    # In ImmuneML the MLApplication's config_path reference the zip output of prior TrainMLModel
    # in the same run. Using a relative path based on the standard output layout:
    zip_ref = f"./{train_name}/optimal_case_control_other_latest/zip/ml_settings_case_control_other_latest.zip"
    yaml_str = TEMPLATE.format(
        repo=REPO, train_meta=train_meta, test_meta=test_meta,
        train_name=train_name, apply_name=apply_name, zip_ref=zip_ref)
    with open(outdir / f"{out_name}.yaml", 'w') as f:
        f.write(yaml_str)


def main():
    outdir = CONFIGS / "compairr_rerun"
    outdir.mkdir(exist_ok=True)
    count = 0

    # S2S: 12 pairs
    for train_study, test_study in permutations(STUDIES, 2):
        make(train_study, test_study, outdir,
             f"compairr_{train_study}_to_{test_study}")
        count += 1

    # LOSO: 4 pairs
    for study in STUDIES:
        train_meta = META / f"metadata_pid_{study}_loso.csv"
        test_meta = META / f"metadata_pid_{study}.csv"
        make(f"not_{study}", study, outdir, f"compairr_LOSO_{study}",
             train_meta=train_meta, test_meta=test_meta)
        count += 1

    print(f"Generated {count} combined CompAIRR configs in {outdir}")


if __name__ == "__main__":
    main()
