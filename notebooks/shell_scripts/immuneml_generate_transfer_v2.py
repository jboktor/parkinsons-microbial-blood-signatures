#!/usr/bin/env python3
"""Generate v2 transfer-study configs with proper MLApplication + rich diagnostics."""

import os
from pathlib import Path
from itertools import permutations

BASE = Path("/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML")
CONFIGS = BASE / "configs"
META = BASE / "metadata"
META_FILT = META / "deeprc_filtered"
RESULTS = BASE / "results"
REPO = BASE / "data_input" / "participant_repertoires"

STUDIES = ["pdbp", "ppmi", "biofind", "hbs"]

# =========================================================================
# MLApplication-only YAMLs (to apply existing trained models to new studies)
# =========================================================================

APPLY_TEMPLATE = """definitions:
  datasets:
    apply_dataset:
      format: AIRR
      params:
        path: {repo}
        metadata_file: {meta}
        is_repertoire: true
        import_productive: true
        import_with_stop_codon: false
        import_out_of_frame: false
        import_illegal_characters: false
        import_empty_nt_sequences: true
        import_empty_aa_sequences: false
        region_type: IMGT_CDR3

instructions:
  {name}:
    type: MLApplication
    dataset: apply_dataset
    config_path: {zip_path}
    number_of_processes: 4
    metrics:
      - accuracy
      - balanced_accuracy
      - auc
      - precision
      - recall
      - f1_micro
      - f1_macro
      - confusion_matrix
"""


def make_apply_configs():
    """For each already-trained model, create a YAML that applies it to the target study."""
    out_dir = CONFIGS / "apply_v2"
    out_dir.mkdir(exist_ok=True)

    models = ["kmer_logreg", "kmer_svm", "kmer_rf", "pubclone", "compairr"]
    count = 0

    # S2S: apply each trained model to its target study (not the source)
    for train_study, test_study in permutations(STUDIES, 2):
        for model in models:
            # Find the existing zip from the S2S run
            src_instr_name = f"S2S_{train_study}_to_{test_study}_{model}"
            zip_path = RESULTS / "S2S" / src_instr_name / src_instr_name / "optimal_case_control_other_latest" / "zip" / "ml_settings_case_control_other_latest.zip"
            if not zip_path.exists():
                continue
            # The zip was built while the trainer used train_study data; now apply to REAL test_study
            test_meta = META / f"metadata_pid_{test_study}.csv"
            inst_name = f"apply_{train_study}_to_{test_study}_{model}"
            yaml_str = APPLY_TEMPLATE.format(
                repo=REPO, meta=test_meta, name=inst_name, zip_path=zip_path)
            with open(out_dir / f"{inst_name}.yaml", 'w') as f:
                f.write(yaml_str)
            count += 1

    # LOSO: apply each trained model to its held-out study
    for study in STUDIES:
        for model in models:
            src_instr_name = f"LOSO_{study}_{model}"
            # LOSO instruction name uses "not_X_to_X" internally
            inner_instr = f"LOSO_not_{study}_to_{study}_{model}"
            zip_path = RESULTS / "LOSO" / src_instr_name / inner_instr / "optimal_case_control_other_latest" / "zip" / "ml_settings_case_control_other_latest.zip"
            if not zip_path.exists():
                continue
            test_meta = META / f"metadata_pid_{study}.csv"
            inst_name = f"apply_LOSO_{study}_{model}"
            yaml_str = APPLY_TEMPLATE.format(
                repo=REPO, meta=test_meta, name=inst_name, zip_path=zip_path)
            with open(out_dir / f"{inst_name}.yaml", 'w') as f:
                f.write(yaml_str)
            count += 1

    return count


# =========================================================================
# New DeepRC configs — optimized per deeprc_strategy.md
# =========================================================================

DEEPRC_TEMPLATE = """definitions:
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

  encodings:
    deeprc_enc: DeepRC

  ml_methods:
    deeprc_model:
      DeepRC:
        validation_part: 0.2
        kernel_size: 7
        n_kernels: 32
        n_additional_convs: 1
        n_attention_network_layers: 2
        n_attention_network_units: 64
        n_output_network_units: 32
        learning_rate: 5.0e-5
        n_updates: 200000
        evaluate_at: 5000
        sample_n_sequences: 1000
        add_positional_information: true
        consider_seq_counts: true
        sequence_counts_scaling_fn: log
        sequence_reduction_fraction: 0.1
        reduction_mb_size: 100
        training_batch_size: 16
        n_workers: 4
        n_torch_threads: 4
        l2_weight_decay: 0.001
        l1_weight_decay: 1.0e-5
        keep_dataset_in_ram: true
        pytorch_device_name: cuda:0

  reports:
    roc_summary: ROCCurveSummary
    pr_curve: PrecisionRecallCurveSummary
    ml_settings: MLSettingsPerformance
    conf_matrix: ConfusionMatrix
    perf_by_depth:
      PerformancePerLabel:
        alternative_label: clonal_volume_binned
        metric: balanced_accuracy
    perf_by_sex:
      PerformancePerLabel:
        alternative_label: sex
        metric: balanced_accuracy
    perf_by_age:
      PerformancePerLabel:
        alternative_label: age_at_baseline_binned
        metric: balanced_accuracy

instructions:
  {name}:
    type: TrainMLModel
    dataset: train_dataset
    labels:
      - case_control_other_latest:
          positive_class: Case
    settings:
      - encoding: deeprc_enc
        ml_method: deeprc_model
    assessment:
      split_strategy: stratified_k_fold
      split_count: 5
      reports:
        models:
          - conf_matrix
          - perf_by_depth
    selection:
      split_strategy: stratified_k_fold
      split_count: 5
    strategy: GridSearch
    optimization_metric: balanced_accuracy
    metrics: [auc, balanced_accuracy, precision, recall, f1_micro]
    reports:
      - roc_summary
      - pr_curve
      - ml_settings
      - perf_by_sex
      - perf_by_age
    refit_optimal_model: true
    number_of_processes: 4
"""


def make_deeprc_configs():
    """6 training configs: 4 LOSO + 2 S2S sanity checks (PPMI, PDBP)."""
    out_dir = CONFIGS / "deeprc_v2"
    out_dir.mkdir(exist_ok=True)
    count = 0

    # 4 LOSO DeepRC — use filtered metadata (depth >= 1000)
    for study in STUDIES:
        train_meta = META_FILT / f"metadata_pid_{study}_loso.csv"
        name = f"deeprc_v2_LOSO_{study}"
        yaml_str = DEEPRC_TEMPLATE.format(
            repo=REPO, train_meta=train_meta, name=name)
        with open(out_dir / f"{name}.yaml", 'w') as f:
            f.write(yaml_str)
        count += 1

    # 2 S2S sanity checks from PPMI and PDBP (the larger cohorts)
    for source in ['ppmi', 'pdbp']:
        train_meta = META_FILT / f"metadata_pid_{source}.csv"
        name = f"deeprc_v2_S2S_{source}"
        yaml_str = DEEPRC_TEMPLATE.format(
            repo=REPO, train_meta=train_meta, name=name)
        with open(out_dir / f"{name}.yaml", 'w') as f:
            f.write(yaml_str)
        count += 1

    return count


# =========================================================================
# DeepRC MLApplication — apply each trained DeepRC model to target studies
# =========================================================================

DEEPRC_APPLY_TEMPLATE = """definitions:
  datasets:
    apply_dataset:
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

instructions:
  {name}:
    type: MLApplication
    dataset: apply_dataset
    config_path: {zip_path}
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


def make_deeprc_apply_configs():
    """For each trained DeepRC config × seed, apply to its target(s)."""
    out_dir = CONFIGS / "deeprc_apply_v2"
    out_dir.mkdir(exist_ok=True)
    count = 0

    # LOSO: train on not_X, test on X (filtered metadata)
    for study in STUDIES:
        train_name = f"deeprc_v2_LOSO_{study}"
        test_meta = META_FILT / f"metadata_pid_{study}.csv"
        for seed in range(5):
            zip_path = RESULTS / "deeprc_v2" / f"{train_name}_seed{seed}" / train_name / "optimal_case_control_other_latest" / "zip" / "ml_settings_case_control_other_latest.zip"
            inst_name = f"apply_LOSO_{study}_seed{seed}"
            yaml_str = DEEPRC_APPLY_TEMPLATE.format(
                repo=REPO, test_meta=test_meta, name=inst_name, zip_path=zip_path)
            with open(out_dir / f"{inst_name}.yaml", 'w') as f:
                f.write(yaml_str)
            count += 1

    # S2S: train on PPMI/PDBP, test on the other 3 studies
    for source in ['ppmi', 'pdbp']:
        train_name = f"deeprc_v2_S2S_{source}"
        for target in STUDIES:
            if target == source:
                continue
            test_meta = META_FILT / f"metadata_pid_{target}.csv"
            for seed in range(5):
                zip_path = RESULTS / "deeprc_v2" / f"{train_name}_seed{seed}" / train_name / "optimal_case_control_other_latest" / "zip" / "ml_settings_case_control_other_latest.zip"
                inst_name = f"apply_S2S_{source}_to_{target}_seed{seed}"
                yaml_str = DEEPRC_APPLY_TEMPLATE.format(
                    repo=REPO, test_meta=test_meta, name=inst_name, zip_path=zip_path)
                with open(out_dir / f"{inst_name}.yaml", 'w') as f:
                    f.write(yaml_str)
                count += 1

    return count


if __name__ == "__main__":
    n_apply = make_apply_configs()
    n_deeprc = make_deeprc_configs()
    n_deeprc_apply = make_deeprc_apply_configs()
    print(f"Generated {n_apply} MLApplication configs (CPU models → transfer)")
    print(f"Generated {n_deeprc} DeepRC v2 training configs")
    print(f"Generated {n_deeprc_apply} DeepRC v2 application configs")
