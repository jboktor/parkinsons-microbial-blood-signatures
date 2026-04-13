#!/usr/bin/env python3
"""Generate all S2S and LOSO ImmuneML YAML configs for the transfer-study matrix."""

import os
from pathlib import Path
from itertools import permutations

CONFIGS_DIR = Path("/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs")
META_DIR = Path("/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/metadata")
REPO_DIR = Path("/central/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/data_input/participant_repertoires")

STUDIES = ["pdbp", "ppmi", "biofind", "hbs"]

DATASET_BLOCK = """    {name}:
      format: AIRR
      params:
        path: {repo_dir}
        metadata_file: {meta_file}
        is_repertoire: true
        import_productive: true
        import_with_stop_codon: false
        import_out_of_frame: false
        import_illegal_characters: false
        import_empty_nt_sequences: true
        import_empty_aa_sequences: false
        region_type: IMGT_CDR3"""

LABEL_BLOCK = """    labels:
      - case_control_other_latest:
          positive_class: Case"""

# ── Model templates ──────────────────────────────────────────────────────────

MODELS = {
    "kmer_logreg": {
        "encodings": """
    kmer3_all:
      KmerFrequency:
        k: 3
        sequence_encoding: continuous_kmer
        normalization_type: relative_frequency
        reads: all
        scale_to_unit_variance: true
        scale_to_zero_mean: true
        sequence_type: amino_acid
        region_type: IMGT_CDR3""",
        "ml_methods": """
    log_reg_l1:
      LogisticRegression:
        penalty: l1
        C: [0.01, 0.1, 1.0]
        max_iter: 5000
      model_selection_cv: true
      model_selection_n_folds: 3
    log_reg_l2:
      LogisticRegression:
        penalty: l2
        C: [0.01, 0.1, 1.0]
        max_iter: 5000
      model_selection_cv: true
      model_selection_n_folds: 3""",
        "settings": """
      - encoding: kmer3_all
        ml_method: log_reg_l1
      - encoding: kmer3_all
        ml_method: log_reg_l2""",
        "reports_def": """
    coefficients:
      Coefficients:
        coefs_to_plot: [n_largest]
        n_largest: [50]
    settings_perf: MLSettingsPerformance
    roc_summary: ROCCurveSummary""",
        "assessment_reports": """
        models:
          - coefficients""",
        "instruction_reports": """
      - settings_perf
      - roc_summary""",
        "metrics": "[auc, balanced_accuracy, precision, recall]",
    },

    "kmer_svm": {
        "encodings": """
    kmer3_all:
      KmerFrequency:
        k: 3
        sequence_encoding: continuous_kmer
        normalization_type: relative_frequency
        reads: all
        scale_to_unit_variance: true
        scale_to_zero_mean: true
        sequence_type: amino_acid
        region_type: IMGT_CDR3""",
        "ml_methods": """
    svm_linear:
      SVM:
        kernel: linear
        C: [0.01, 0.1, 1.0]
      model_selection_cv: true
      model_selection_n_folds: 3""",
        "settings": """
      - encoding: kmer3_all
        ml_method: svm_linear""",
        "reports_def": """
    roc_summary: ROCCurveSummary""",
        "assessment_reports": "",
        "instruction_reports": """
      - roc_summary""",
        "metrics": "[auc, balanced_accuracy, precision, recall]",
    },

    "kmer_rf": {
        "encodings": """
    kmer3_all:
      KmerFrequency:
        k: 3
        sequence_encoding: continuous_kmer
        normalization_type: relative_frequency
        reads: all
        scale_to_unit_variance: true
        scale_to_zero_mean: true
        sequence_type: amino_acid
        region_type: IMGT_CDR3""",
        "ml_methods": """
    random_forest:
      RandomForestClassifier:
        n_estimators: 100
        class_weight: balanced""",
        "settings": """
      - encoding: kmer3_all
        ml_method: random_forest""",
        "reports_def": """
    coefficients:
      Coefficients:
        coefs_to_plot: [n_largest]
        n_largest: [50]
    roc_summary: ROCCurveSummary""",
        "assessment_reports": """
        models:
          - coefficients""",
        "instruction_reports": """
      - roc_summary""",
        "metrics": "[auc, balanced_accuracy, precision, recall]",
    },

    "pubclone": {
        "encodings": """
    seq_abundance:
      SequenceAbundance:
        p_value_threshold: 0.1""",
        "ml_methods": """
    prob_binary:
      ProbabilisticBinaryClassifier:
        max_iterations: 200
        update_rate: 0.01""",
        "settings": """
      - encoding: seq_abundance
        ml_method: prob_binary""",
        "reports_def": """
    roc_summary: ROCCurveSummary""",
        "assessment_reports": "",
        "instruction_reports": """
      - roc_summary""",
        "metrics": "[auc, balanced_accuracy]",
    },

    "compairr": {
        "encodings": """
    compairr_1mm:
      CompAIRRSequenceAbundance:
        compairr_path: /central/groups/MazmanianLab/jboktor/software/compairr/src/compairr
        p_value_threshold: 0.001
        ignore_genes: false""",
        "ml_methods": """
    log_reg_compairr:
      LogisticRegression:
        penalty: l1
        C: [0.1, 1.0]
      model_selection_cv: true
      model_selection_n_folds: 3""",
        "settings": """
      - encoding: compairr_1mm
        ml_method: log_reg_compairr""",
        "reports_def": """
    roc_summary: ROCCurveSummary""",
        "assessment_reports": "",
        "instruction_reports": """
      - roc_summary""",
        "metrics": "[auc, balanced_accuracy]",
    },
}


def make_yaml(train_name, train_meta, test_name, test_meta, model_key, model_cfg, design):
    """Generate a complete YAML config string."""
    train_ds = DATASET_BLOCK.format(name="train_dataset", repo_dir=REPO_DIR, meta_file=train_meta)
    instruction_name = f"{design}_{train_name}_to_{test_name}_{model_key}"

    assessment_reports_block = ""
    if model_cfg["assessment_reports"]:
        assessment_reports_block = f"""
      reports:{model_cfg['assessment_reports']}"""

    yaml = f"""definitions:
  datasets:
{train_ds}

  encodings:{model_cfg['encodings']}

  ml_methods:{model_cfg['ml_methods']}

  reports:{model_cfg['reports_def']}

instructions:
  {instruction_name}:
    type: TrainMLModel
    dataset: train_dataset
{LABEL_BLOCK}
    settings:{model_cfg['settings']}
    assessment:
      split_strategy: random
      split_count: 1
      training_percentage: 0.7{assessment_reports_block}
    selection:
      split_strategy: k_fold
      split_count: 3
    strategy: GridSearch
    optimization_metric: balanced_accuracy
    metrics: {model_cfg['metrics']}
    reports:{model_cfg['instruction_reports']}
    refit_optimal_model: true
    number_of_processes: 4
"""
    return yaml


def main():
    os.makedirs(CONFIGS_DIR / "S2S", exist_ok=True)
    os.makedirs(CONFIGS_DIR / "LOSO", exist_ok=True)
    count = 0

    # S2S: all pairwise combinations (12 pairs × 5 models = 60 configs)
    for train_study, test_study in permutations(STUDIES, 2):
        train_meta = META_DIR / f"metadata_pid_{train_study}.csv"
        test_meta = META_DIR / f"metadata_pid_{test_study}.csv"
        for model_key, model_cfg in MODELS.items():
            fname = f"S2S_{train_study}_to_{test_study}_{model_key}.yaml"
            yaml_content = make_yaml(train_study, train_meta, test_study, test_meta,
                                      model_key, model_cfg, "S2S")
            outpath = CONFIGS_DIR / "S2S" / fname
            with open(outpath, 'w') as f:
                f.write(yaml_content)
            count += 1

    # LOSO: 4 pairs × 5 models = 20 configs
    for study in STUDIES:
        train_meta = META_DIR / f"metadata_pid_{study}_loso.csv"
        test_meta = META_DIR / f"metadata_pid_{study}.csv"
        for model_key, model_cfg in MODELS.items():
            fname = f"LOSO_{study}_{model_key}.yaml"
            yaml_content = make_yaml(f"not_{study}", train_meta, study, test_meta,
                                      model_key, model_cfg, "LOSO")
            outpath = CONFIGS_DIR / "LOSO" / fname
            with open(outpath, 'w') as f:
                f.write(yaml_content)
            count += 1

    print(f"Generated {count} YAML configs")
    print(f"  S2S:  {len(list(permutations(STUDIES, 2))) * len(MODELS)} in {CONFIGS_DIR / 'S2S'}")
    print(f"  LOSO: {len(STUDIES) * len(MODELS)} in {CONFIGS_DIR / 'LOSO'}")


if __name__ == "__main__":
    main()
