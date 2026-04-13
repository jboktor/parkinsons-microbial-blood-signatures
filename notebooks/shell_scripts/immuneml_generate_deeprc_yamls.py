#!/usr/bin/env python3
"""Generate DeepRC YAML configs for all S2S and LOSO pairs."""

import os
from pathlib import Path
from itertools import permutations

CONFIGS_DIR = Path("/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/configs")
META_DIR = Path("/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/metadata")
REPO_DIR = Path("/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/airr/ImmuneML/data_input/participant_repertoires")

STUDIES = ["pdbp", "ppmi", "biofind", "hbs"]

TEMPLATE = """definitions:
  datasets:
    train_dataset:
      format: AIRR
      params:
        path: {repo_dir}
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
        kernel_size: 5
        n_kernels: 16
        n_additional_convs: 1
        n_attention_network_layers: 2
        n_attention_network_units: 32
        n_output_network_units: 32
        learning_rate: 0.00005
        n_updates: 10000
        evaluate_at: 2000
        sample_n_sequences: 2000
        l2_weight_decay: 0.001
        pytorch_device_name: cuda:0

  reports:
    roc_summary: ROCCurveSummary

instructions:
  {instruction_name}:
    type: TrainMLModel
    dataset: train_dataset
    labels:
      - case_control_other_latest:
          positive_class: Case
    settings:
      - encoding: deeprc_enc
        ml_method: deeprc_model
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
    number_of_processes: 4
"""


def main():
    os.makedirs(CONFIGS_DIR / "S2S", exist_ok=True)
    os.makedirs(CONFIGS_DIR / "LOSO", exist_ok=True)
    count = 0

    # S2S: 12 pairs
    for train_study, test_study in permutations(STUDIES, 2):
        train_meta = META_DIR / f"metadata_pid_{train_study}.csv"
        instruction_name = f"S2S_{train_study}_to_{test_study}_deeprc"
        fname = f"S2S_{train_study}_to_{test_study}_deeprc.yaml"
        yaml_content = TEMPLATE.format(
            repo_dir=REPO_DIR, train_meta=train_meta,
            instruction_name=instruction_name)
        with open(CONFIGS_DIR / "S2S" / fname, 'w') as f:
            f.write(yaml_content)
        count += 1

    # LOSO: 4 pairs
    for study in STUDIES:
        train_meta = META_DIR / f"metadata_pid_{study}_loso.csv"
        instruction_name = f"LOSO_not_{study}_to_{study}_deeprc"
        fname = f"LOSO_{study}_deeprc.yaml"
        yaml_content = TEMPLATE.format(
            repo_dir=REPO_DIR, train_meta=train_meta,
            instruction_name=instruction_name)
        with open(CONFIGS_DIR / "LOSO" / fname, 'w') as f:
            f.write(yaml_content)
        count += 1

    print(f"Generated {count} DeepRC YAML configs (12 S2S + 4 LOSO)")


if __name__ == "__main__":
    main()
