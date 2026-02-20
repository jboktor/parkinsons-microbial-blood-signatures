#!/usr/bin/env python
# coding: utf-8

import os
import sonia
from sonnia.sonnia import SoNNia
from sonnia.sonia import Sonia
from sonia.plotting import Plotter
from sonia.evaluate_model import EvaluateModel
from sonia.sequence_generation import SequenceGeneration
from sonnia.processing import Processing
import numpy as np
import pandas as pd
import argparse
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)

# Add this block at the top, after imports
logger.info("Parsing arguments")
parser = argparse.ArgumentParser(description="Run SoNNia pipeline with specified input data file.")
parser.add_argument('--input_data', type=str, required=True, 
                    help='Full path to the input csv file, should contain exact columns: amino_acid, v_gene, j_gene')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Directory where results and figures will be saved')
parser.add_argument('--output_csv', type=str, required=True,
                    help='Name of the output csv file (will be saved in output_dir)')
args = parser.parse_args()

# Ensure output directory and figures subdirectory exist
os.makedirs(args.output_dir, exist_ok=True)
figures_dir = os.path.join(args.output_dir, 'figures')
os.makedirs(figures_dir, exist_ok=True)

logger.info(f"Reading input data from {args.input_data}")
data_seqs = pd.read_csv(args.input_data)

# preprocess data
logger.info("Preprocessing data")
processor=Processing(pgen_model='humanTRB')
filtered=processor.filter_dataframe(data_seqs)

logger.info("Converting filtered data to list of sequences")
data_seqs = list(filtered.values.astype(str))
logger.info(f"First 3 sequences: {data_seqs[:3]}")

# initialize model
logger.info("Initializing SoNNia model")
qm = SoNNia(data_seqs=data_seqs,pgen_model='humanTRB')
# qm_linear=Sonia(data_seqs=data_seqs,pgen_model='humanTRB')

# add generated sequences (you can add them from file too, more is better.)
logger.info("Adding generated sequences")
qm.add_generated_seqs(int(5e5)) 
# qm_linear.add_generated_seqs(int(5e5))

#train model
logger.info("Training model")
qm.infer_selection(epochs=30,batch_size=int(1e4))
# qm_linear.infer_selection(epochs=30,batch_size=int(1e4))

# plotting deep learning model
import matplotlib.pyplot as plt
plot_sonia=Plotter(qm)
logger.info("Plotting model learning curve")
fig = plot_sonia.plot_model_learning()
plt.savefig(os.path.join(figures_dir, 'model_learning.png'))
plt.close()
# logger.info("Plotting logQ")
# fig = plot_sonia.plot_logQ()
# plt.savefig(os.path.join(figures_dir, 'logQ.png'))
# plt.close()
# plot_sonia.plot_vjl()  # Disabled due to previous errors

# # Generate sequences
logger.info("Generating sequences from pgen")
pre_seqs=qm.generate_sequences_pre(int(1e4))

logger.info("Generating sequences from ppost")
post_seqs=qm.generate_sequences_post(int(1e4))

# # Evaluate sequences
logger.info("Evaluating sequences")
Q_data,pgen_data,ppost_data=qm.evaluate_seqs(qm.data_seqs[:int(1e4)])
Q_gen,pgen_gen,ppost_gen=qm.evaluate_seqs(pre_seqs)
Q_model,pgen_model,ppost_model=qm.evaluate_seqs(post_seqs)
logger.info(f"Q_model[:3]: {Q_model[:3]}")
logger.info(f"pgen_model[:3]: {pgen_model[:3]}")
logger.info(f"ppost_model[:3]: {ppost_model[:3]}")

# Save probability plots as PNGs
try:
    logger.info("Plotting and saving P_{pre}")
    fig = plot_sonia.plot_prob(data=pgen_data,gen=pgen_gen,model=pgen_model,ptype='P_{pre}')
    plt.savefig(os.path.join(figures_dir, 'P_pre.png'))
    plt.close()
    logger.info("Plotting and saving P_{post}")
    fig = plot_sonia.plot_prob(ppost_data,ppost_gen,ppost_model,ptype='P_{post}')
    plt.savefig(os.path.join(figures_dir, 'P_post.png'))
    plt.close()
    logger.info("Plotting and saving Q")
    fig = plot_sonia.plot_prob(Q_data,Q_gen,Q_model,ptype='Q',bin_min=-4,bin_max=2)
    plt.savefig(os.path.join(figures_dir, 'Q.png'))
    plt.close()
except Exception as e:
    logger.error(f"Plotting failed: {e}")

# Evaluating full dataset and saving results
logger.info("Evaluating full dataset")
Q_data,pgen_data,ppost_data=qm.evaluate_seqs(qm.data_seqs)

# save results, including input data columns
output_csv_path = os.path.join(args.output_dir, args.output_csv)
logger.info(f"Saving results to {output_csv_path}")
results_df = pd.DataFrame({
    'amino_acid': filtered['amino_acid'].values,
    'v_gene': filtered['v_gene'].values,
    'j_gene': filtered['j_gene'].values,
    'Q_data': Q_data,
    'pgen_data': pgen_data,
    'ppost_data': ppost_data
})

# save results
results_df.to_csv(output_csv_path, index=False)

# save model
logger.info("Saving model")
qm.save_model(os.path.join(args.output_dir, 'model'))

logger.info("Pipeline completed successfully")
