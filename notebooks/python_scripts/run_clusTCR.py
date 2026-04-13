#!/usr/bin/env python
# coding: utf-8

# This script should be run with the clustcr mamba environment

from clustcr import Clustering
import pandas as pd
import pickle

wkdir = '/resnick/groups/MazmanianLab/jboktor/PDMBS/pdairr/data/interim/clusTCR'
# clustering = Clustering()
clustering = Clustering(n_cpus=31)

# from clustcr import datasets
# cdr3 = datasets.test_cdr3()

# reading in TCRB data
cdr3_data = pd.read_table(f'{wkdir}/input/amppd_tcrb_2025-08-04.tsv')
# cdr3_data_test = cdr3_data.head(10000)
# print(cdr3_data_test.head())

output_with_vgene = clustering.fit(cdr3_data, include_vgene=True, cdr3_col='CDR3b_aa', v_gene_col='TRBV')
print(output_with_vgene.clusters_df.head())

# Save pickle file
with open(f'{wkdir}/results/clusTCR_output_vgene_restricted_31cpus.pkl', 'wb') as f:
    pickle.dump(output_with_vgene, f)

# Save csv files
output_with_vgene.summary().to_csv(f'{wkdir}/results/clusTCR_summary_vgene_restricted_31cpus.csv', index=True)
output_with_vgene.clusters_df.to_csv(f'{wkdir}/results/clusTCR_output_vgene_restricted_31cpus.csv', index=True)

# output = clustering.fit(cdr3_data_test, include_vgene=False, cdr3_col='CDR3b_aa')
# print(output.clusters_df.head())
# output.summary().to_csv(f'{wkdir}/results/clusTCR_summary.csv', index=True)
# output.clusters_df.to_csv(f'{wkdir}/results/clusTCR_output.csv', index=True)
