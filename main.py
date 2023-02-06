# This file demonstrates the usage of the analyses in function.py file

# NOTE: some of the analyses use a specific set of genes - "DE genes".
# these genes are differentially expressed (q-value<0.01 in the DE analysis in response to dsRNA)
# in at least 1 species among the 6 species: Human, Macaque, Mouse, Rat, Rousettus , Pipistrellus

import numpy as np
import pandas as pd
#import import_ipynb
from functions import *
import matplotlib.pyplot as plt

# Get fully vectorized PDFs, incl. font:
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Define set of "DE genes": q-value in the DE analysis is less than 0.01 in at least 1 of the 6 species
metadata_one_to_one = pd.read_excel(r'Supporting_Tables.xlsx', sheet_name='Table S9', index_col=0, header=2)
species_list = [i.strip(' FC') for i in metadata_one_to_one.columns if ' FC' in i]
DE_genes_df = metadata_one_to_one.loc[[g for g in metadata_one_to_one.index if min(
    [metadata_one_to_one.loc[g, '{} Q-value'.format(s)] for s in species_list]) <= 0.01]]
# Remove genes with no DE analysis results
DE_genes_df = DE_genes_df.dropna()

# ---------------------------------------------------------------------------
# Analysis#1: calculating spearman correlation
# ---------------------------------------------------------------------------
# Calculating spearman correlation of FC values between each 2 species in "DE genes"
conbinations = [('Homo sapiens FC', 'Macaca mulatta FC'),
                ('Rattus norvegicus FC', 'Mus musculus FC'),
                ('Homo sapiens FC', 'Mus musculus FC'),
                ('Rousettus aegyptiacus FC', 'Pipistrellus kuhlii FC')
                ]
for i in conbinations:
    scatter_correlation(species1_fc=i[0], species2_fc=i[1], df=DE_genes_df)


# ---------------------------------------------------------------------------
# Analysis#2: features distribution in 3 divergent groups (high, medium, low)
# ---------------------------------------------------------------------------
# In this analysis some numeric biological and evolutionary features's distribution are analysed in 3 divergent groups of genes.
# The groups are defined from the 'DE_genes', and split according to:
# measure_of_divergence- 'distance score' (column name='distance_rous-pip_non-normalized') between Rousettus aegyptiacus and Pipistrellus kuhlii with absolute ('absolute'=True) value.
# where 'distance score' is the measure of divergence in gene expression between these 2 species and it defined as: (Rousettus aegyptiacus FC)-(Pipistrellus kuhlii FC).
# The FC values are from the DE analysis of each species between the control and stimulated with dsRNA condition
# The side of Mann-Whitney test for each feature is defined from biological and evolutionary perspective

# process the full table of all features data - subset 'DE_genes'
metadata_rous = pd.read_excel(r'Supporting_Tables.xlsx', sheet_name='Table S7', index_col=0, header=2)
human_to_rous_gene_names = dict(zip(metadata_rous['ortho_homo_gene'], metadata_rous.index))
DE_genes_df.rename(human_to_rous_gene_names, inplace=True)
metadata_rous_de_genes = metadata_rous.loc[DE_genes_df.index]

# prepare features to analyse and Mann-Whitney test side:
features_mann_dict = {'P_val_gene_gain_and_loss_ENSEMBL': 'less',
                      # rate of gene gain and loss statistics (-log(p-value gene gain and loss)) is lower in 'low divergence group' than the 'high divergence group'
                      'pnas_dNdS_P': 'less',
                      # rate of dn\ds ratio statistics (-log(p-value dn\ds)) is lower in 'low divergence group' than the 'high divergence group'
                      'seq_similarity_rous-pip_per_identity': 'greater',
                      # percentage of identity in the coding sequence between rous and pip is higher in 'low divergence group' than the 'high divergence group'
                      'DM_in_single-cells_LF_RA': 'two-sided',  # cell-to-cell variability (DM)
                      'DM_in_single-cells_PIC_RA': 'two-sided',
                      'DM_in_single-cells_LF_HS': 'two-sided',
                      'DM_in_single-cells_PIC_HS': 'two-sided'
                      }

for feature, mann_whitney_side in features_mann_dict.items():
    feature_in_3_divergent_groups(df=metadata_rous_de_genes,
                                  feature=feature,
                                  mann_whitney_side=mann_whitney_side)
    # default values: res_dir_name='feature_in_3_divergent_groups_results', measure_of_divergence='distance_rous-pip_non-normalized', absolute=True


# ---------------------------------------------------------------------------
# Analysis#3: Resistance and Tolerance measures in Rousettous and Pipistrellus specific genes
# ---------------------------------------------------------------------------
# This analysis explores whether Rousettus and Pipistrellus specific genes (up-regulated in response to dsRNA in one of them compared to the other)
# are more associated with expression programs that underlie ‘disease tolerance’ and ‘disease resistance’.
# The specific genes are defined according to 'distance score' between the 2 bats, as top/bottom 20% scores (frac=0.2).
# The distribution of 4 measures (2 for resistance and 2 for tolerance) is compared between both groups of bats specific genes (Mann-Whitney test).
# This analysis is done once on 'DE_genes', and once on all 1-1 orthologs genes between these bats.

# pre-process data: prepare dicts and index renaming (to human ENSEMBL name)
# humanENSEMBL to general-gene-name dict
metadata = pd.read_csv(r'1-1_human_mouse_metadata.csv', index_col=0)
metadata_1_1 = metadata[metadata['Mouse homology type'] == 'ortholog_one2one']
metadata_1_1["Gene name human"] = metadata_1_1["Gene name human"].str.lower()
gene_name_to_eng_dict_1_1 = dict(zip(metadata_1_1['Gene name human'], metadata_1_1.index))
# renaming Rousettus to humanENSEMBL
metadata_rous = pd.read_excel(r'Supporting_Tables.xlsx', sheet_name='Table S7', index_col=0, header=2)
rous_to_human_gene_names = dict(zip(metadata_rous.index, metadata_rous['ortho_homo_gene']))
metadata_rous.set_index('ortho_homo_gene', inplace=True)
# renaming T&R table index to humanENSEMBL name
tol_res_df = pd.read_excel(r'Supporting_Tables.xlsx', sheet_name='Table S16', index_col=0, header=2)
tol_res_df.index = tol_res_df.index.str.lower()
tol_res_df.rename(index=gene_name_to_eng_dict_1_1, inplace=True)

# merge data
a = metadata_rous[
    ['distance_rous-pip_non-normalized']].dropna()
merged = a.merge(tol_res_df, how='left', left_on=a.index, right_on=tol_res_df.index)
merged = merged.dropna(
    subset=['key_0'])
merged.set_index('key_0', inplace=True)

# process data
# example 1 - all 1-1 genes between 2 bats
tolerance_and_resistance_in_2_bats(df_DE_gene=merged,
                                   name_DE_genes='all_1-1_genes',
                                   frac=0.2
                                   )
# example 2 - 'DE_genes'
DE_genes = [rous_to_human_gene_names[i] for i in DE_genes_df.index]
tolerance_and_resistance_in_2_bats(df_DE_gene=merged.loc[DE_genes],
                                   name_DE_genes='DE_genes_in_at_least_1_of_6species',
                                   frac=0.2
                                   )


# ---------------------------------------------------------------------------
# Analysis#4: Gene expression in comparative tissues of species
# ---------------------------------------------------------------------------
# this analysis explores the up-regulation level (FC values) in response to dsRNA
# and basal expression level (normalized log(TPM)) in comparative tissues between human, mouse and Rousettus.
# it is done on specific groups of genes from different viral RNA-sensing pathway: 'RLRs', 'TLRs'.
# create heatmaps and conduct paired t-tests for the basal expression level in tissues between bat and human and bat and mouse.

# pre-process data:
# tissues df
metadata_one_to_one = pd.read_excel(r'Supporting_Tables.xlsx', sheet_name='Table S9', index_col=0, header=2)
metadata_one_to_one.set_index('Rousettus aegyptiacus_ortho_gene', inplace=True)
tissues_df = metadata_one_to_one.loc[:, metadata_one_to_one.columns.str.contains('_bat|_mouse|_human')]
# FC df
FC_df = metadata_one_to_one.loc[:,
        metadata_one_to_one.columns.str.contains('FC') & metadata_one_to_one.columns.str.contains('Mus|Homo|Rous')]
FC_df.rename(columns={'Homo sapiens FC': 'human',
                      'Mus musculus FC': 'mouse',
                      'Rousettus aegyptiacus FC': 'bat'
                      }, inplace=True)

# analyse RLRs pathway genes
speciel_genes = ['IRF3', 'DDX58', 'IFIH1', 'DHX58', 'MAVS', 'NFKB2', 'PPP1CA', 'PPP1CC', 'PRKCB', 'PRKCA', 'TRIM25',
                 'USP21', 'TBK1']
tissues_basal_expression_heatmaps(tissues_df=tissues_df, FC_df=FC_df, genes_list=speciel_genes, genes_list_name='RLRs')
# noramlization=True (default)
# 'ttests'=[('bat','human'),('bat','mouse')] (default)

# analyse TLRs pathway genes
speciel_genes = ['TLR3', 'TRAF3', 'TBK1', 'IKBKE', 'IRF3', 'RIPK1', 'CHUK', 'IKBKB', 'NFKB2', 'TICAM1']
tissues_basal_expression_heatmaps(tissues_df=tissues_df, FC_df=FC_df, genes_list=speciel_genes, genes_list_name='TLRs')
# noramlization=True (default)
# 'ttests'=[('bat','human'),('bat','mouse')] (default)
