import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
import statsmodels.stats.multitest as sm
import scipy
from collections import defaultdict
import os


# get fully vectorized PDFs, incl. font:
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# ---------------------------------------------------------------------------
# Analysis#1: calculating spearman correlation
# ---------------------------------------------------------------------------

def scatter_correlation(species1_fc, species2_fc, df, is_FC=True, res_dir_name='scatter_correlation_results'):
    """
    This function calculates spearman correlation between two numeric variables and plots a scatter plot.

    :param species1_fc: string. column name of the variable #1 in 'df' param.
    :param species2_fc: string. column name of the variable #2 in 'df' param.
    :param df: pandas dataframe. containing species1_fc and species2_fc columns (numeric values).
    NOTE: if the analyzed variables are FC values of two species the species1_fc,species2_fc params
    should be in format: '{species name} FC'
    :param is_FC: boolean.
    =True: (default) the analysed variables are FC values of 2 species &
    species1_fc and species2_fc params are in format '{species name} FC'
    =False: else.
    :param res_dir_name: string. name of the results directory.
    :return: pdf file with spearman correlation results and plot saved to results directory.
    """

    # define species names
    if is_FC:
        species1_name = species1_fc.strip(' FC')
        species2_name = species2_fc.strip(' FC')
    else:
        species1_name = species1_fc
        species2_name = species2_fc

    # calculating spearman correlation
    spearman_corr_DE, pval_corr_DE = scipy.stats.spearmanr(a=df[species1_fc], b=df[species2_fc])
    corr_text = f'correlation: {round(spearman_corr_DE, 3)}\np-value:{pval_corr_DE}'

    # Plot data and a linear regression model fit.
    sns.set_style("ticks", {"axes.facecolor": ('white'), 'patch.edgecolor': 'black'})
    plt.figure(figsize=[8, 5])
    plt.title('FC correlation {} - {} \n{}'.format(species1_name, species2_name, corr_text), size=14)
    g = sns.regplot(data=df, x=species1_fc, y=species2_fc, scatter=True, scatter_kws={'color': '#838383', 's': 25},
                    color='#191919')
    plt.xlabel(species1_fc, size=12)
    plt.ylabel(species2_fc, size=12)
    # add y=x line
    x0, x1 = g.get_xlim()
    y0, y1 = g.get_ylim()
    lims = [max(x0, y0), min(x1, y1)]
    g.plot(lims, lims, linestyle='--', color='black', lw=0.5, scalex=False, scaley=False)
    sns.despine(offset=10, trim=True)

    # save results
    results_dir = r'results\{}'.format(res_dir_name)
    os.makedirs(results_dir, exist_ok=True)
    plt.savefig(r'{}\{}_vs_{}_FC_scatter_plot.pdf'.format(results_dir, species1_name, species2_name))

    print(f'scatter_correlation for {species1_fc}-{species2_fc} - DONE!')
# ---------------------------------------------------------------------------
# Analysis#2: features distribution in 3 divergent groups (high, medium, low)
# ---------------------------------------------------------------------------

def splitting_to_3_divergent_groups(df, measure_of_divergence, absolute):

    """
    # This function splits the data into 3 groups – high, medium and low according to specific numeric measure - measure_of_divergence
    # it can work with any measure that differentiate high, medium and low values

    :param df: pandas dataframe containing the relevant rows (originally genes) and column of the measure_of_divergence (numeric values)
    :param measure_of_divergence: string. measure of divergence column name.
    :param absolute: boolean. =True - get and use absolute values of measure_of_divergence for splitting the data; =Falue - else.
    # NOTE: the 'measure_of_divergence' and the 'absolute' arguments should be correlated and meet the biological logic behind it.
    :return: df-data as pandas dataframe containing 'divergence_group' column indicating the 3 divergence levels
    """

    # initialize measure_of_divergence
    if absolute:
        df['measure_of_divergence'] = df[measure_of_divergence].abs()
    else:
        df['measure_of_divergence'] = df[measure_of_divergence]

    # data processing
    distance_file_sorted = df.sort_values(by='measure_of_divergence')
    splitted = np.array_split(distance_file_sorted, 3)
    # define 3 divergence groups
    low_v_genes = splitted[0]
    medium_v_genes = splitted[1]
    high_v_genes = splitted[2]
    low_v_genes['divergence_group'] = 'low_v_genes'
    medium_v_genes['divergence_group'] = 'medium_v_genes'
    high_v_genes['divergence_group'] = 'high_v_genes'

    total_groups = pd.concat([low_v_genes, medium_v_genes, high_v_genes], axis=0)

    return total_groups


def feature_in_3_divergent_groups(feature, df, mann_whitney_side,
                                  measure_of_divergence='distance_rous-pip_non-normalized', absolute=True,
                                  res_dir_name='feature_in_3_divergent_groups_results'):
    """
    # This function generates boxplots of a feature distribution for 3 divergence groups according to measure_of_divergence
    # with significance annotations of Mann-Whitney test between the low and high divergent groups.

    :param feature: string. column name of the feature to be analyzed.
    :param df: pandas dataframe containing the relevant rows (genes mostly), the columns of the 'feature' and 'measure_of_divergence'
    # NOTE: both 'feature' & 'measure_of_divergence' values should be numeric
    #      'feature' data should follow Mann-Whitney test assumptions.
    #      all features provided must be processed in the df (except for 'pnas_dNdS_P' &'P_val_gene_gain_and_loss_ENSEMBL' which are processed in the function.
    :param mann_whitney_side: string indicating the side state of the test:
    # ='less': the low divergent group is lower than the 'high'
    # ='greater': the low divergent group is greater than the 'high'
    # ='two-sided': the test is two-sided
    # NOTE: defined from biological and evolutionary perspective
    :param measure_of_divergence: string. measure of divergence column name.
    :param absolute: boolean. =True - get and use absolute values of measure_of_divergence for splitting the data; =Falue - without absolute.
    NOTE: by default - measure_of_divergence is the 'distance score' (name='distance_rous-pip_non-normalized') between Rousettus aegyptiacus and Pipistrellus kuhlii with absolute ('absolute'=True) value.
    # where 'distance score' is the measure of divergence in gene expression between these 2 species and it defined as: (Rousettus aegyptiacus FC)-(Pipistrellus kuhlii FC).
    # the FC values are from the DE analysis of each species between the control and stimulated with dsRNA condition
    :param res_dir_name: string. name of the results directory.
    :return: plot's pdf file saved to results directory.
    """

    # split to 3 divergence groups
    total_groups = splitting_to_3_divergent_groups(df=df, measure_of_divergence=measure_of_divergence,
                                                   absolute=absolute)
    # pre-process:
    # remove genes with no 'feature' data
    total_groups = total_groups.dropna(subset=[feature])
    if feature == 'seq_similarity_rous-pip_per_identity':
        total_groups = total_groups[total_groups['seq_similarity_rous-pip_per_identity'] != '-']
        total_groups = total_groups.astype(
            {'seq_similarity_rous-pip_per_identity': float})  # convert column type to numeric

    # log transformation for specific features
    if (feature == 'pnas_dNdS_P') | (feature == 'P_val_gene_gain_and_loss_ENSEMBL'):
        feature_new = f'minus_log({feature})'
        total_groups[feature_new] = -np.log10(total_groups[feature] + 0.000001)
    else:
        feature_new = feature

    # mann-whitney test labels
    if mann_whitney_side == 'less':
        mann_test = 'Mann-Whitney-ls'
        mann_test_short = 'less'
    elif mann_whitney_side == 'greater':
        mann_test = 'Mann-Whitney-gt'
        mann_test_short = 'greater'
    else:
        mann_test = 'Mann-Whitney'
        mann_test_short = 'two-sided'


    # plot results:
    order_boxplot = ['high_v_genes', 'medium_v_genes', 'low_v_genes']
    mypalette = {'high_v_genes': '#5e43a5', 'medium_v_genes': '#6785be', 'low_v_genes': '#9ebbc9'}
    fig = plt.figure(figsize=[6, 9])
    box = sns.boxplot(data=total_groups, x='divergence_group', y=feature_new, palette=mypalette, showfliers=False,
                      order=order_boxplot
                      , width=0.7)
    sns.despine(offset=10, trim=False, fig=fig)
    box.set_xticklabels(['high', 'medium', 'low'], fontsize=14)
    box.set_yticklabels([round(num, 2) for num in box.get_yticks()], fontsize=14)
    box.set_xlabel('Divergence group', size=18)
    box.set_ylabel(feature_new, size=18)
    # significance annotations in figure - Mann-Whitney p-value for the relevant test side
    add_stat_annotation(box, data=total_groups, x='divergence_group', y=feature_new,
                        box_pairs=[("low_v_genes", "high_v_genes")],
                        test=mann_test, text_format='star', loc='outside', verbose=2)

    # print feature name and Mann-Whitney results
    mann_w = scipy.stats.mannwhitneyu(
        x=total_groups[total_groups['divergence_group'] == 'low_v_genes'][feature_new],
        y=total_groups[total_groups['divergence_group'] == 'high_v_genes'][feature_new],
        alternative=mann_test_short)
    print(feature_new)
    print(mann_w)

    # save results
    results_path = 'results\{}'.format(res_dir_name)
    os.makedirs(results_path, exist_ok=True)
    fig_path = r'{}\{}.pdf'.format(results_path, feature_new)
    plt.savefig(fig_path)

    print(f'feature_in_3_divergent_groups for {feature} - DONE!')



# ---------------------------------------------------------------------------
# Analysis#3: Resistance and Tolerance measures in Rousettous and Pipistrellus specific genes
# ---------------------------------------------------------------------------

def tolerance_and_resistance_in_2_bats(df_DE_gene, name_DE_genes, frac=0.2,
                                       res_dir_name='tolerance&resistance_in_bats'):
    """
    # This function defines Rousettus and Pipistrellus specific genes (up-regulated in response to dsRNA in one of them compared to the other):
        #  based on 'distance score' ((Rousettus_FC)-(Pipitrellus_FC)) which is the measure of divergence in gene expression between these 2 bats
        #  specific genes are defined as having top or bottom distance scores given by the 'frac' input
        #  top scores are more Rousettus up-regulated genes in response to dsRNA compared to Pipistrellus, and bottom scores are the other way around
    # For each one of the Tolerance/Resistance axis (different methods) Mann-Whitney test is conducted between Rousettus and Pipistrellus specific genes
    # as a result it generates boxplots for Tolerance/Resistance measures values in Rousettus and Pipistrellus specific genes.

    :param df_DE_gene: pandas dataframe. Containing genes of interest as rows and columns for 'distance score' (name: distance_rous-pip_non-normalized) and tolerance&resistance axis values.
    :param name_DE_genes: string. name of group of genes to analyse, for the results path
    :param frac: 0<float<1 - 0.2 as default, the fraction of genes among all the genes to be defined as top (more specific to Rousettus) or bottom (more specific to Pipistrellus)
    :param res_dir_name: string. name of the results directory.
    :return: plots' pdf files for each tolerance&resistance measures saved to results directory.
    """

    # define Rousettus and Pipistrellus specific genes based on 'distance score' and 'frac' :
    relevant_col = 'distance_rous-pip_non-normalized'
    distance_file_sorted_minus_top = df_DE_gene.sort_values(by=relevant_col)
    top_pip = distance_file_sorted_minus_top.head(round(len(df_DE_gene) * frac))
    distance_file_sorted_plus_top = df_DE_gene.sort_values(by=relevant_col, ascending=False)
    top_rous = distance_file_sorted_plus_top.head(round(len(df_DE_gene) * frac))
    df_DE_gene.loc[top_pip.index, 'species'] = 'Pipistrellus'
    df_DE_gene.loc[top_rous.index, 'species'] = 'Rousettus'

    data = df_DE_gene.dropna() # remove genes with no tolerance&resistance values

    my_pal = {"Pipistrellus": '#b96363',
              "Rousettus": '#a83d3d'
              }
    results_dir = r'results\{}\{}'.format(res_dir_name, name_DE_genes)
    os.makedirs(results_dir, exist_ok=True)
    # plot results for each one of the tolerance\resistance measures
    for T_or_R in [c for c in data.columns if 'axis' in c or '-to-' in c]:
        fig = plt.figure(figsize=(7, 7))
        ax = sns.boxplot(data=data, x='species', y=T_or_R, order=['Pipistrellus', 'Rousettus'], palette=my_pal,
                         width=0.7, showfliers=False)
        test_results = add_stat_annotation(ax, data=data, x='species', y=T_or_R,
                                           box_pairs=[("Rousettus", "Pipistrellus")],
                                           test='Mann-Whitney', text_format='star',
                                           loc='outside', verbose=2)
        sns.despine(offset=10, trim=False, fig=fig)
        plt.tick_params(axis='x', which='both', labelsize=14, labelrotation=45)
        plt.ylabel(T_or_R, fontsize=20)
        plt.tick_params(axis='y', labelsize=14)

        # save results
        plt.savefig(r'{}\{}_frac-{}.pdf'.format(results_dir, T_or_R, frac))
        print(f'tolerance_and_resistance_in_2_bats for {name_DE_genes} - DONE!')

# ---------------------------------------------------------------------------
# Analysis#4: Gene expression in comparative tissues of species
# ---------------------------------------------------------------------------

def normalize_and_transform_tissues_gene_expression(tissues_df, norm_or_not):
    """
    # This function normalizes and log transforms data as followed:
    # if norm_or_not==True each column values are divided by the total sum of the column and multiplied by 10^6, then log transformed

    :param tissues_df: pandas dataframe with numeric values.
    :param norm_or_not: Boolean.
    # =False - no normalization, only log transformation
    # =True - perform the normalization before log
    :return: normalized processed pandas dataframe
    """

    trans_tissues_df = tissues_df.dropna()
    # normalization
    if norm_or_not == True:
        trans_tissues_df = tissues_df.div(tissues_df.sum(axis=0))
        trans_tissues_df = trans_tissues_df * 1000000

    # log transformation
    trans_tissues_df = np.log10(trans_tissues_df)
    trans_tissues_df.replace([np.inf, -np.inf], 0, inplace=True)

    return trans_tissues_df


def tissues_basal_expression_heatmaps(tissues_df, FC_df, genes_list, genes_list_name,
                                      ttests=[('bat', 'human'), ('bat', 'mouse')], norm_or_not=True,
                                      res_dir_name='basal_expression_in_tissues'):
    """
    # This function generates heatmaps of basal level-normalized log(TPM) values in comparative tissues & FC values in the analysed species
    # for each gene in the set of genes provided in 'genes_list'
    # For each two species provided in 'ttests' param paired-T-test is conducted between tissues normalized log(TPM) values.

    :param tissues_df: pandas dataframe of TPM values. containing genes names as rows (must be 1-1 orthologs genes between the analysed species),
    # columns as comparative tissues names in different species.
    # NOTE: - columns names are in format of '{tissues name}_{species name}'
    #       - should contain data on at least 2 species and 2 tissues
    #       - the tissues are comparable between species- meaning each species in the analysis should have data on each tissue
    :param FC_df: pandas dataframe of FC values. containing genes as rows (match the gene names in tissues_df) and columns with speceis names.
    # NOTE: - the species names here should match the species names in tissues_df
    #       - FC values are taken from the DE analysis between control and stimulation with dsRNA
    :param genes_list: list of strings. names of genes to analyse. should match the gene names in tissues_df and FC_df.
    :param genes_list_name: string. name of group of genes to analyse, for the results path.
    :param ttests: list of tuples (or lists). each tuple contains two species names to compare in paired-T-test.
    # should match the species names in tissues_df and FC_df.
    # by default- [('bat', 'human'), ('bat', 'mouse')]
    :param norm_or_not: Boolean.
    # =False - no normalization, only log transformation
    # =True - perform the normalization before log
    :param res_dir_name: string. name of the results directory.
    :return:
    # - csv file containing all paired t-test results according to 'ttests' input.
    #       for each pair of species:
            - p-value result of one sided tests for both species1-tissues-log(TPM)-values > species2-tissues-log(TPM)-values,
    #       and the other way around
    #       - FDR correction for chosen p-values (which are the minimum between the 2 paired t-tests for each analyzed gene)
    # - plot's pdf file saved to results directory
    """

    # pre-processing
    tissues_df = tissues_df.reindex(sorted(tissues_df.columns), axis=1)
    FC_df = FC_df.reindex(sorted(FC_df.columns), axis=1)

    trans_tissues_df = normalize_and_transform_tissues_gene_expression(tissues_df, norm_or_not=norm_or_not) # normalizing tpm and log transform

    FC_df = FC_df.loc[(FC_df.index).intersection(trans_tissues_df.index)]

    FC_df = FC_df.loc[set(genes_list).intersection(FC_df.index)]
    trans_tissues_df = trans_tissues_df.loc[set(genes_list).intersection(trans_tissues_df.index)]

    # preparing dfs results
    ttest_df_results = pd.DataFrame()
    p_valuse_fo_FDR = defaultdict(dict)

    # plot results
    mycmap = sns.diverging_palette(200, 320, s=40, as_cmap=True)  # create cmap for heatmaps

    figsize = (8, len(genes_list) * 2)
    fig, axes = plt.subplots(len(genes_list), 2, sharex=False, figsize=figsize,
                             gridspec_kw={'width_ratios': [9, 1], 'wspace': 0.1, 'hspace': 0.75})

    for j in range(len(genes_list)):  # for each gene:
        # prepare tissues data for heatmap (rows as species, columns as tissues)
        df_per_gene = trans_tissues_df.loc[[genes_list[j]]]
        df_per_gene = df_per_gene.melt()
        df_per_gene[['tissue', 'species']] = df_per_gene.variable.str.split('_', expand=True)
        df_per_gene = df_per_gene.pivot_table(values='value', index='species', columns='tissue')
        # plot TPM tissues heatmap
        map0 = sns.heatmap(ax=axes[j][0], data=df_per_gene, cmap=mycmap, center=1, annot=False, square=True)
        map0.set_yticklabels(map0.get_ymajorticklabels(), rotation=0)
        map0.set_xticklabels(map0.get_xmajorticklabels(), rotation=90)
        map0.set_ylabel(genes_list[j], fontsize=14)
        map0.set_xlabel('')

        # prepare FC data for heatmap (rows as species, column-FC)
        df_per_gene_FC = FC_df.loc[[genes_list[j]]]
        df_per_gene_FC = df_per_gene_FC.T
        # plot FC heatmap
        map1 = sns.heatmap(ax=axes[j][1], data=df_per_gene_FC, cmap=mycmap, center=0, annot=False, square=True, vmin=-1,
                           xticklabels=False)
        map1.set_ylabel('')
        map1.set_xlabel('')
        map1.set_yticklabels(map1.get_ymajorticklabels(), rotation=0)

        # paired t-tests for gene for each pair in 'ttests' input
        gene = genes_list[j]
        for pair in ttests:
            species1 = pair[0]
            species2 = pair[1]
            species1_tpms = df_per_gene.loc[species1]
            species2_tpms = df_per_gene.loc[species2]

            # one sided paired-t-tests - both sides
            s1_greater_s2_pval = scipy.stats.ttest_rel(species1_tpms, species2_tpms, alternative='greater').pvalue
            s1_less_s2_pval = scipy.stats.ttest_rel(species1_tpms, species2_tpms, alternative='less').pvalue
            ttest_df_results.at[gene, f'{species1}_greater_{species2}'] = s1_greater_s2_pval
            ttest_df_results.at[gene, f'{species1}_less_{species2}'] = s1_less_s2_pval
            # get the minimum p-value (from both paired t-tests)
            min_pval = min(s1_greater_s2_pval, s1_less_s2_pval)
            p_valuse_fo_FDR[f'{species1}_vs_{species2}'][
                gene] = min_pval

    # organize to ttest_df_results
    ttest_df_results[''] = [''] * len(genes_list)
    for pair in ttests:
        species1 = pair[0]
        species2 = pair[1]
        pair_name = f'{species1}_vs_{species2}'
        genes = p_valuse_fo_FDR[pair_name].keys()
        ttest_df_results[f'gene_name_{pair_name}'] = genes
        pvalues = p_valuse_fo_FDR[pair_name].values()
        ttest_df_results[f'min_pval_{pair_name}'] = pvalues
        cors_pvals = sm.multipletests(method='fdr_bh', pvals=list(pvalues))[1]  # FDR correction for p-values
        ttest_df_results[f'FDR_{pair_name}'] = cors_pvals

    # save results
    if norm_or_not == True:
        normalization = 'norm'
    else:
        normalization = 'not_norm'
    results_dir = r'results\{}\{}'.format(res_dir_name, genes_list_name)
    os.makedirs(results_dir, exist_ok=True)
    plt.savefig(r'{}\{}_tpm_and_FC_heatmap.pdf'.format(results_dir, normalization))
    ttest_df_results.to_csv(r'{}\{}_tpm_paired_t_tests_results.csv'.format(results_dir, normalization))

    print(f'tissues_basal_expression_heatmaps for {genes_list_name} - DONE!')