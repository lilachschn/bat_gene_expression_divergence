# bat_gene_expression_divergence

This project explores the divergence of gene expression in response to stimulation to dsRNA and in basal level tissues, between two bats species - Rousettus aegyptiacus and Pipistrellus kuhlii and between rousettus and human and mouse. 

## Simply run main.py ! ! !

This project contains: 
- functions.py - functions related to the different analyses
- main.py - examples of the usage of the different analyses in functions.py. 
- results - directory with all results produced by running main.py (located in sub-directory for each analysis)
- Supporting_Tables.xlsx - file conaining the example data used in main.py.  
- 1-1_human_mouse_metadata.csv - file used for preprocessing the data in one of the analysis.
- DE_analysis_edgeR.R - DE analysis function written in R (optional - not used in main)
 

## DATA in Supporting_Tables.xlsx: 
---------------------------------
### Table S7: 
full Rousettus genes data with features values
- 'P_val_gene_gain_and_loss_ENSEMBL':    
The significance (P-value) of a gene family to have undergone a high rate of gene duplication and loss (contraction) over the course of vertebrate evolution, in comparison to other gene families, was obtained from ENSEMBL.
- 'pnas_dNdS_P':    
The dN/dS ratio (non-synonymous to synonymous codon substitutions) values of Rousettus genes across bats and their statistical significance were obtained from a previous studies that used orthologous genes from 18 bat species (Hawkins et al., 2019)
- 'seq_similarity_rous-pip_per_identity':     
level of similarity between orthologs obtained by BLAST of Rousettus versus Pipistrellus. selected the longest CDS for each pair of one-to-one orthologs (and the lowest TSL, in case two transcripts had the same length) of the two species and ran blast2seq on these the CDS sequences. 
- DM values:    
cell-to-cell variability of genes, used the DM (Distance to Median) approach.

### Table S9: 
data on 1-1 ortohlogs genes between the 6 species: Huamn, Macaque, Mouse, Rat, Rousettus, Pipistrellus
- gene expression-tpm data for 8 comparable tissues between human, mouse and rousettus-
Obtained from Gtex, BodyMap and Pavlovich,2020 respectively, were mapped and quantified. 
- DE analysis resutls (FC, p-value, q-value) for the 6 species.

### Table S16:
- tolerance and resistance values:   
Genes primarily associated with tolerance and resistance expression programs were taken from a previous study, which compared the transcriptional response and physiological outcomes of infection across a large set of mouse strains. There are 2 measures for Resistance and 2 for Tolerance. 


### Obtaining quantification of gene expression and DE analysis results:
----------------------------------------------------------------------
- Mapping:   
Reads were mapped and gene expression was quantified using Salmon with the following command parameters for each sample: ‘salmon quant -i [index_file_directory] -l IU -1 [left_read_library/lane1] [left_read_library/lane2] -2 [right_read_library/lane1] [right_read_library/lane2] --geneMap [transcript_to_gene_file] --seqBias --gcBias -q --numBootstrap 100 --threads 8 --validateMappings -o [output_directory]’.   
Each sample was mapped to its respective species' annotated transcriptome
- DE analysis:   
differential gene expression between dsRNA-treatment and control (for each species separately) by using edgeR with rounded estimated counts from Salmon.   
DE analysis produced FC value, p-value and corrected p-values for multiple testing by estimating FDR (q-value).   
DE analysis was done using the edgeR analysis general function in DE_analysis_edgeR.R script.   


## Scripts: 
--------
- DE_analysis_edgeR.R - a general function for DE analysis using edgeR package.   
USAGE:   
edgeR_analysis(samples_info, count_values, pair, out_filename) 


### functions.py
- scatter_correlation:  
USAGE:  
scatter_correlation(species1_fc,species2_fc,df,is_FC,res_dir_name)
- splitting_to_3_divergent_groups:   
USAGE:  
results_df=(df, measure_of_divergence, absolute)
- feature_in_3_divergent_groups:   
USAGE:  
feature_in_3_divergent_groups(feature, df, mann_whitney_side, measure_of_divergence, absolute, res_dir_name)
- tolerance_and_resistance_in_2_bats:   
USAGE:  
tolerance_and_resistance_in_2_bats(df_DE_gene, name_DE_genes, frac, res_dir_name)
- normalize_and_transform_tissues_gene_expression:   
USAGE:  
normalize_and_transform_tissues_gene_expression(tissues_df, norm_or_not)
- tissues_basal_expression_heatmaps:   
USAGE:   
tissues_basal_expression_heatmaps(tissues_df, FC_df, genes_list, genes_list_name, ttests, norm_or_not, res_dir_name)

### main.py
explanation in file.
