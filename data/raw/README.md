Raw Data
This folder contains the raw data files used in the analysis of Neisseria gonorrhoeae genomes, pilE, and pilS sequences.
File Descriptions
##

#Statistical Analysis of SmaCla Genomic Distribution
##File: v08_SmaCla_Position_Output_N65.xlsx
Input file for R scripts analyzing SmaCla genomic distribution
Used for Figure 3 (SmaCla repetitive sequences' relative genomic positions in the finished N. gonorrhoeae genomes, n = 65), Table 1 (SmaCla sequences analysis in N. gonorrhoeae finished genomes, n = 65)
Associated Script: v08_SmaClaStat_RScript.R (Contains statistical data for SmaCla genomic distribution analysis)

##
#pilE Nucleotide Variation Analysis
##File: v06_Ntds_pilE_VarIndex_data.xlsx
Input data for the pilE nucleotide variation and statistics. pilE gene sequence variation indices manually curated from multiple sequence alignments using Geneious Prime
Used for Figure 4 (Nucleotide variation of the pilE gene in N. gonorrhoeae genomes, n = 1558)
Associated Script: v26_pilE_Var_Ind_Stat_RScript_single_panel.R

##
#PilE Structure Variation Mapping
##Raw file 1: rv_02_pilE_Pt_PubMLST_variation_score.xlsx
Source data for the PilE Structure Variation Mapping Algorithm. PilE sequence variation indices manually curated from multiple sequence alignments using Geneious Prime software
Used for Figure 5 (PilE structure of N. gonorrhoeae FA1090 showing sequence variability and major truncation sites)
Associated Script: generate_heatmap_with_gradient_v09_5.py
##Raw file 2: pilE_FA1090_Pt_b470a_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb
Used for a Three-dimensional structural model for the PilE Structure Variation Mapping Algorithm, Figure 5
FA1090 PilE structure predicted using ColabFold 1.5.5 as described by Kim G, Lee S, Levy Karin E, et al. Easy and accurate protein structure prediction using ColabFold. Nat Protoc. 2025; 20(3):620–642.
Associated Script: generate_heatmap_with_gradient_v09_5.py

#
All scripts for processing this raw data are located in the main repository folder:
v08_SmaClaStat_RScript.R - SmaCla statistical analysis
v08_RScript_SmaCla_Pos_95Sens_pilE_dnaA.R - SmaCla-dnaA relative distance  algorithm
v26_pilE_Var_Ind_Stat_RScript_single_panel.R - pilE variation statistics
generate_heatmap_with_gradient_v09_5.py - PilE structure variation mapping algorithm

Notes
All files are in their original format as used for the published analysis.
