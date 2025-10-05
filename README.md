# Neisseria_gonorrhoeae_genomes_pilE_pilS
This repository contains scripts used in the Manuscript "Publicly available Neisseria gonorrhoeae genomes predominantly represent in vitro-derived nonpiliated variants".

The repository provides computational tools for analyzing SmaCla repetitive sequences and pilE gene/PilE protein variation in Neisseria gonorrhoeae genomes.

Scripts Description

1. SmaCla_Pos_95Sens_pilE_dnaA_RScript.R (SmaCla-dnaA relative distance  algorithm)
Analyzes genomic positions of SmaCla repetitive sequences and target sequences (dnaA, garP, Leader sequences of pilE). Identifies pilE-associated SmaCla repeats and calculates relative positions to dnaA start sites for genomic context analysis.
Use Protocol_PubMLST_in_silico.docx to upload finished genomes as input for Script 1 [Boiko, I., Metaane, S., & Seifert, H. S. (2025). Publicly available Neisseria gonorrhoeae genomes predominantly represent in vitro-derived nonpiliated variants (Version 2) [Data set]. Prism. Galter Health Sciences Library. Northwestern University. https://doi.org/10.18131/k42dy-w4a75]

3. SmaClaStat_RScript.R
Performs statistical analysis of SmaCla repetitive sequences distribution in N. gonorrhoeae genomes. Compares SmaCla element counts and genomic positions between genomes with and without pilE genes using statistical tests (t-tests, Wilcoxon tests, and effect size calculations).

4. generate_heatmap_with_gradient.py (PilE structure variation mapping algorithm)
Generates a 3D protein structure heatmap visualization of amino acid variation indices for the pilE protein. Creates PyMol-compatible scripts for protein structure visualization with color-coded variation patterns mapped to the protein structure.

5. pilE_Var_Ind_Stat_RScript.R
Performs statistical analysis and visualization of pilE gene variation indices across different structural components. Creates boxplots with pairwise statistical comparisons and generates statistical reports.

6.pilS_copies_generator_with_report
Identification and extraction pilS loci based on the recognition of the specific motifs cys2 and SmaCla

Data Requirements
Input files (located in data/raw/ folder):
v08_SmaCla_Position_Output_N65.xlsx: SmaCla position data for statistical analysis;
rv_02_pilE_Pt_PubMLST_variation_score.xlsx: Amino acid variation index data;
pilE_FA1090_Pt_b470a_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb: PDB structure file;
v06_Ntds_pilE_VarIndex_data.xlsx: pilE variation index data for statistical analysis.

Software Requirements

R Dependencies
install.packages(c("readxl", "dplyr", "writexl", "tidyr", "car", "effsize", 
                   "Biostrings", "openxlsx", "stringr", "ggplot2", "gridExtra", 
                   "grid", "gtable", "rstatix", "ggpubr", "stats"))

Python Dependencies
pip install pandas biopython matplotlib

See script-specific README files in the /scripts folder.

When using these scripts, please cite the  Manuscript: [Citation will be provided upon publishing the Manuscript]

This work is dual-licensed under the MIT License for software components and the Creative Commons Attribution 4.0 International License for research data and documentation. See the LICENSE file for details.

For questions regarding usage or licensing, please contact the Correspondence author.
