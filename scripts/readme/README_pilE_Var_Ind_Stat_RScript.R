README_pilE_Var_Ind_Stat_RScript.R

This R script performs statistical analysis and visualization of pilE gene variation indices across its different structural components in Neisseria gonorrhoeae.
The script creates boxplots with pairwise statistical comparisons and generates statistical reports.

The script requires an Excel file containing the pilE variation index data:
Input file: v06_Ntds_pilE_VarIndex_data.xlsx (located in data/raw/ folder);
The file should contain columns: Reference_Position and Variation_Index;
The data represent amino acid variation indices across pilE gene positions.

Install required missing R packages (if any) using: install.packages(c("ggplot2", "gridExtra", "grid", "gtable", "openxlsx", "readxl", "dplyr", "rstatix", "ggpubr", "stats"))

Analysis Steps
Reads pilE variation index data from an Excel file;
Categorizes amino acid positions into structural components (Conserved, Semi-variable, cys1, Hypervariable loop, cys2, Hypervariable tail);
Filters groups with sufficient data (minimum 3 observations) for statistical analysis;
Performs descriptive statistics (mean, median, SD, IQR, confidence intervals);
Conducts Wilcoxon rank-sum tests for pairwise comparisons between structural components;
Applies Bonferroni correction for multiple testing adjustment;
Tests the data normality using the Shapiro-Wilk test.

Creates boxplots with statistical significance annotations

Output Files
Creates Output_AB_Ntds_pilE_Var_Ind_Stat/ folder containing:
v26_Ntds_pilE_Var_Ind_Stat.tif: High-resolution boxplot (900 DPI, TIFF format);
v26_Ntds_pilE_Var_Ind_Stat.png: Boxplot in PNG format;
v26_Ntds_pilE_Var_Ind_Stat.jpg: Boxplot in JPEG format;
v26_Ntds_pilE_Var_Ind_Stat_statistics.xlsx: Statistical results with multiple sheets including basic statistics, pairwise comparisons, normality tests, and statistical methods documentation.

Run this script from R Studio.
