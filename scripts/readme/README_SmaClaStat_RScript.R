This R script performs statistical analysis of SmaCla repetitive sequences distribution in Neisseria gonorrhoeae genomes. The script compares SmaCla element counts and genomic positions between genomes with and without pilE genes.
The script requires an Excel file containing SmaCla position data:
Input file: v08_SmaCla_Position_Output_N65.xlsx (located in data/raw/ folder)
File should contain columns: Finished_genome_ID, Feature, Position, Relative_position_to_dnaA, Strand, Genome_Length
Install required missing R packages (if any) using: install.packages(c("readxl", "dplyr", "writexl", "tidyr", "car", "effsize"))
Analysis Steps
Reads SmaCla position data from Excel file
Excludes reference genomes (FA1090 and MS11) and counts only tested genomes
Calculates descriptive statistics for SmaCla counts per genome
Identifies genomes with and without pilE genes based on presence of Leader_sequence_Class_I and garP features
Performs statistical comparisons between groups using t-tests, Wilcoxon tests, and effect size calculations
Analyzes genomic position distributions using Kolmogorov-Smirnov tests
Examines strand distribution patterns with chi-square tests
Output Files
Creates SmaCla_Stat_Output/ folder containing:
SmaCla_stat.xlsx: Statistical results with multiple sheets
Genomes_ID_with_pilE_gene.xlsx: List of genomes containing pilE genes
Genomes_ID_without_pilE_gene.xlsx: List of genomes lacking pilE genes
Run this script from R Studio. The script will automatically create output directories and generate statistical summary files. All paths are relative to the script location.
