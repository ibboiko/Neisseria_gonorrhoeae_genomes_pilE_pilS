README_SmaCla_Pos_95Sens_pilE_dnaA_RScript.R

This R script analyzes the genomic positions of SmaCla repetitive sequences and other target sequences (dnaA, garP, and the Leader sequences of pilE) in Neisseria gonorrhoeae genomes. The script detects SmaCla and pilE-associated SmaCla repeats near Leader sequences.

The script requires standard GenBank files (.gb format) containing genomic sequences.

Input folder: Finished_genome_N65/ (located in data/raw/ folder); you need to download the tested finished genomes following the Manuscript (see section data availability at https://prism.northwestern.edu/).

Script handles special reference genomes FA1090 and MS11.

Install required missing R packages (if any) using: install.packages(c("Biostrings", "openxlsx", "stringr"))

Analysis Steps
Reads GenBank files and extracts sequences
Searches for target sequences (dnaA, SmaCla, garP, Leader_sequence_Class_I) using approximate pattern matching
Uses a 95% similarity threshold for SmaCla detection and 80% for other targets
Identifies pilE-associated SmaCla repeats within 1000bp of Leader sequences
Calculates relative positions to dnaA start site for the genomic context
Handles dnaA detection with trimmed sequences if the complete sequence is not found
Creates output file v08_SmaCla_Position_Output_N65.xlsx containing:
Genome ID, length, and feature positions for all analyzed genomes
Strand orientation (forward/reverse) and similarity scores
Relative positions calculated from the dnaA start site
Notes on trimmed sequences and distance measurements
pilE-associated SmaCla repeat positions with distance to Leader sequences

Run this script from R Studio.
