#!/usr/bin/env Rscript

# =====================================================================
# SmaCla Position Analysis Script (Version 03) - Statistics Only
# =====================================================================
# This script analyzes SmaCla element distribution and performs statistical
# analysis of SmaCla elements between genomes with and without pilE genes.
# Output includes only statistical tables, no figures.
# =====================================================================

# Load required libraries
library(readxl)        # For reading Excel files
library(dplyr)         # For data manipulation
library(writexl)       # For writing Excel files
library(tidyr)         # For data tidying
library(car)           # For statistical tests (Levene's test)
library(effsize)       # For calculating effect sizes

# =====================================================================
# Define input and output paths
# =====================================================================
input_file <- "./v08_SmaCla_Position_Output_N65.xlsx"
output_folder <- "./SmaCla_Stat_Output"

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# =====================================================================
# Read data
# =====================================================================
# Read the Excel file
cat("Reading input file:", input_file, "\n")
data <- read_excel(input_file)

# Display basic information about the dataset
cat("\nBasic dataset information:\n")
cat("Total number of rows:", nrow(data), "\n")
cat("Number of unique genomes:", length(unique(data$Finished_genome_ID)), "\n")

# =====================================================================
# Task 1: Descriptive statistics of SmaCla
# =====================================================================

# Remove reference genomes (FA1090 and MS11)
filtered_data <- data %>%
  filter(!Finished_genome_ID %in% c("FA1090", "MS11"))

cat("\nAfter excluding reference genomes:\n")
cat("Number of genomes:", length(unique(filtered_data$Finished_genome_ID)), "\n")

# Create a dataframe for SmaCla counts per genome
smacla_counts <- filtered_data %>%
  filter(Feature == "SmaCla") %>%
  group_by(Finished_genome_ID) %>%
  summarize(
    SmaCla_Count = n(),
    Genome_Length = Genome_Length[1]  # Get first genome length value
  )

# Calculate descriptive statistics for SmaCla counts
smacla_stats <- data.frame(
  Statistic = c("Min", "1st Quartile", "Median", "Mean", "3rd Quartile", "Max", "SD", "N"),
  Value = c(
    min(smacla_counts$SmaCla_Count),
    quantile(smacla_counts$SmaCla_Count, 0.25),
    median(smacla_counts$SmaCla_Count),
    mean(smacla_counts$SmaCla_Count),
    quantile(smacla_counts$SmaCla_Count, 0.75),
    max(smacla_counts$SmaCla_Count),
    sd(smacla_counts$SmaCla_Count),
    length(smacla_counts$SmaCla_Count)
  )
)

cat("\nSmaCla Count Statistics:\n")
print(smacla_stats)

# =====================================================================
# Task 2: Identify genomes with and without pilE gene
# =====================================================================

# Identify genomes with pilE gene (both "Leader_sequence_Class_I" and "garP" are found)
genome_features <- filtered_data %>%
  group_by(Finished_genome_ID) %>%
  summarize(
    has_leader = any(Feature == "Leader_sequence_Class_I" & 
                       !grepl("not found", Position, fixed = TRUE)),
    has_garP = any(Feature == "garP" & 
                     !grepl("not found", Position, fixed = TRUE))
  ) %>%
  mutate(has_pilE = has_leader & has_garP)

# Create lists of genomes with and without pilE
genomes_with_pilE <- genome_features %>%
  filter(has_pilE) %>%
  select(Finished_genome_ID)

genomes_without_pilE <- genome_features %>%
  filter(!has_pilE) %>%
  select(Finished_genome_ID)

cat("\nGenomes with pilE gene:", nrow(genomes_with_pilE), "\n")
cat("Genomes without pilE gene:", nrow(genomes_without_pilE), "\n")

# =====================================================================
# Statistical analysis: SmaCla differences between genomes with/without pilE
# =====================================================================

# Add pilE status to the SmaCla counts
smacla_counts_with_status <- smacla_counts %>%
  left_join(
    genome_features %>% select(Finished_genome_ID, has_pilE),
    by = "Finished_genome_ID"
  ) %>%
  mutate(pilE_status = factor(ifelse(has_pilE, "With pilE", "Without pilE")))

# Calculate descriptive statistics by pilE status
smacla_stats_by_pilE <- smacla_counts_with_status %>%
  group_by(pilE_status) %>%
  summarize(
    n = n(),
    min = min(SmaCla_Count),
    q1 = quantile(SmaCla_Count, 0.25),
    median = median(SmaCla_Count),
    mean = mean(SmaCla_Count),
    q3 = quantile(SmaCla_Count, 0.75),
    max = max(SmaCla_Count),
    sd = sd(SmaCla_Count)
  )

cat("\nSmaCla Count Statistics by pilE Status:\n")
print(smacla_stats_by_pilE)

# Test for normality (Shapiro-Wilk test)
shapiro_with_pilE <- shapiro.test(
  smacla_counts_with_status$SmaCla_Count[smacla_counts_with_status$has_pilE]
)
shapiro_without_pilE <- shapiro.test(
  smacla_counts_with_status$SmaCla_Count[!smacla_counts_with_status$has_pilE]
)

cat("\nNormality tests:\n")
cat("With pilE: W =", shapiro_with_pilE$statistic, "p-value =", shapiro_with_pilE$p.value, "\n")
cat("Without pilE: W =", shapiro_without_pilE$statistic, "p-value =", shapiro_without_pilE$p.value, "\n")

# Test for homogeneity of variance (Levene's test)
levene_result <- leveneTest(SmaCla_Count ~ pilE_status, data = smacla_counts_with_status)
cat("\nLevene's test for homogeneity of variance:\n")
print(levene_result)

# Perform appropriate statistical tests
t_test_result <- t.test(SmaCla_Count ~ pilE_status, data = smacla_counts_with_status)
wilcox_result <- wilcox.test(SmaCla_Count ~ pilE_status, data = smacla_counts_with_status)

cat("\nGroup comparison tests:\n")
cat("t-test result:\n")
print(t_test_result)
cat("\nWilcoxon rank-sum test result:\n")
print(wilcox_result)

# Calculate effect size (Cohen's d)
cohen_d_result <- cohen.d(smacla_counts_with_status$SmaCla_Count, 
                          smacla_counts_with_status$pilE_status)
cat("\nEffect size (Cohen's d):", cohen_d_result$estimate, 
    "Interpretation:", cohen_d_result$magnitude, "\n")

# =====================================================================
# Analyze genomic positions of SmaCla
# =====================================================================

# Gather position data for SmaCla elements
smacla_positions <- filtered_data %>%
  filter(Feature == "SmaCla") %>%
  select(Finished_genome_ID, Relative_position_to_dnaA, Strand) %>%
  mutate(Relative_position_to_dnaA = as.numeric(as.character(Relative_position_to_dnaA))) %>%
  filter(!is.na(Relative_position_to_dnaA)) %>%
  left_join(
    genome_features %>% select(Finished_genome_ID, has_pilE),
    by = "Finished_genome_ID"
  ) %>%
  mutate(pilE_status = factor(ifelse(has_pilE, "With pilE", "Without pilE")))

# Helper function to safely calculate statistics
safe_stats <- function(x) {
  if(length(x) == 0 || all(is.na(x))) {
    return(list(
      min = NA_real_,
      q1 = NA_real_,
      median = NA_real_,
      mean = NA_real_,
      q3 = NA_real_,
      max = NA_real_,
      sd = NA_real_
    ))
  }
  
  x <- x[!is.na(x)]
  
  list(
    min = if(length(x) > 0) min(x) else NA_real_,
    q1 = if(length(x) > 0) as.numeric(quantile(x, 0.25)) else NA_real_,
    median = if(length(x) > 0) median(x) else NA_real_,
    mean = if(length(x) > 0) mean(x) else NA_real_,
    q3 = if(length(x) > 0) as.numeric(quantile(x, 0.75)) else NA_real_,
    max = if(length(x) > 0) max(x) else NA_real_,
    sd = if(length(x) > 1) sd(x) else NA_real_
  )
}

# Calculate position statistics by pilE status
position_stats_by_pilE <- smacla_positions %>%
  group_by(pilE_status) %>%
  summarize(
    n = sum(!is.na(Relative_position_to_dnaA)),
    min_pos = safe_stats(Relative_position_to_dnaA)$min,
    q1_pos = safe_stats(Relative_position_to_dnaA)$q1,
    median_pos = safe_stats(Relative_position_to_dnaA)$median,
    mean_pos = safe_stats(Relative_position_to_dnaA)$mean,
    q3_pos = safe_stats(Relative_position_to_dnaA)$q3,
    max_pos = safe_stats(Relative_position_to_dnaA)$max,
    sd_pos = safe_stats(Relative_position_to_dnaA)$sd
  )

cat("\nSmaCla Position Statistics by pilE Status:\n")
print(position_stats_by_pilE)

# Test for normality of position data
shapiro_pos_with_pilE <- tryCatch({
  positions_with_pilE <- smacla_positions$Relative_position_to_dnaA[smacla_positions$has_pilE]
  if(length(positions_with_pilE) >= 3) {
    shapiro.test(positions_with_pilE)
  } else {
    list(statistic = NA, p.value = NA)
  }
}, error = function(e) {
  list(statistic = NA, p.value = NA)
})

shapiro_pos_without_pilE <- tryCatch({
  positions_without_pilE <- smacla_positions$Relative_position_to_dnaA[!smacla_positions$has_pilE]
  if(length(positions_without_pilE) >= 3) {
    shapiro.test(positions_without_pilE)
  } else {
    list(statistic = NA, p.value = NA)
  }
}, error = function(e) {
  list(statistic = NA, p.value = NA)
})

cat("\nNormality tests for positions:\n")
cat("With pilE: W =", shapiro_pos_with_pilE$statistic, "p-value =", shapiro_pos_with_pilE$p.value, "\n")
cat("Without pilE: W =", shapiro_pos_without_pilE$statistic, "p-value =", shapiro_pos_without_pilE$p.value, "\n")

# Test for position differences
ks_test_result <- tryCatch({
  positions_with_pilE <- smacla_positions$Relative_position_to_dnaA[smacla_positions$has_pilE]
  positions_without_pilE <- smacla_positions$Relative_position_to_dnaA[!smacla_positions$has_pilE]
  
  if(length(positions_with_pilE) > 0 && length(positions_without_pilE) > 0) {
    ks.test(positions_with_pilE, positions_without_pilE)
  } else {
    list(statistic = NA, p.value = NA)
  }
}, error = function(e) {
  list(statistic = NA, p.value = NA)
})

cat("\nKolmogorov-Smirnov test for position distributions:\n")
print(ks_test_result)

# Strand analysis
strand_analysis <- smacla_positions %>%
  mutate(Strand = as.character(Strand)) %>%
  mutate(Strand = ifelse(is.na(Strand), "Unknown", Strand)) %>%
  group_by(pilE_status, Strand) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(pilE_status) %>%
  mutate(percentage = count / sum(count) * 100)

cat("\nStrand distribution by pilE status:\n")
print(strand_analysis)

# Chi-square test for strand distribution
strand_contingency <- table(smacla_positions$pilE_status, smacla_positions$Strand)
chisq_result <- tryCatch({
  if(sum(strand_contingency) > 0 && min(strand_contingency) > 0) {
    chisq.test(strand_contingency)
  } else {
    fisher_result <- fisher.test(strand_contingency)
    list(statistic = NA, p.value = fisher_result$p.value)
  }
}, error = function(e) {
  tryCatch({
    fisher_result <- fisher.test(strand_contingency)
    list(statistic = NA, p.value = fisher_result$p.value)
  }, error = function(e2) {
    list(statistic = NA, p.value = NA)
  })
})

cat("\nChi-square test for strand distribution:\n")
print(chisq_result)

# =====================================================================
# Create output files
# =====================================================================

# Create a list to hold all the data for the output file
smacla_stat_sheets <- list(
  "Overall_Stats" = smacla_stats,
  "Stats_By_pilE" = smacla_stats_by_pilE,
  "Position_Stats" = position_stats_by_pilE,
  "Statistical_Tests" = data.frame(
    Test = c("Shapiro-Wilk (With pilE)", "Shapiro-Wilk (Without pilE)", 
             "Levene's Test", "t-test", "Wilcoxon Test", "Cohen's d",
             "KS Test (Positions)", "Chi-square (Strand)"),
    Statistic = c(shapiro_with_pilE$statistic, shapiro_without_pilE$statistic,
                  levene_result$`F value`[1], t_test_result$statistic,
                  wilcox_result$statistic, cohen_d_result$estimate,
                  ks_test_result$statistic, chisq_result$statistic),
    p_value = c(shapiro_with_pilE$p.value, shapiro_without_pilE$p.value,
                levene_result$`Pr(>F)`[1], t_test_result$p.value,
                wilcox_result$p.value, NA,
                ks_test_result$p.value, chisq_result$p.value),
    Interpretation = c(
      ifelse(shapiro_with_pilE$p.value < 0.05, "Non-normal", "Normal"),
      ifelse(shapiro_without_pilE$p.value < 0.05, "Non-normal", "Normal"),
      ifelse(levene_result$`Pr(>F)`[1] < 0.05, "Unequal variance", "Equal variance"),
      ifelse(t_test_result$p.value < 0.05, "Significant difference", "No significant difference"),
      ifelse(wilcox_result$p.value < 0.05, "Significant difference", "No significant difference"),
      cohen_d_result$magnitude,
      ifelse(ks_test_result$p.value < 0.05, "Different distributions", "Similar distributions"),
      ifelse(chisq_result$p.value < 0.05, "Different strand distributions", "Similar strand distributions")
    )
  ),
  "SmaCla_Counts_Per_Genome" = smacla_counts_with_status,
  "Strand_Analysis" = strand_analysis
)

# Write the output file
write_xlsx(smacla_stat_sheets, 
           path = file.path(output_folder, "SmaCla_stat.xlsx"))
cat("\nCreated output file:", file.path(output_folder, "SmaCla_stat.xlsx"), "\n")

# Write the genomes with pilE gene
write_xlsx(list("Genomes_with_pilE" = genomes_with_pilE), 
           path = file.path(output_folder, "Genomes_ID_with_pilE_gene.xlsx"))
cat("Created output file:", file.path(output_folder, "Genomes_ID_with_pilE_gene.xlsx"), "\n")

# Write the genomes without pilE gene
write_xlsx(list("Genomes_without_pilE" = genomes_without_pilE), 
           path = file.path(output_folder, "Genomes_ID_without_pilE_gene.xlsx"))
cat("Created output file:", file.path(output_folder, "Genomes_ID_without_pilE_gene.xlsx"), "\n")

cat("\nAnalysis complete! Only statistical outputs were generated.\n")
