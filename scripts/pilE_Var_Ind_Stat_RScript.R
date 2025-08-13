# Load required libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(openxlsx)
library(readxl)
library(dplyr)
library(rstatix)
library(ggpubr)
library(stats)

###################
## File Paths
###################
stat_input <- "./v06_Ntds_pilE_VarIndex_data.xlsx"
output_dir <- "./Output_AB_Ntds_pilE_Var_Ind_Stat"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define output files
base_filename <- "v26_Ntds_pilE_Var_Ind_Stat"
tif_output <- file.path(output_dir, paste0(base_filename, ".tif"))
png_output <- file.path(output_dir, paste0(base_filename, ".png"))
jpg_output <- file.path(output_dir, paste0(base_filename, ".jpg"))
xlsx_output <- file.path(output_dir, paste0(base_filename, "_statistics.xlsx"))

###################
## Common Parameters
###################
COMMON_TEXT_SIZE <- 8
COMMON_AXIS_TEXT_SIZE <- 11
COMMON_AXIS_TITLE_SIZE <- 11
DOT_SIZE <- 1.5
DOT_ALPHA <- 0.6
PLOT_MARGIN <- unit(c(2, 2, 2, 2), "pt")
PANEL_SPACING <- unit(2, "pt")  # Exactly 2px as per journal requirement
BAR_WIDTH <- 0.6

# Part labels and orders
PART_ORDER <- c("Conserved", "Semi-variable", "cys1", "Hypervariable loop", "cys2", "Hypervariable tail")
PART_LABELS <- c("C", "SV", expression(italic("cys1")), "HVL", expression(italic("cys2")), "HVT")
names(PART_LABELS) <- PART_ORDER

###################
## Common Theme
###################
common_theme <- theme_minimal() +
  theme(
    text = element_text(size = COMMON_TEXT_SIZE, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
    axis.text = element_text(color = "black", size = COMMON_AXIS_TEXT_SIZE),
    axis.title = element_text(color = "black", size = COMMON_AXIS_TITLE_SIZE),
    axis.title.x = element_text(margin = margin(t = 6)),  # Reduced space above x-axis title
    plot.background = element_rect(fill = "white", color = NA)
  )

###################
## Statistical Analysis Functions
###################
calculate_statistics <- function(data) {
  # Filter out parts with insufficient data
  valid_parts <- data %>%
    group_by(Part) %>%
    filter(n() >= 3) %>%
    pull(Part) %>%
    unique()
  
  filtered_data <- data %>% filter(Part %in% valid_parts)
  
  # Basic statistics
  basic_stats <- filtered_data %>%
    group_by(Part) %>%
    summarise(
      N = n(),
      Mean = mean(Variation_Index, na.rm = TRUE),
      SD = sd(Variation_Index, na.rm = TRUE),
      Median = median(Variation_Index, na.rm = TRUE),
      Min = min(Variation_Index, na.rm = TRUE),
      Max = max(Variation_Index, na.rm = TRUE),
      IQR = IQR(Variation_Index, na.rm = TRUE),
      CI_lower = ifelse(n() >= 3, t.test(Variation_Index)$conf.int[1], NA),
      CI_upper = ifelse(n() >= 3, t.test(Variation_Index)$conf.int[2], NA)
    ) %>%
    ungroup()
  
  # Pairwise comparisons only if we have at least 2 groups with sufficient data
  if(length(valid_parts) >= 2) {
    stat_test <- filtered_data %>%
      wilcox_test(
        Variation_Index ~ Part,
        p.adjust.method = "bonferroni",
        detailed = TRUE
      ) %>%
      add_significance("p") %>%
      mutate(
        Significance = case_when(
          p < 0.0001 ~ "****",
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
  } else {
    stat_test <- data.frame(Note = "Insufficient data for pairwise comparisons")
  }
  
  # Normality test
  normality_test <- filtered_data %>%
    group_by(Part) %>%
    summarise(
      shapiro_stat = ifelse(n() >= 3, shapiro.test(Variation_Index)$statistic, NA),
      shapiro_p = ifelse(n() >= 3, shapiro.test(Variation_Index)$p.value, NA)
    )
  
  return(list(
    basic_stats = basic_stats,
    pairwise_tests = stat_test,
    normality_test = normality_test,
    filtered_parts = valid_parts
  ))
}

###################
## Panel Creation
###################
# Read and process data
stat_data <- read_excel(stat_input)

# Add Part column based on Reference_Position
stat_data$Part <- with(stat_data, 
                       case_when(
                         Reference_Position <= 150 ~ "Conserved",
                         Reference_Position <= 360 ~ "Semi-variable",
                         Reference_Position <= 399 ~ "cys1",
                         Reference_Position <= 444 ~ "Hypervariable loop",
                         Reference_Position <= 477 ~ "cys2",
                         TRUE ~ "Hypervariable tail"
                       )
)

# Convert Part to factor with correct order
stat_data$Part <- factor(stat_data$Part, levels = PART_ORDER)

# Get valid parts with sufficient data
valid_parts <- stat_data %>%
  group_by(Part) %>%
  filter(n() >= 3) %>%
  pull(Part) %>%
  unique()

# Define colors for boxplot
box_colors <- c(
  "Conserved" = "green4",
  "Semi-variable" = "blue",
  "cys1" = "purple",
  "Hypervariable loop" = "grey64",
  "cys2" = "purple",
  "Hypervariable tail" = "orange"
)

# Define comparisons only for valid parts
comparisons <- list(
  c("Conserved", "Hypervariable loop"),
  c("Conserved", "Hypervariable tail"),
  c("Semi-variable", "Hypervariable loop"),
  c("Semi-variable", "Hypervariable tail"),  
  c("cys1", "Hypervariable loop"),
  c("cys1", "Hypervariable tail"),
  c("cys2", "Hypervariable loop"),
  c("cys2", "Hypervariable tail")
) %>% purrr::keep(~ all(.x %in% valid_parts))

# Create main plot with adjusted significance bars
main_plot <- ggplot(stat_data, aes(x = Part, y = Variation_Index, fill = Part)) +
  geom_boxplot(
    width = BAR_WIDTH,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  scale_fill_manual(values = alpha(box_colors, 0.6)) +
  scale_x_discrete(labels = PART_LABELS) +
  scale_y_continuous(
    limits = c(0, 305),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  labs(
    x = expression(paste(italic("pilE"), " gene structural components")),
    y = expression(paste(italic("pilE"), " variation index, %"))
  ) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  ) +
  stat_compare_means(
    comparisons = comparisons,
    label = "p.signif",
    method = "wilcox.test",
    step.increase = 0.285,
    tip.length = 0.02,
    size = 3.5,
    vjust = 0
  )

###################
## Create Compact Legend (2 rows) - Modified Version
###################
legend_items <- list(
  labels = list("Conserved (C)", "Semi-variable (SV)", 
                "Hypervariable loop (HVL)", "Hypervariable tail (HVT)"),
  colors = c("#333333", "blue", "grey64", "orange"),
  x_pos = c(1.5, 2.7, 1.5, 2.7),  # Centered positions for 2 columns in 2 rows
  y_pos = c(1.8, 1.8, 1.0, 1.0)        # Two rows
)

# Create compact legend plot
legend_plot <- ggplot() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

# Add legend items with tight spacing
for (i in 1:4) {
  legend_plot <- legend_plot +
    annotate("point", x = legend_items$x_pos[i], y = legend_items$y_pos[i], 
             color = legend_items$colors[i], size = 2.0) +
    annotate("text", x = legend_items$x_pos[i], y = legend_items$y_pos[i] + 0.35,
             label = legend_items$labels[[i]], hjust = 0.5, size = 2.8,
             fontface = "plain")
}

# Set tight plot limits
legend_plot <- legend_plot +
  xlim(0.8, 3.2) +  # Adjusted for better centering
  ylim(0.7, 2.2)    # Tight vertical limits

###################
## Combine Plot and Legend
###################
final_plot <- arrangeGrob(
  main_plot,
  legend_plot,
  nrow = 2,
  heights = c(0.88, 0.12),  # More space for main plot, less for legend
  padding = unit(2, "pt")    # Minimal spacing
)

###################
## Save Statistical Results
###################
# Perform statistical analysis
stat_results <- calculate_statistics(stat_data)

# Create workbook
wb <- createWorkbook()

# Add basic statistics sheet
addWorksheet(wb, "Basic Statistics")
writeData(wb, "Basic Statistics", stat_results$basic_stats, rowNames = TRUE)

# Add pairwise comparisons sheet
addWorksheet(wb, "Pairwise Comparisons")
writeData(wb, "Pairwise Comparisons", stat_results$pairwise_tests, rowNames = TRUE)

# Add normality test sheet
addWorksheet(wb, "Normality Tests")
writeData(wb, "Normality Tests", stat_results$normality_test, rowNames = TRUE)

# Add methods description sheet
addWorksheet(wb, "Statistical Methods")
writeData(wb, "Statistical Methods", data.frame(
  Section = c("Basic Statistics", "Pairwise Comparisons", "Multiple Testing", "Significance Levels", "Normality Test"),
  Description = c(
    "Descriptive statistics including mean, median, standard deviation (SD), minimum, maximum, interquartile range (IQR), and 95% confidence intervals were calculated for each group.",
    "Wilcoxon rank-sum test was used for pairwise comparisons due to non-normal distribution of the data.",
    "Bonferroni correction was applied to adjust for multiple comparisons.",
    "Significance levels: **** p<0.0001, *** p<0.001, ** p<0.01, * p<0.05, ns p≥0.05",
    "Shapiro-Wilk test was used to assess normality of the data distribution in each group."
  )
))

# Save workbook
saveWorkbook(wb, xlsx_output, overwrite = TRUE)

###################
## Save Functions
###################
save_plot <- function(filename, plot, width, height, dpi = 900) {
  tryCatch({
    if (file.exists(filename)) file.remove(filename)
    
    if (grepl("\\.jpg$", filename)) {
      jpeg(filename, width = width * dpi, height = height * dpi,
          res = dpi, type = "cairo")
      grid.draw(plot)
      dev.off()
    } else if (grepl("\\.tif$", filename)) {
      tiff(filename, width = width * dpi, height = height * dpi, 
           res = dpi, compression = "lzw", type = "cairo")
      grid.draw(plot)
      dev.off()
    } else if (grepl("\\.png$", filename)) {
      png(filename, width = width * dpi, height = height * dpi, 
          res = dpi, type = "cairo")
      grid.draw(plot)
      dev.off()
    }
    cat("Successfully saved:", filename, "\n")
  }, error = function(e) {
    cat("Error saving", filename, ":", e$message, "\n")
  })
}

###################
## Save Outputs
###################
PLOT_WIDTH <- 3.425  # inches
PLOT_HEIGHT <- 4.2   # inches (adjusted to be more compact)

save_plot(tif_output, final_plot, PLOT_WIDTH, PLOT_HEIGHT)
save_plot(png_output, final_plot, PLOT_WIDTH, PLOT_HEIGHT)
save_plot(jpg_output, final_plot, PLOT_WIDTH, PLOT_HEIGHT)

cat("Analysis complete. Plot and statistics saved in:", output_dir, "\n")
