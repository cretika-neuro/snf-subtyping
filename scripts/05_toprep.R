####### METASNF: Step 5 ######
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metasnf)
library(SNFtool)
library(grid)
library(BiocManager)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

rm(list = ls())

# Set working directory here (left blank on purpose): <project root>/snf/

run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)
select <- dplyr::select

# Load R objects
input_dl <- readRDS(file.path(run, "R_vars/01_input_dl.rds"))

my_sconfig <- readRDS(file.path(run, "R_vars/02_my_sconfig.rds"))     # SNF config file for all solutions
sol_df <- readRDS(file.path(run, "R_vars", "02_SNF_soldf.rds"))       # Participant cluster assignment for all solutions
annotation_data <- readRDS(file.path(run, "R_vars/03_annotation_data.rds"))

pairwise_aris <- readRDS(file.path(run, "R_vars/04_pairwise_aris.rds"))
rep_solutions <- readRDS(file.path(run, "R_vars/03_rep_solutions.rds")) # Participant cluster assignment per representative solution
rep_data <- readRDS(file.path(run, "R_vars/03_rep_sil.rds"))          # Additional information for each representative solution
sol_rep <- readRDS(file.path(run, "R_vars/04_sol_rep.rds"))

coclustering_results <- readRDS(file.path(run, "R_vars/04_coclustering_results.rds"))

solution_cocluster <- read_csv(file.path(run, "04_SNF_subsample_coclusterfrac.csv")) # Coclustering Fraction dataframe
solution_aris <- read_csv(file.path(run, "04_SNF_subsample_pairwisearis.csv")) # Pairwise ARIs dataframe


################################################
########## IDENTIFY THE TOP SOLUTION: #########
# Identify the top solution from coclustering fraction
cocluster_top_row <- solution_cocluster[which.max(solution_cocluster$avg_cocluster_frac), ] # pull the row with the highest coclustering fraction)

cocluster_top_row_number <- cocluster_top_row$row
cocluster_top_metacluster_letter <- solution_cocluster$metacluster[match(cocluster_top_row_number, solution_cocluster$row)]

cocluster_matched_row <- rep_data[rep_data$metacluster == cocluster_top_metacluster_letter, ]
cocluster_top_solution_number <- cocluster_matched_row[, "solution"]
cocluster_top_solution_clusters <- cocluster_matched_row[, "nclust"]

print(paste0("Subject Stability top solution is Metacluster ", cocluster_top_row_number, cocluster_top_metacluster_letter, " Solution # ", cocluster_top_solution_number, " with ", cocluster_top_solution_clusters, " identified clusters"))

# Identify the top solution from pairwise ARIs
aris_top_row <- solution_aris[which.max(solution_aris$mean_ari), ]

aris_top_row_number <- aris_top_row$row
aris_top_metacluster_letter <- solution_aris$metacluster[match(aris_top_row_number, solution_aris$row)]

aris_matched_row <- rep_data[rep_data$metacluster == aris_top_metacluster_letter, ]
aris_top_solution_number <- aris_matched_row[, "solution"]
aris_top_solution_clusters <- aris_matched_row[, "nclust"]

print(paste0("Solution Stability top solution is Metacluster ", aris_top_row_number, aris_top_metacluster_letter, " Solution # ", aris_top_solution_number, " with ", aris_top_solution_clusters, " identified clusters"))

if (cocluster_top_solution_number == aris_top_solution_number) {
  print("Coclustering fraction and pairwise ARIs agree.")
  
  top_row_number <- cocluster_top_row_number
  top_metacluster_letter <- cocluster_top_metacluster_letter
  top_solution_number <- cocluster_top_solution_number
  print(paste0("Top solution is Metacluster ", top_row_number, top_metacluster_letter, " Solution # ", top_solution_number))
} else {
  print("Discrepancies between coclustering fraction and pairwise ARIs. Please inspect manually")
}

################################################
########## IDENTIFY THE TOP SOLUTION: ##########
# Note: because of the way the representative solutions are vectorized chronologically, the metaclusters letter may not map directly.
# Make sure you check this with rep_data!

# USER INPUT REQUIRED: set the top solution printed above
# (row in rep_data, metacluster letter, solution number).
top_row_number <- as.numeric("xxx")
top_metacluster_letter <- "xxx"
top_solution_number <- as.numeric("xxx")

toprep_subassignments <- annotation_data[annotation_data$solution == top_solution_number, ] # <- this is hard coded by the user

# Mean Raw ARIs across 100x100 pairwise matrix for TOP REPRESENTATIVE SOLUTION
cumulative_means <- cummean(pairwise_aris$raw_aris[[paste0("s", top_row_number)]])
toprep_pairwisearis <- ggplot(data.frame(x = seq_along(cumulative_means), y = cumulative_means),
                              aes(x = x, y = y)) +
  geom_point(color = "steelblue", size = 0.5) +
  labs(
    x = "Pairwise Comparison #",
    y = "Cumulative Pairwise ARI",
    title = paste0("Cumulative Pairwise ARI for Solution ", top_metacluster_letter, " #", top_solution_number)
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title  = element_text(size = 12),
    axis.text   = element_text(size = 12)
  )

########################
### DENSITY PLOTS ######
# Density plot: Assess co-clustering fraction
solutions <- sol_rep$solution
cocluster_dfs <- coclustering_results$cocluster_dfs

toprep_coclusterhm <- cocluster_heatmap(cocluster_dfs[[top_row_number]],  
                                        dl = input_dl)

# Be careful here. Representative solutions use metaclusters A-L with corresponding #s from original 500 solutions. 
# Subsampling takes the vector index and assigns 1-2
sol_df_order <- rep_data %>% arrange(metacluster) %>% 
  mutate(row = rank(solution)) %>% 
  rename(idx = row,
         label = metacluster) %>%
  select(-avg_silhouette) %>%
  select(idx, label, solution, nclust)

plot_list <- list()

for (i in 1:nrow(sol_df_order)) {
  
  idx   <- sol_df_order$idx[i] # the row index to pull cocluster_dfs
  label <- sol_df_order$label[i] # metacluster alphabet
  sol <- sol_df_order$solution[i] # the actual solution number
  
  # Find the original row corresponding to the actual solution number to pull cocluster_dfs
  # which() finds the positions in solutions vector where the condition is true
  original_row <- which(solutions == sol)
  
  plot_list[[i]] <- cocluster_density(cocluster_dfs[[original_row]]) +
    ggtitle(paste0("Solution ", label, " (#", sol, ")")) + 
    coord_cartesian(ylim = c(0, 7))
}

density_all <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4, nrow = 3)

toprep_density <- cocluster_density(cocluster_dfs[[top_row_number]]) + ggtitle(paste0("Solution ", top_metacluster_letter, " #", top_solution_number)) + #  <- hard coded by the user
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


#############################
# SAVE NECESSARY DATAFRAMES #

# Inter-subsample ARI heatmap for the top solution
png(file.path(run, "figures/05_SNF_toprep_intersampleARI.png"), width = 12, height = 12, units = "in", res = 300)
ComplexHeatmap::Heatmap(
  pairwise_aris$raw_aris[[paste0("s", top_row_number)]],
  heatmap_legend_param = list(
    color_bar = "continuous",
    title = paste0("Inter-Subsample ARI: Solution ", top_metacluster_letter),
    at = c(0, 0.5, 1)
  ),
  show_column_names = FALSE,
  show_row_names = FALSE
)
dev.off()

# Participant cluster assignments for top solution: in 06_analysis_features.R

# Pairwise ARI over subsample comparisons for top solution
ggsave(file.path(run, "figures/05_SNF_toprep_cumulativeari.png"), toprep_pairwisearis, width = 8, height = 6, units = "in", dpi = 300)

# Coclustering fraction for top solution
png(file.path(run, "figures/05_SNF_toprep_coclusterfrac.png"), width = 12, height = 12, units = "in", res = 300)
toprep_coclusterhm
dev.off()

# Density plot for top solution
ggsave(file.path(run, "figures/05_SNF_toprep_density.png"), toprep_density, width = 8, height = 6, units = "in", dpi = 300)

# Density plot for all representative solutions
ggsave(file.path(run, "figures/04_SNF_repsol_density.png"), density_all, width = 16, height = 12, units = "in", dpi = 300)
