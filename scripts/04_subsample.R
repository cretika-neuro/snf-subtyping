####### METASNF: Step 4 ######

# NOTE: this script can run automatically, ~15 minutes
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metasnf)
library(SNFtool)
library(ComplexHeatmap)

rm(list = ls())

# Set working directory here (left blank on purpose): <project root>/snf/

run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)
select <- dplyr::select

# Load R objects
input_dl <- readRDS(file.path(run, "R_vars/01_input_dl.rds"))

my_sconfig <- readRDS(file.path(run, "R_vars/02_my_sconfig.rds"))     # SNF config file for all solutions
sol_df <- readRDS(file.path(run, "R_vars", "02_SNF_soldf.rds"))       # Participant cluster assignment for all solutions
annotation_data <- readRDS(file.path(run, "R_vars/03_annotation_data.rds"))

rep_solutions <- readRDS(file.path(run, "R_vars/03_rep_solutions.rds")) # Participant cluster assignment per representative solution
rep_data <- readRDS(file.path(run, "R_vars/03_rep_sil.rds"))          # Additional information for each representative solution

################################
##### STABILITY MEASURES #######

# Calculate resampling-based stability measures
my_dl_subsamples <- subsample_dl(input_dl,
                                 n_subsamples = 100,        # number of subsamples per representative solution
                                 subsample_fraction = 0.80) # fraction of original sample

rep_solution <- rep_solutions$solution
num_repsol <- length(rep_solution)
row_idx <- which(my_sconfig$settings_df$solution %in% rep_solution)

# Here, we create a new config file for only the representative solutions to account for the weights matrix. 
# filtered_config is the same as config_rep

filtered_config <- my_sconfig # make a copy
filtered_config$settings_df <- my_sconfig$settings_df[row_idx, ]
filtered_config$weights_matrix <- my_sconfig$weights_matrix[row_idx, ]

##################################################
############### SUBSAMPLING ######################
# batch_snf_subsamples: run SNF clustering pipeline on a list of subsampled data lists
# returns a one-element list (cluster_solutions) which is a list of cluster solution data frames corresponding to each data list subsamples
batch_subsample_results <- batch_snf_subsamples(my_dl_subsamples,
                                                filtered_config,
                                                processes = "max",
                                                return_sim_mats = FALSE)

##################################################
################# SOLUTION STABILITY #############
# Calculate ARIs between cluster solutions across subsamples
# pairwise_aris(): list containing ARIs between subsamples for each row of dataframe
pairwise_aris <- subsample_pairwise_aris(batch_subsample_results)
solution_aris <- data.frame(pairwise_aris$ari_summary) 
solution_aris <- data.frame(pairwise_aris$ari_summary %>% rename(row = solution))

# Be careful here. pairwise_aris gives the CHRONOLOGICAL order of the representative solutions, not the letters.
# Assign solution and metacluster mapping with rep_data
rep_data$row <- rank(rep_data$solution)
rep_data <- rep_data %>% select(row, solution, everything())

solution_aris <- solution_aris %>% left_join(rep_data, by="row")

subsample_pairwisearis <- ggplot(solution_aris, aes(x = metacluster, y = mean_ari)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(
    title = paste0("Solution Stability per Representative Solution"),
    x = "Solution",
    y = "Mean Pairwise ARI"
  ) +
  geom_hline(yintercept = mean(solution_aris$mean_ari),
             color = "red", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title  = element_text(size = 14, color = "black"),
    axis.text   = element_text(size = 14, color = "black")
  )

##################################################
################# SUBJECT STABILITY #############
# Calculate how pair of observations cluster together across subsamples
sol_rep <- sol_df[rep_solution, ] # Participant cluster assignments per representative solution

coclustering_results <- calculate_coclustering(batch_subsample_results, sol_rep)

# coclustering_results(): dataframe per cluster solution showing 
# cocluster_dfs: 
# cocluster_ss_mats: 
# cocluster_sc_mats: number of times every pair of obs occurred in same cluster
# cocluster_cf_mats: fraction of times every pair of obs occurred in same cluster
# cocluster_summary: among pairs of obs that cluster in original solution, mean fraction of pairs that remain clustered in subsample iterations

# Coclustering Average per Representative Solution
# Be careful here, this is also CHRONOLOGICAL order: updating rank order with rep_data
solution_cocluster <- coclustering_results$cocluster_summary %>% left_join(rep_data, by = "row")

cocluster_fraction <- ggplot(solution_cocluster, aes(x = metacluster, y = avg_cocluster_frac)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_hline(yintercept = mean(solution_cocluster$avg_cocluster_frac), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Subject Stability per Representative Solution",
    x = "Solution",
    y = "Mean Coclustering Fraction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title  = element_text(size = 14, color = "black"),
    axis.text   = element_text(size = 14, color = "black")
  )

####################################
##### SAVE NECESSARY DATAFRAMES ####
# Subsample configs
saveRDS(my_dl_subsamples, file.path(run, "R_vars/04_my_dl_subsamples.rds"))

# Subsample results
saveRDS(batch_subsample_results, file.path(run, "R_vars/04_batch_subsample_results.rds"))

# Mean Pairwise ARIs across subsamples per representative solution
write_csv(solution_aris, file.path(run, "04_SNF_subsample_pairwisearis.csv"))
saveRDS(pairwise_aris, file.path(run, "R_vars/04_pairwise_aris.rds"))
ggsave(file.path(run, "figures/04_SNF_subsample_pairwisearis.png"), subsample_pairwisearis, width = 8, height = 6, units = "in", dpi = 300)

# Coclustering fraction per representative solution
saveRDS(sol_rep, file.path(run, "R_vars/04_sol_rep.rds"))
saveRDS(coclustering_results, file.path(run, "R_vars/04_coclustering_results.rds"))
write_csv(solution_cocluster, file.path(run, "04_SNF_subsample_coclusterfrac.csv"))
ggsave(file.path(run, "figures/04_SNF_subsample_coclusterfrac.png"), cocluster_fraction, width = 8, height = 6, units = "in", dpi = 300)


# SAVE REPRESENTATIVE SOLUTIONS TOGETHER
#coclustering results: coclustering_results$cocluster_summary
cocluster <- coclustering_results$cocluster_summary

# pairwise aris: pairwise_aris$ari_summary$mean_ari
pairwise <- pairwise_aris$ari_summary %>% rename(row="solution") %>% 
  dplyr::select(mean_ari, ari_sd, row)

repsol_info <- rep_data %>% left_join(cocluster, by="row") %>% left_join(pairwise, by = "row")
repsol_info <- repsol_info %>% dplyr::select(row, metacluster, solution, nclust, avg_cocluster_frac, mean_ari, ari_sd, avg_silhouette)

write_csv(repsol_info, file.path(run, "04_SNF_repsol_info.csv"))
