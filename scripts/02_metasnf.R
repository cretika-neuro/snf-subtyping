####### METASNF: Step 2 ######
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metasnf)
library(SNFtool)

rm(list = ls())

# Load R objects from previous script
# Set working directory here (left blank on purpose): <project root>/snf/

run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)

input_dl <- readRDS(file.path(run, "R_vars/01_input_dl.rds"))

################################
##### GENERATE SNF CONFIG ######
# The SNF Config stores all necessary information to convert raw data into final cluster solution
# Settings df: tracks one set of hyperparameters per row (alpha, k, t, scheme, clustering and distance function, component dropout)
# Weights matrix: tracks one set of feature weights per row
# Distance functions list: stores all functions used to convert to distance matrices
# Clustering algorithms list: stores all functions used to convert networks to cluster solutions

set.seed(42)
my_sconfig <- snf_config(
  dl = input_dl,
  n_solutions = 500,
  min_k = 10,
  max_k = 40,
  min_alpha = 0.3,
  max_alpha = 0.8
)
my_sconfig

####################################
############# RUN SNF ##############

# This step will take a while to load
# batchsnf(): integrates datalist using all iterated hyperparameters and functions in SNF config
# produces: solution index, cluster solution (nclust), metacluster (mc), and clusters each subject was assigned to
sol_df <- batch_snf(input_dl, my_sconfig, return_sim_mats = TRUE)

# Transposed format
sub_cluster_solutions <- t(sol_df)

# Plot solution distributions of cluster #. 
snf_nclust <- ggplot(sol_df, aes(x = factor(nclust))) +
  geom_bar(fill = "#4d87b0", color = "black", linewidth = 0.3) +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.5, size = 4.5) +
  labs(
    title = paste0(nrow(sol_df), " Solutions: Distribution of Clusters"),
    x = "nclust",
    y = "Count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "plain"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line  = element_line(color = "black"),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

####################################
##### SAVE NECESSARY DATAFRAMES ####

# SNF Config File
write_csv(my_sconfig$settings_df, file.path(run, "02_SNF_solsettings.csv"))
saveRDS(my_sconfig, file.path(run, "R_vars", "02_my_sconfig.rds"))

# SNF Solutions
saveRDS(sol_df, file.path(run, "R_vars", "02_SNF_soldf.rds"))

# Number of clusters per solution
ggsave(file.path(run, "figures/02_SNF_sol_nclust.png"), snf_nclust, width = 8, height = 6, units = "in", dpi = 300)
