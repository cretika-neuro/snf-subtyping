####### METASNF: Step 3 ######
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

# USER INPUT REQUIRED: inspect the metacluster heatmap with shiny_annotator()
# and set split_vec below. Everything else runs automatically.

rm(list = ls())

# Set working directory here (left blank on purpose): <project root>/snf/

run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)

select <- dplyr::select

# Load R objects from previous script
my_sconfig <- readRDS(file.path(run, "R_vars", "02_my_sconfig.rds"))
sol_df <- readRDS(file.path(run, "R_vars", "02_SNF_soldf.rds"))

################################
######## METACLUSTERS ##########

# Cluster the clustering solutions with Adjusted Rand Index (ARI): how similar they are to each other)
# This step takes a bit to run!
sol_aris <- calc_aris(sol_df)
head(sol_aris, 3)

# Hierarchical clustering-based ordering of rows in ARIs
meta_cluster_order <- get_matrix_order(sol_aris)

# Create metacluster heatmap
ari_heatmap <- meta_cluster_heatmap(sol_aris, order = meta_cluster_order)

# Inspect interactive heat map
shiny_annotator(ari_heatmap)

# USER INPUT: row indices at which to split the heatmap into metaclusters.
# Populates the metacluster column (mc) in the solutions dataframe.
split_vec <- c(xxx) # left blank on purpose, refer to shiny_annotator hetamap

# Split cluster solution heatmap. Visually inspect.
ari_mc_hm <- meta_cluster_heatmap(sol_aris, 
                                  order = meta_cluster_order, 
                                  split_vector = split_vec)
ari_mc_hm

mc_sol_df <- label_meta_clusters(sol_df, order = meta_cluster_order, split_vector = split_vec)

###################################
######## QUALITY INDICES ##########
##### SILHOUETTE SCORES
# Calculate silhouette scores: similaritity to own cluster than other clusters (range: -1 to 1)
silhouette_scores <- calculate_silhouettes(sol_df) 

# Average the silhouette scores across all participants in each solution. 
mean_silhouettes <- data.frame(solution = 1:nrow(sol_df),
                               mean_silhouette = sapply(silhouette_scores, function(x) mean (x[,"sil_width"])))

annotation_data <- as.data.frame(mc_sol_df) %>%
  left_join(mean_silhouettes, by = "solution") %>%
  dplyr::select(solution, nclust, mc, mean_silhouette, everything())

# Silhouette scores per solution, per metacluster
solution_silhouette_plot <- ggplot(annotation_data, aes(x = mc, y = mean_silhouette)) +
  geom_boxplot(fill = "steelblue", alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "steelblue4") +
  #geom_hline(yintercept = mean(annotation_data$mean_silhouette),
             #color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Silhouette Score Distribution of Each Metacluster",
    x = "Metacluster",
    y = "Silhouette Score (per solution)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title  = element_text(size = 14, color = "black"),
    axis.text   = element_text(size = 14, color = "black")
  )

# Silhouette scores averaged per meta cluster
silhouette_mc_df <- annotation_data %>%
  group_by(mc) %>%
  summarise(
    mean_sil = mean(mean_silhouette),
    min_sil  = min(mean_silhouette),
    max_sil  = max(mean_silhouette)
  ) %>%
  arrange(mc)

n_mc <- nrow(silhouette_mc_df)

metacluster_silhouette <- ggplot(silhouette_mc_df, aes(x = factor(LETTERS[1:n_mc]), y = mean_sil)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_hline(yintercept = mean(silhouette_mc_df$mean_sil), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Mean Silhouette Score per Metacluster",
    x = "Metacluster",
    y = "Mean Silhouette Score"
  ) +
  scale_x_discrete(labels = LETTERS[1:n_mc]) +
  scale_y_continuous(limits = c(0, 0.085)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title  = element_text(size = 14, color = "black"),
    axis.text   = element_text(size = 14, color = "black")
  )

# annotated heatmap
annotated_hm <- meta_cluster_heatmap(
  sol_aris,
  order = meta_cluster_order,
  split_vector = split_vec,
  data = annotation_data,
  bottom_bar = list("Silhouette Score" = "mean_silhouette",
                    "# of Clusters" = "nclust")
)

# View Participant Cluster Assignments of representative solutions
# get_representative_solutions() : extracts one cluster solution per meta cluster based on highest average ARI
rep_solutions <- get_representative_solutions(sol_aris, annotation_data)
print(paste0("Number of representative solutions: ", length(rep_solutions$solution)))
print(paste0(rep_solutions$solution))

########################################
#### SELECT REPRESENTATIVE SOLUTIONS ####
# Silhouette scores of representative data
rep_sil <- data.frame(
  solution = rep_solutions$solution,
  nclust = rep_solutions$nclust,
  metacluster = rep_solutions$mc,
  avg_silhouette = rep_solutions$mean_silhouette
) 

# Solution #'s in alphabetical order of metacluster 
solution_order <- rep_sil %>% arrange(metacluster) %>% pull(solution)
  
repsol_silhouette <- ggplot(rep_sil, aes(x = factor(solution, levels = solution_order), y = avg_silhouette)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_hline(yintercept = mean(rep_sil$avg_silhouette),
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Mean Silhouette Score per Metacluster Representative Solution",
    x = "Metacluster Representative Solution",
    y = "Mean Silhouette Score"
  ) +
  scale_y_continuous(limits = c(0, 0.12)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title  = element_text(size = 14, color = "black"),
    axis.text   = element_text(size = 14, color = "black")
  )

# HYPERPARAMETERS of representative solutions
config_rep <- my_sconfig # make a copy of the original config
config_rep$settings_df <- my_sconfig$settings_df %>% # "trim" down to the representative solutions
  filter(solution %in% rep_solutions$solution)

config_rep_mclabel <- config_rep$settings_df %>%
  left_join(rep_sil %>% select(solution, metacluster), by = "solution") 

# Heatmap of representative solution hyperparameters in metacluster alphabetical order
# The ORDER of row_positions vector calls the row that will be plotted first.
# For example, row_positions <- c(8,9,2,1,3,7,10,6,4,5) will plot row 8, then row 9, etc.
row_positions <- match(solution_order, config_rep$settings_df$solution)

config_hm <- config_heatmap(config_rep$settings_df, 
                            hide_fixed = TRUE,
                            include_weights = TRUE,
                            include_settings = TRUE,
                            order = row_positions)

# Participant Cluster Assignments per representative solutions
sol_rep <- mc_sol_df %>%
  filter(solution %in% rep_solutions$solution)

# Get the number of participants per cluster, per cluster solution
plot_list <- list()

rep_data_unique <- rep_sil[order(rep_sil$metacluster),]

solution_indices <- rep_data_unique$solution
solution_labels  <- rep_data_unique$metacluster

for (i in seq_along(solution_indices)) {
  idx   <- solution_indices[i]
  label <- solution_labels[i]
  
  cluster_assignments <- as.numeric(mc_sol_df[mc_sol_df$solution == idx,
                                              grep("^uid_", colnames(mc_sol_df))])
  
  plot_df <- as.data.frame(table(cluster_assignments))
  colnames(plot_df) <- c("Cluster", "Count")
  
  p <- ggplot(plot_df, aes(x = Cluster, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "white") +
    geom_text(aes(label = Count), vjust = -0.5, size = 3) +
    labs(
      title = paste0("Solution ", label, " (#", idx, ")"),
      x = "Cluster",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title  = element_text(size = 9),
      axis.text   = element_text(size = 8)
    )
  
  plot_list[[i]] <- p
}

# Combine figures
combined_clusternum <- gridExtra::grid.arrange(grobs = plot_list, ncol = 4, nrow = 3)

############################################
######### FINAL SIMILARITY MATRICES ########

# Visualize final fused similarity matrix from selected cluster solution
similarity_matrices <- sim_mats_list(mc_sol_df)

simmats_list <- list()

for (i in seq_along(solution_indices)) {
  idx   <- solution_indices[i]
  label <- solution_labels[i]
  
  similarity_matrix <- similarity_matrices[[idx]]
  cluster_solution  <- t(mc_sol_df)
  
  similarity_matrix_hm <- similarity_matrix_heatmap(
    similarity_matrix = similarity_matrix,
    cluster_solution  = cluster_solution,
    heatmap_height    = grid::unit(10, "cm"),
    heatmap_width     = grid::unit(10, "cm"),
    column_title      = paste0("Final Similarity Matrix: Solution ", label, " (#", idx, ")")
  )
  
  simmats_list[[i]] <- grid.grabExpr(ComplexHeatmap::draw(similarity_matrix_hm))
}

# Combine and export
simmats_combined <- gridExtra::grid.arrange(grobs = simmats_list, ncol = 4, nrow = 3)

####################################
##### SAVE NECESSARY DATAFRAMES ####
# Solutions ordered by pairwise ARIs in heatmap
png(file.path(run, "figures/02_SNF_sol_heatmap.png"), width = 8, height = 8, units = "in", res = 300)
ari_heatmap
dev.off()

# Pairwise ARIs heatmap split by metacluster 
png(file.path(run, "figures/02_SNF_sol_heatmapsplit.png"), width = 8, height = 8, units = "in", res = 300)
ari_mc_hm
dev.off()

# Solution dataframe with metacluster label
saveRDS(mc_sol_df, file.path(run, "R_vars/02_mc_sol_df.rds"))

################ SILHOUETTE SCORE
# Silhouette score distribution of each metacluster
ggsave(file.path(run, "figures/03_SNF_metacluster_soldis_silhouette.png"), solution_silhouette_plot, width = 8, height = 6, units = "in", dpi = 300)

# Silhouette scores averaged per metacluster
write_csv(silhouette_mc_df, file.path(run, "03_SNF_metacluster_silhouette.csv"))
ggsave(file.path(run, "figures/03_SNF_metacluster_meansilhouette.png"), metacluster_silhouette, width = 8, height = 6, units = "in", dpi = 300)

# Heatmap with silhouette score annotations
png(file.path(run, "figures/03_SNF_metacluster_heatmapannotated.png"), width = 8, height = 8, units = "in", res = 300)
annotated_hm
dev.off()

# Solution dataframe with metacluster label and mean silhouette score per solution
saveRDS(annotation_data, file.path(run, "R_vars/03_annotation_data.rds"))
write_csv(annotation_data, file.path(run, "03_annotation_data.csv"))

############### REPRESENTATIVE SOLUTIONS
# Participant cluster assignment per representative solution
write_csv(rep_solutions, file.path(run, "03_SNF_repsol_participantclusters.csv"))
saveRDS(rep_solutions, file.path(run, "R_vars/03_rep_solutions.rds"))

# Mean silhouette score of ONLY each representative solution per metacluster
write_csv(rep_sil, file.path(run, "03_SNF_repsol_meansilhoutte.csv"))
saveRDS(rep_sil, file.path(run, "R_vars/03_rep_sil.rds"))
ggsave(file.path(run, "figures/03_SNF_repsol_meansilhoutte.png"), repsol_silhouette, width = 8, height = 6, units = "in", dpi = 300)

# Hyperparameters per representative solution
write_csv(config_rep_mclabel, file.path(run, "03_SNF_repsol_settings.csv"))

# Heatmap of hyperparameters per representative solution in metacluster order
png(file.path(run, "figures/03_SNF_repsol_settings_mcorder.png"), width = 6, height = 8, units = "in", res = 300)
config_hm
dev.off()

# Number of participants per cluster per representative solution in metacluster order
ggsave(file.path(run, "figures/03_SNF_repsol_nclust.png"), combined_clusternum, width = 12, height = 12, units = "in", dpi = 300)

# All similarity matrices combined
png(file.path(run, "figures/03_SNF_repsol_simmats.png"), width = 20, height = 15, units = "in", res = 500)
gridExtra::grid.arrange(grobs = simmats_list, ncol = 4, nrow = 3)
dev.off()
