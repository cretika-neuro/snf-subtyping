####### METASNF Step 6: Cluster and group feature analysis ######
library(readr)
library(readxl)
library(ggplot2)
library(rstatix)
library(dunn.test)
library(ggpubr)
library(plotly)
library(effectsize)
library(MASS)
library(tidyverse)
library(emmeans)
library(lmerTest)
library(lme4)
library(broom.mixed)
library(MuMIn)
library(caret)
library(car)

rm(list = ls())
select <- dplyr::select
options(scipen=999)

# Set working directory here (left blank on purpose): <project root>/

run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)

##############################
#### IDENTIFIED TOP SOLUTION #
##############################
# USER INPUT REQUIRED: top solution identified in 05_toprep.R

top_row_number <- as.numeric("xxx")
top_metacluster_letter <- "xxx"
top_solution_number <- as.numeric("xxx")

top_cluster <- paste0("s",top_solution_number)

rep_solutions <- readRDS(file.path("snf", run, "R_vars/03_rep_solutions.rds")) # Participant cluster assignment per representative solution

top_rep_assignment <- rep_solutions %>%
  select(solution, starts_with("uid_")) %>%            # drop nclust, mc, silhouette
  pivot_longer(-solution, names_to = "ACE_ID", values_to = "cluster") %>%
  mutate(ACE_ID = sub("^uid_", "", ACE_ID)) %>%               # uid_ASD070 -> ASD070
  pivot_wider(names_from = solution, values_from = cluster,
              names_prefix = "sol_")    
top_rep_assignment <- top_rep_assignment %>% select(ACE_ID, paste0("sol_", top_solution_number)) %>%
  rename(clusters=paste0("sol_", top_solution_number))
write_csv(top_rep_assignment, file.path("snf", run, "05_SNF_toprep_subassignments.csv"))

##############################
#### DATA FOR ANALYSIS #######
##############################
# Participant cluster labels
clusterlabels <- read_csv(file.path("snf", run, "05_SNF_toprep_subassignments.csv"))
table(clusterlabels$clusters)

analysis_df <- read_csv("snf_inputs/SNF_MASTER_ALL.csv") %>% 
  mutate(ACE_ID = str_remove_all(as.character(ACE_ID), "\\s+"), ACE_ID = str_remove(ACE_ID, "\\.00$")) %>% # standardize ID strings
  left_join(clusterlabels, by="ACE_ID") %>% relocate(clusters, .after="ACE_ID") %>%
  mutate(clusters=replace_na(as.character(clusters), "4")) 

# Participant IDs to exclude (set locally; kept out of the public repo)
excluded_ids <- character(0)

demo_data <- read_csv("static/SNF_demos.csv") %>% filter(!ACE_ID %in% excluded_ids) %>%
  mutate(ACE_ID = str_remove_all(as.character(ACE_ID), "\\s+"), ACE_ID = str_remove(ACE_ID, "\\.00$")) %>% # standardize ID strings
  left_join(clusterlabels, by="ACE_ID") %>% relocate(clusters, .after="ACE_ID") %>%
  mutate(clusters=replace_na(as.character(clusters), "4"))
  
table(analysis_df$Group)
table(analysis_df$clusters)

######################################################
######## SAMPLE DESCRIPTIVES: GROUP ANALYSIS #########
# Test demographics on the full sample
# ASD vs. TD
# Sex per Group
table(demo_data$Sex, demo_data$Group)
chisq.test(table(demo_data$Sex, demo_data$Group))$expected 

sex_group_test <- if (all(chisq.test(table(demo_data$Sex, demo_data$Group))$expected >= 5)) {
  chisq.test(table(demo_data$Sex, demo_data$Group))
} else {
  fisher.test(table(demo_data$Sex, demo_data$Group), simulate.p.value = TRUE, B = 10000)
}
print(sex_group_test)

# Scanner per Group
table(demo_data$Site_Scanner, demo_data$Group) 
chisq.test(table(demo_data$Site_Scanner, demo_data$Group))$expected  

site_group_test <- if (all(chisq.test(table(demo_data$Site_Scanner, demo_data$Group))$expected >= 5)) {
  chisq.test(table(demo_data$Site_Scanner, demo_data$Group))
} else {
  fisher.test(table(demo_data$Site_Scanner, demo_data$Group), simulate.p.value = TRUE, B = 10000)
}
print(site_group_test)

# Age per Group
demo_data$Age_Years <- as.numeric(demo_data$Age_Years)
tapply(demo_data$Age_Years, demo_data$Group, summary)
tapply(demo_data$Age_Years, demo_data$Group, sd)

age_test_group <- if (shapiro.test(demo_data$Age_Years)$p.value > 0.05) {
  t.test(Age_Years ~ Group, data = demo_data)
} else {
  wilcox.test(Age_Years ~ Group, data = demo_data)
}
print(age_test_group)

# FSIQ per Group
demo_data$FSIQ <- as.numeric(demo_data$FSIQ)
tapply(demo_data$FSIQ, demo_data$Group, summary)
tapply(demo_data$FSIQ, demo_data$Group, sd)

fsiq_test_group <- if (shapiro.test(demo_data$FSIQ)$p.value > 0.05) {
  t.test(FSIQ ~ Group, data = demo_data)
} else {
  wilcox.test(FSIQ ~ Group, data = demo_data)
}
print(fsiq_test_group)

###########################################################
################## CLUSTER DESCRIPTIVES #######################
# Test demographics across clusters (cluster 4 = TD)
table(demo_data$clusters) 
demo_12 <- demo_data %>% filter(clusters %in% c(1,2))
demo_1td <- demo_data %>% filter(clusters %in% c(1,4))
demo_2td <- demo_data %>% filter(clusters %in% c(2,4))

ggplot(demo_data, aes(x=factor(clusters))) + 
  geom_bar(fill = "steelblue") +
  geom_text(stat="count", aes(label = after_stat(count)), vjust = -0.5) + 
  labs(title = paste0("Solution ", top_metacluster_letter, " #", top_solution_number),
       x = NULL,
       y = "Count") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_blank())

# Sex per cluster
table(demo_data$Sex, demo_data$clusters)
chisq.test(table(demo_12$Sex, demo_12$clusters))$expected 
chisq.test(table(demo_1td$Sex, demo_1td$clusters))$expected 
chisq.test(table(demo_2td$Sex, demo_2td$clusters))$expected 

chisq.test(table(demo_data$Sex, demo_data$clusters))

demos <- list(demo_12 = demo_12, demo_1td = demo_1td, demo_2td = demo_2td)

for (demo in names(demos)) {
  cat("Sex Chi-square test for: ", demo)
  df <- demos[[demo]]
  
  sex_test <- chisq.test(table(df$Sex, df$clusters))
  print(sex_test)
}

# Scanner per cluster
scannertab <- table(demo_data$Site_Scanner, demo_data$clusters) 
chisq.test(scannertab)$expected
chisq.test(scannertab)
chisq.test(scannertab, simulate.p.value = TRUE, B = 1e5)
fisher.test(scannertab, simulate.p.value = TRUE, B = 1e5)

for (demo in names(demos)) {
  cat("Scanner test for : ", demo)
  
  df <- demos[[demo]]
  
  tbl <- table(df$Site_Scanner, df$clusters)
  
  chisq_res <- chisq.test(tbl)
  
  sitescanner_test <- if (all(chisq_res$expected >= 5)) {
    chisq_res
  } else {
    fisher.test(tbl, simulate.p.value = TRUE, B = 10000)
  }
  print(sitescanner_test)
}

# Age per cluster
demo_data$Age_Years <- as.numeric(demo_data$Age_Years)
tapply(demo_data$Age_Years, demo_data$clusters, summary)
tapply(demo_data$Age_Years, demo_data$clusters, sd)

fit <- aov(Age_Years ~ clusters, data = demo_data)
summary(fit)

for (demo in names(demos)) {
  cat("Age test for : ", demo)
  
  df <- demos[[demo]]
  df$clusters <- factor(df$clusters)
  
  # Check normality within each cluster group
  p_vals <- tapply(df$Age_Years, df$clusters, function(x) shapiro.test(x)$p.value)
  
  age_test <- if (all(p_vals > 0.05, na.rm = TRUE)) {
    t.test(Age_Years ~ factor(clusters), data = df)
  } else {
    wilcox.test(Age_Years ~ factor(clusters), data = df)
  }
  print(age_test)  
}

# FSIQ per cluster (ANOVA, then pairwise tests)
demo_data$FSIQ <- as.numeric(demo_data$FSIQ)
tapply(demo_data$FSIQ, demo_data$clusters, summary)
tapply(demo_data$FSIQ, demo_data$clusters, sd)

fit <- aov(FSIQ ~ clusters, data = demo_data)
summary(fit)

for (demo in names(demos)) {
  cat("FSIQ test for : ", demo)
  
  df <- demos[[demo]]
  df$clusters <- factor(df$clusters)
  
  # Check normality within each cluster group
  p_vals <- tapply(df$FSIQ, df$clusters, function(x) shapiro.test(x)$p.value)
  
  fsiq_test <- if (all(p_vals > 0.05, na.rm = TRUE)) {
    t.test(FSIQ ~ factor(clusters), data = df)
  } else {
    wilcox.test(FSIQ ~ factor(clusters), data = df)
  }
  print(fsiq_test)  
}

#########################################################
########### CLUSTER ANALYSIS: PAIRWISE LINEAR MODELS ###############
#########################################################

# Pairwise linear models (feature ~ cluster + FSIQ): Cluster 1 vs. Cluster 2
features_clusters <- analysis_df %>% select(-ACE_ID, -Group, -clusters, -KP_ID, -NDAR_ID, -Site, -Scanner, -Site_Scanner, -FSIQ, -Age_Months, -Age_Years, -Sex) %>% names()

# remove reward-network rsfMRI features
features_clusters <- features_clusters %>% .[!str_detect(., "Amygdala|Thalamus|NucleusAccumbens")] 

# remove diffusion
features_clusters <- features_clusters[!str_starts(features_clusters, "FA_")]
length(features_clusters)

# Format covariates (FSIQ is z-scored in place; run this block once per session)
analysis_df_final <- analysis_df %>% 
  mutate(
    clusters = factor(clusters),
    FSIQ = as.numeric(scale(FSIQ)),
    Site_Scanner = factor(Site_Scanner)
  ) %>%
  relocate(FSIQ, .after = clusters)

data_c1_c2 <- analysis_df_final %>%
  filter(clusters %in% c("1", "2")) %>%
  mutate(clusters = factor(clusters, levels = c("2", "1")))   # reference = Cluster 2; estimate = Cluster 1 - Cluster 2

data_c1_td <- analysis_df_final %>%
  filter(clusters %in% c("1", "4")) %>%
  mutate(clusters = factor(clusters, levels = c("4", "1")))   # reference = TD; estimate = Cluster 1 - TD

data_c2_td <- analysis_df_final %>%
  filter(clusters %in% c("2", "4")) %>%
  mutate(clusters = factor(clusters, levels = c("4", "2")))   # reference = TD; estimate = Cluster 2 - TD

results_c1_c2 <- lapply(features_clusters, function(col) {
  formula <- as.formula(paste0("`", col, "` ~ clusters + FSIQ"))
  model <- lm(formula, data = data_c1_c2)
  
  # R2 of the full model (cluster + FSIQ)
  rsq <- summary(model)$r.squared
  
  # Tidy coefficients, keep the cluster row
  broom::tidy(model) %>%
    filter(str_starts(term, "clusters")) %>%
    mutate(
      feature = col,
      comparison = "Cluster 1 vs. Cluster 2",
      R_square = rsq
    )
}) %>% bind_rows() %>%
  mutate(fdr_pvalue = p.adjust(p.value, method = "fdr")) %>%
  select(comparison, feature, fdr_pvalue, p.value, everything()) %>%
  arrange(fdr_pvalue)

write_csv(results_c1_c2, file.path("snf", run, "analysis_outputs/analysis_cluster1cluster2lm_features.csv"))

sig_pairwise_cluster_results <- results_c1_c2 %>% filter(fdr_pvalue<0.05)
write_csv(sig_pairwise_cluster_results, file.path("snf", run, "analysis_outputs/analysis_cluster1cluster2lm_features_fdrsig.csv"))

sig_features_clusters <- sig_pairwise_cluster_results$feature
sort(sig_features_clusters)

############### TD ANALYSIS
# Comparison 1: Cluster 1 vs. TD (estimate = Cluster 1 - TD)
results_c1_td <- lapply(sig_features_clusters, function(col) {
  formula <- as.formula(paste0("`", col, "` ~ clusters + FSIQ"))
  model <- lm(formula, data = data_c1_td)
  
  # R2 of the full model (cluster + FSIQ)
  rsq <- summary(model)$r.squared
  
  # Tidy coefficients, keep the cluster row
  broom::tidy(model) %>%
    filter(str_starts(term, "clusters")) %>%
    mutate(
      feature = col,
      comparison = "Cluster 1 vs. TD",
      R_square = rsq
    )
}) %>% bind_rows()

# Comparison 2: Cluster 2 vs. TD (estimate = Cluster 2 - TD)
results_c2_td <- lapply(sig_features_clusters, function(col) {
  formula <- as.formula(paste0("`", col, "` ~ clusters + FSIQ"))
  model <- lm(formula, data = data_c2_td)
  
  # R2 of the full model (cluster + FSIQ)
  rsq <- summary(model)$r.squared
  
  # Tidy coefficients, keep the cluster row
  broom::tidy(model) %>%
    filter(str_starts(term, "clusters")) %>%
    mutate(
      feature = col,
      comparison = "Cluster 2 vs. TD",
      R_square = rsq
      )
}) %>% bind_rows()

# FDR correction across both TD comparisons (Cluster 1 vs. TD, Cluster 2 vs. TD)
pairwise_clustertd_results <- bind_rows(results_c1_td, results_c2_td) %>%
  mutate(fdr_pvalue = p.adjust(p.value, method = "fdr")) %>%
  select(comparison, feature, fdr_pvalue, p.value, everything()) %>%
  arrange(fdr_pvalue)

write_csv(pairwise_clustertd_results, file.path("snf", run, "analysis_outputs/analysis_cluster1tdcluster2td_features.csv"))

sig_pairwise_clustertd_results <- pairwise_clustertd_results %>% filter(fdr_pvalue<0.05)
write_csv(sig_pairwise_clustertd_results, file.path("snf", run, "analysis_outputs/analysis_cluster1tdcluster2td_features_fdrsig.csv"))


################################################
############# GROUP ANALYSIS ##################
###############################################
# ASD vs. TD analysis
# Linear model per significant feature (feature ~ Group + FSIQ)
# Format covariates (FSIQ is z-scored in place; run this block once per session)
analysis_df_group <- analysis_df %>% 
  mutate(
    Group = factor(Group),
    FSIQ = as.numeric(scale(FSIQ))
  ) %>%
  relocate(FSIQ, Site_Scanner, .after = clusters)

analysis_df_group <- analysis_df_group %>%
  mutate(Group = relevel(factor(Group), ref = "TD"))

# Dependent variables
sig_features_forgroup <- sig_features_clusters
results_group <- lapply(sig_features_forgroup, function(col) {
  dd    <- analysis_df_group %>% filter(!is.na(.data[[col]]), !is.na(FSIQ))
  model <- lm(as.formula(paste0("`", col, "` ~ Group + FSIQ")), data = dd)
  
  broom::tidy(model, conf.int = TRUE) %>%
    filter(str_starts(term, "Group")) %>%
    mutate(feature = col, comparison = "ASD vs. TD",
           n_TD = sum(dd$Group == "TD"), n_ASD = sum(dd$Group == "ASD"),
           R_square = summary(model)$r.squared)
}) %>%
  bind_rows() %>%
  mutate(FDR_pvalue = p.adjust(p.value, method = "BH")) %>%
  arrange(FDR_pvalue) %>%
  select(feature, comparison, term, n_TD, n_ASD, estimate, conf.low, conf.high,
         p.value, FDR_pvalue, R_square) %>%
  relocate(FDR_pvalue, p.value, .after="feature")


write_csv(results_group, file.path("snf", run, "analysis_outputs/analysis_groupfeatures.csv"))
