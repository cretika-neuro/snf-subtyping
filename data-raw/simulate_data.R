####### Simulate a small multimodal dataset for the SNF pipeline #######
# Generates synthetic data with the same file names, folders, and column
# names the numbered scripts expect. No real participant data are used.
#
# Usage (from the repository root):
#   Rscript data-raw/simulate_data.R
# Output: ./simulated_project/  (use this folder as <project root> in scripts/)

set.seed(2026)

run        <- "snf9"
out_root   <- "simulated_project"
n_asd      <- 60
n_td       <- 70
n_clusters <- 2

# ---- folders ----------------------------------------------------------------
dirs <- c("snf_inputs", "static",
          file.path("snf", run, "R_vars"),
          file.path("snf", run, "figures"),
          file.path("snf", run, "analysis_outputs"))
for (d in dirs) dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)

# ---- feature names (match the names used in scripts 06 and 07) --------------
rsfmri_features <- c("CEN_rdlPFC*SN_ACC", "CEN_lPPC*SN_ACC", "CEN_rPPC*SN_ACC",
                     "SN_ACC*DMN_PCC", "SN_AI*DMN_mPFC", "CEN_rdlPFC*DMN_PCC",
                     "CEN_lPPC*DMN_mPFC", "SN_AI*CEN_lPPC")

genetics_features <- c("PGS_MDD", "PGS_ADHD", "PGS_SCZ", "PGS_EA", "PGS_ASD")

clinical_features <- c(
  "BRIEF_Inhibit", "BRIEF_Shift", "BRIEF_EmoControl", "BRIEF_Initiate",
  "BRIEF_WorkingMemory", "BRIEF_PlanOrganize", "BRIEF_Monitor", "BRIEF_OrganizationMaterials",
  "CBCL_AnxDepress", "CBCL_Withdrawn", "CBCL_Somatic", "CBCL_Rulebreaking",
  "CBCL_Agg", "CBCL_SocialProbs", "CBCL_ThoughtProbs", "CBCL_AttentionProbs",
  "SRS_Awareness", "SRS_Cog", "SRS_Comm", "SRS_Motivation", "SRS_RRB",
  "RBSR_Comp", "RBSR_Ritual", "RBSR_Same", "RBSR_Restrict", "RBSR_Stereo"
)

# ---- participants -----------------------------------------------------------
asd_ids <- sprintf("ASDsim%03d", seq_len(n_asd))
td_ids  <- sprintf("TDsim%03d",  seq_len(n_td))

# latent subtype of each ASD participant (the structure SNF should recover)
latent <- sample(seq_len(n_clusters), n_asd, replace = TRUE, prob = c(0.6, 0.4))

# Standard-normal noise, plus a mean shift for latent cluster 2 on a subset of
# features. Controls (TD) get no shift.
sim_block <- function(feature_names, n, latent = NULL, shift_idx = integer(0), shift = 1) {
  m <- matrix(rnorm(n * length(feature_names)), nrow = n,
              dimnames = list(NULL, feature_names))
  if (!is.null(latent) && length(shift_idx) > 0) {
    m[latent == 2, shift_idx] <- m[latent == 2, shift_idx] + shift
  }
  as.data.frame(m, check.names = FALSE)
}

asd_rsfmri   <- sim_block(rsfmri_features,   n_asd, latent, shift_idx = 1:3, shift = 1.0)
asd_genetics <- sim_block(genetics_features, n_asd, latent, shift_idx = 1,   shift = 0.8)
asd_clinical <- sim_block(clinical_features, n_asd, latent, shift_idx = c(1, 2, 9, 17, 22), shift = 0.8)

td_rsfmri    <- sim_block(rsfmri_features,   n_td)
td_genetics  <- sim_block(genetics_features, n_td)
td_clinical  <- sim_block(clinical_features, n_td)

# ---- SNF input files (ASD participants only, one file per modality) ---------
write_input <- function(df, ids, file) {
  write.csv(cbind(ACE_ID = ids, df), file.path(out_root, "snf_inputs", file),
            row.names = FALSE)
}
write_input(asd_rsfmri,   asd_ids, "input_rsfmri_networks.csv")
write_input(asd_genetics, asd_ids, "input_genetics.csv")
write_input(asd_clinical, asd_ids, "input_clinical.csv")

# ---- master + demographics (ASD and TD) -------------------------------------
n_all <- n_asd + n_td
all_ids <- c(asd_ids, td_ids)
group   <- rep(c("ASD", "TD"), c(n_asd, n_td))

age_years <- round(runif(n_all, 8, 17), 1)
site      <- sample(c("SiteA", "SiteB"), n_all, replace = TRUE)
scanner   <- sample(c("Scan1", "Scan2"), n_all, replace = TRUE)

demo <- data.frame(
  ACE_ID       = all_ids,
  KP_ID        = sprintf("KPsim%03d", seq_len(n_all)),
  NDAR_ID      = sprintf("NDARsim%05d", seq_len(n_all)),
  Site         = site,
  Scanner      = scanner,
  Site_Scanner = paste(site, scanner, sep = "_"),
  Group        = group,
  FSIQ         = round(c(rnorm(n_asd, 105, 15), rnorm(n_td, 115, 12))),
  Age_Months   = round(age_years * 12),
  Age_Years    = age_years,
  Sex          = sample(c("M", "F"), n_all, replace = TRUE, prob = c(0.65, 0.35)),
  stringsAsFactors = FALSE
)

features <- rbind(cbind(asd_rsfmri, asd_genetics, asd_clinical),
                  cbind(td_rsfmri,  td_genetics,  td_clinical))

master <- cbind(demo, features)
write.csv(master, file.path(out_root, "snf_inputs", "SNF_MASTER_ALL.csv"),
          row.names = FALSE)

write.csv(demo[, c("ACE_ID", "Group", "Sex", "Site_Scanner", "Age_Years", "FSIQ")],
          file.path(out_root, "static", "SNF_demos.csv"), row.names = FALSE)

cat("Simulated", n_asd, "ASD and", n_td, "TD participants in", out_root, "\n")
