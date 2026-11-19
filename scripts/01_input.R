####### METASNF Step 1: Assemble SNF inputs ######
# One cleaned datasheet per modality (same participants, no missing values),
# combined into a metasnf data_list.
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(metasnf)

rm(list = ls())

select <- dplyr::select

########################################################
########### DATASHEET OF ALL MODALITIES, CLEANED ######
########################################################
# Datasheets must have no missing values and have the same participants
run <- "xxx" # set snf run here (exm: snf9, left blank on purpose)
# Set working directory here (left blank on purpose): <project root>/snf_inputs/


# Read in original datasheets
snf_rsfmri_data <- read.csv("input_rsfmri_networks.csv") %>%
  mutate(ACE_ID = str_remove_all(as.character(ACE_ID), "\\s+"), ACE_ID = str_remove(ACE_ID, "\\.00$")) # standardize ID strings

snf_genetics_data <- read.csv("input_genetics.csv")
snf_behavior_data <- read.csv("input_clinical.csv") 

############################
# Put everything into a list and check dimensions
snf_dflist <- list(snf_rsfmri_data, snf_genetics_data, snf_behavior_data)
lapply(snf_dflist, dim)

# Check if any values are missing. If complete, skip imputation
lapply(snf_dflist,
       function(x) {
         any(is.na(x))
       })

# List of subjects with complete data
snf_completeid <- get_complete_uids(snf_dflist, uid = "ACE_ID")
print(length(snf_completeid))

###############################
###### GENERATE DATA_LIST #####
# Structured list of data frames (complete clustering inputs) with metadatas. Check 
# to see if other data will do a good job separating observations
# (data, name, domain, type)

input_dl <- data_list(
  list(snf_rsfmri_data, "rsfmri", "neuroimaging", "continuous"),
  list(snf_genetics_data, "genetics", "genetics", "continuous"),
  list(snf_behavior_data, "behavior", "behavior", "continuous"),
  uid = "ACE_ID"
)

summary(input_dl)

####################################
##### SAVE NECESSARY DATAFRAMES ####
# Set working directory here (left blank on purpose): <project root>/snf/

# Datalist with all input dataframes for SNF
saveRDS(input_dl, file.path(run, "R_vars", "01_input_dl.rds"))
