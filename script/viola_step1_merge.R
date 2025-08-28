#!/usr/bin/Rscript

library(dplyr)
library(tidyr)
library(optparse)

#################
## ARGUMENTS ####
#################

option_list <- list(
  make_option(c("-v", "--vep_input"), dest = "vep_input",
              help = "Path to VEP input file"),
  make_option(c("-c", "--cadd_input"), dest = "cadd_input",
              help = "Path to CADD input file"),
  make_option(c("-l", "--clinvar_input"), dest = "clinvar_input",
              help = "Path to ClinVar input file"),
  make_option(c("-o", "--output_path"), dest = "output_path",
              help = "Path to output directory")
)

args <- parse_args(OptionParser(option_list = option_list))

vep_input <- args$vep_input
cadd_input <- args$cadd_input
clinvar_input <- args$clinvar_input
output_path <- args$output_path

id_patient <- strsplit(basename(vep_input), "_")[[1]][1]

########################################################################
######################## MAIN ##########################################
########################################################################

## Load ClinVar database
message(date(), " : Loading ClinVar database")
clinvar <- read.delim(clinvar_input)

######################################
######## LOAD VEP RESULTS ############
######################################

message(date(), " : Loading VEP results")
df_vep <- read.delim(vep_input)

# Required VEP columns
required_vep_cols <- c("CHROM", "POS", "Gene", "Feature", "AD", "GT", "DP")
missing_vep <- setdiff(required_vep_cols, colnames(df_vep))
if (length(missing_vep) > 0) {
  stop("Missing required columns in VEP input: ", paste(missing_vep, collapse = ", "))
}

df_vep <- df_vep %>%
  mutate(
    CHROM = gsub("chr", "", CHROM),
    uniqueID = paste(CHROM, POS, Gene, Feature, sep = "_")
  ) %>%
  separate(AD, into = c("AD_REF", "AD_ALT1", "AD_ALT2"), sep = ",", convert = TRUE) %>%
  mutate(
    VAF_ALT1 = AD_ALT1 / DP,
    VAF_ALT2 = AD_ALT2 / DP,
    VAF_filter = case_when(
      GT == "1/1" & VAF_ALT1 >= 0.8 ~ "PASS",
      GT == "0/1" & VAF_ALT1 >= 0.3 & VAF_ALT1 <= 0.7 ~ "PASS",
      GT == "1/2" & VAF_ALT1 >= 0.3 & VAF_ALT1 <= 0.7 & VAF_ALT2 >= 0.3 & VAF_ALT2 <= 0.7 ~ "PASS",
      TRUE ~ "FAIL"
    )
  )

######################################
######## LOAD CADD RESULTS ###########
######################################

message(date(), " : Loading CADD results")
df_cadd <- read.delim(cadd_input, skip = 1)

# Required CADD columns
required_cadd_cols <- c("X.Chrom", "Pos", "GeneID", "FeatureID")
missing_cadd <- setdiff(required_cadd_cols, colnames(df_cadd))
if (length(missing_cadd) > 0) {
  stop("Missing required columns in CADD input: ", paste(missing_cadd, collapse = ", "))
}

df_cadd <- df_cadd %>%
  mutate(
    uniqueID = paste(X.Chrom, Pos, GeneID, FeatureID, sep = "_"),
    uniqueID_2 = paste(X.Chrom, Pos, sep = "_")
  )

######################################
######## MERGE DATASETS ##############
######################################

df_merge <- df_vep %>%
  inner_join(df_cadd, by = "uniqueID") %>%
  left_join(clinvar, by = "uniqueID_2")

# Save output
output_file <- file.path(output_path, paste0(id_patient, "_rare_variant_cadd_input_viola.csv"))
write.csv(df_merge, output_file, row.names = FALSE)

message(date(), " : Finished. Output saved to ", output_file)
