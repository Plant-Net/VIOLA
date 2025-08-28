#!/usr/bin/env Rscript

## ===============================
## Required Libraries
## ===============================
library(dplyr)
library(data.table)
library(ontologyIndex)
library(tidyr)
library(stringr)
library(optparse)

## ===============================
## Functions
## ===============================

filter_quality_biotype_vaf <- function(df_dbscan) {
  # Filter variants: keep only protein_coding, quality > 50, and PASS in VAF filter
  df_dbscan_f <- df_dbscan %>%
    filter(BIOTYPE == "protein_coding",
           QUAL > 50,
           VAF_filter == "PASS")
  
  # Exclude homozygous reference genotypes if GT column exists
  if ("GT" %in% names(df_dbscan)) {
    df_dbscan_f <- df_dbscan_f %>%
      filter(GT != "0/0")
  }
  return(df_dbscan_f)
}

jaccard <- function(a, b) {
  # Compute Jaccard index between two vectors
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}

jaccard_for_hpo <- function(hpo_patient, hpo_gene) {
  # Compute Jaccard index for all combinations of patient HPO vs gene HPO
  combinaison <- CJ(hpo_patient, hpo_gene, unique = TRUE)
  jaccard_index <- vector("numeric", nrow(combinaison))
  
  for (i in seq_len(nrow(combinaison))) {
    v <- get_ancestors(hpo, combinaison$hpo_patient[i])
    u <- get_ancestors(hpo, combinaison$hpo_gene[i])
    
    if (length(v) > 0 & length(u) > 0) {
      jaccard_index[i] <- jaccard(v, u)
    } else {
      jaccard_index[i] <- NA
    }
  }
  
  combinaison$jaccard_index <- jaccard_index
  return(combinaison)
}

select_anomaly_variant_with_hpo <- function(dbscan_res, hpo_patient, jaccard_th) {
  # Select anomaly variants based on HPO similarity (Jaccard index)
  jaccard_th <- as.numeric(jaccard_th)
  
  anomaly_unique_gene <- unique(dbscan_res$GeneName)
  anomaly_unique_gene <- anomaly_unique_gene[anomaly_unique_gene != ""]
  
  hpo_simili_score_2_gene <- data.frame(GeneName = anomaly_unique_gene,
                                        jaccard_index = 0)
  
  for (i in seq_along(anomaly_unique_gene)) {
    gene <- anomaly_unique_gene[i]
    if (gene %in% hpo_term_md_related_to_a_gene$gene_symbol) {
      hpo_gene <- hpo_term_md_related_to_a_gene$hpo_id[
        hpo_term_md_related_to_a_gene$gene_symbol == gene
      ]
      simili_gene <- jaccard_for_hpo(hpo_patient, hpo_gene)
      hpo_simili_score_2_gene$jaccard_index[i] <- max(simili_gene$jaccard_index, na.rm = TRUE)
    }
  }
  
  filter_jaccard <- hpo_simili_score_2_gene %>%
    filter(jaccard_index >= jaccard_th)
  
  if (nrow(filter_jaccard) > 0) {
    return(dbscan_res %>% inner_join(filter_jaccard, by = "GeneName"))
  } else {
    return(paste0("No gene has a Jaccard index > ", jaccard_th))
  }
}

compute_vscore <- function(df, df_unique_variant) {
  # Compute VIOLA Score
  max_mahal <- max(df$mahalnobis)
  
  df_vs <- df %>%
    left_join(df_module_wgcna, by = "GeneName") %>%
    mutate(in_enriched_module = ifelse(module_name %in% enriched_module, 1, -1),
           unique_ID_2 = paste(X.Chrom, Pos, sep = "_"),
           is_unique = ifelse(unique_ID_2 %in% df_unique_variant$unique_ID, 1, 0),
           known_md_gene = ifelse(GeneName %in% known_md_genes, 1, 0),
           norm_mahal_score = mahalnobis / max_mahal,
           vs = 0.5 * norm_mahal_score + 0.5 * in_enriched_module +
             0.5 * is_unique + 0.01 * known_md_gene,
           vs_mm = scale(vs, center = min(vs), scale = max(vs) - min(vs))
    ) %>%
    arrange(desc(vs_mm)) %>%
    mutate(rank_vs = rank(-vs_mm, ties.method = "max"))
  
  threshold <- quantile(df_vs$vs_mm, 0.75, na.rm = TRUE)
  df_vs <- df_vs %>% filter(vs_mm > threshold)
  return(df_vs)
}

compute_ar_score <- function(df) {
  # Compute AR Score (autosomal recessive model)
  df_heterozygote <- df %>%
    filter(GT != "1/1") %>%
    group_by(GeneName) %>%
    summarise(n_var = n(), .groups = "drop") %>%
    filter(n_var > 1)
  
  df_ar <- df %>%
    filter(GT == "1/1" | GeneName %in% df_heterozygote$GeneName) %>%
    arrange(desc(vs)) %>%
    mutate(rank_ar = rank(-vs, ties.method = "max"))
  
  return(df_ar)
}

get_unique_variant <- function(vcf_file) {
  # Load unique rare variants from VCF
  df_unique_var <- fread(vcf_file, sep = '\t', header = TRUE, skip = '#CHROM')
  names(df_unique_var)[1] <- "CHROM"
  df_unique_var$CHROM <- gsub("chr", "", df_unique_var$CHROM)
  df_unique_var$unique_ID <- paste(df_unique_var$CHROM, df_unique_var$POS, sep = "_")
  return(df_unique_var)
}

## ===============================
## Main
## ===============================

option_list <- list(
  make_option(c("-f", "--input"), dest = "input_file", help = "Path to input file"),
  make_option(c("-o", "--output_path"), dest = "output_path", help = "Path to output directory"),
  make_option(c("-t", "--hpo_table"), dest = "hpo_table_path", help = "Path to patient HPO table"),
  make_option(c("-p", "--path_unique_var"), dest = "path_unique_var", help = "Path to VCF of unique rare variants"),
  make_option(c("-r", "--resources"), dest = "resources_path", help = "Path to resources folder containing all required files")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

input_file      <- args$input_file
output_path     <- args$output_path
hpo_table_path  <- args$hpo_table_path
path_unique_var <- args$path_unique_var
resources_path  <- args$resources_path

file_prefix <- strsplit(basename(input_file), "_")[[1]][1]

## Build paths to resources
hpo_gene_map_file <- file.path(resources_path, "hpoterms_md_to_gene.txt")
wgcna_modules_file <- file.path(resources_path, "df_all_merged_module_wgcna.csv")
known_md_genes_file <- file.path(resources_path, "Mitochondrial_disease_genes_Prokish.csv")

## Step 1: Quality filter
message(date(), " : Apply filters on quality and biotype")
dbscan_res <- read.csv(input_file)
dbscan_res_f <- filter_quality_biotype_vaf(dbscan_res)
write.csv(dbscan_res_f,
          file = file.path(output_path, paste0(file_prefix, "_res_filter_quality_biotype_vaf.csv")),
          row.names = FALSE)

## Step 2: HPO filter
message(date(), " : Apply HPO filter")
hpo_term_md_related_to_a_gene <- read.delim(hpo_gene_map_file)
hpo <- get_ontology("https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo",
                    extract_tags = "everything")
table_hpo <- read.csv(hpo_table_path)
hpo_patient <- table_hpo$hpo_id[table_hpo$patient_id == file_prefix]

res_HPO <- select_anomaly_variant_with_hpo(dbscan_res_f, hpo_patient, 0.8)
write.csv(res_HPO,
          file = file.path(output_path, paste0(file_prefix, "_res_HPO.csv")),
          row.names = FALSE)

## Step 3: Mahalanobis score
message(date(), " : Compute Mahalanobis score")
features <- c('SIFTval','PolyPhenVal','priPhCons','mamPhCons','verPhCons','priPhyloP','mamPhyloP','verPhyloP',
              'EncodeH3K27ac.max','EncodeH3K27ac.sum','EncodeH3K4me1.max','EncodeH3K4me1.sum','EncodeH3K4me3.sum',
              'EncodeH3K4me3.max','EncodeDNase.sum','EncodeDNase.max','MMSp_acceptorIntron','MMSp_acceptor','MMSp_exon',
              'MMSp_donor','MMSp_donorIntron','Freq100bp','Rare100bp','Sngl100bp','Freq1000bp','Rare1000bp',
              'Sngl1000bp','Freq10000bp','Rare10000bp','Sngl10000bp','RawScore','PHRED')

data <- res_HPO[, features] %>%
  mutate(across(all_of(features), as.double)) %>%
  mutate(across(all_of(features), ~ replace_na(., median(., na.rm = TRUE))))

data$mahalnobis <- mahalanobis(data, colMeans(data), cov(data), tol = 1e-25)
data$pvalue <- pchisq(data$mahalnobis, df = 3, lower.tail = FALSE)
subset_res <- cbind(res_HPO, data[, c("mahalnobis", "pvalue")]) %>%
  filter(pvalue < 0.001)

## Step 4: VIOLA Score
message(date(), " : Compute VIOLA Score and AR ranking")
df_module_wgcna <- read.csv(wgcna_modules_file)
enriched_module <- c("blue","cyan","brown","midnightblue","black","green","salmon","greenyellow")
df_known_md_genes <- read.csv(known_md_genes_file)
known_md_genes <- df_known_md_genes$Gene.name

vcf_file <- file.path(path_unique_var, file_prefix, paste0(file_prefix, "_unique_rare_variant.vcf"))
df_unique_var <- get_unique_variant(vcf_file)

res_vs <- compute_vscore(subset_res, df_unique_var)
write.csv(res_vs, file = file.path(output_path, paste0(file_prefix, "_vrank.csv")), row.names = FALSE)

## Step 5: AR rank
res_vs_genotype <- compute_ar_score(res_vs)
write.csv(res_vs_genotype, file = file.path(output_path, paste0(file_prefix, "_ar_rank.csv")), row.names = FALSE)
