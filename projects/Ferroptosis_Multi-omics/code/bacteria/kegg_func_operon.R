library(data.table)
library(dplyr)
library(zoo)
library(stringr)
library(clusterProfiler)

# Function to extract feature information
extract_feature_info <- function(base_path, diff_path) {
  genome_diff <- fread(file.path(base_path, "genome_differences.txt"), fill = TRUE, skip = 4)
  colnames(genome_diff) <- c("locus", "type", "start", "end", "len", "len_ref", "diff")
  
  gff_file <- list.files(base_path, pattern = "\\.gff$", full.names = TRUE)
  if (length(gff_file) == 0) stop("No GFF file found.")
  gff_lines <- readLines(gff_file[1])
  header_indices <- grep("^##", gff_lines)
  gene_indices <- tail(header_indices, 2)
  gff_data <- gff_lines[(gene_indices[1] + 1):(gene_indices[2] - 1)]
  
  extract_gene_info <- function(data) {
    do.call(rbind, lapply(data, function(line) {
      fields <- strsplit(line, "\t")[[1]]
      locus <- fields[1]
      start <- as.numeric(fields[4])
      end <- as.numeric(fields[5])
      gene_match <- regmatches(line, regexec("gene=([^;]+)", line))
      gene <- ifelse(length(gene_match[[1]]) > 1, gene_match[[1]][2], NA)
      protein_match <- regmatches(line, regexec("locus_tag=([^;]+)", line))
      protein <- ifelse(length(protein_match[[1]]) > 1, protein_match[[1]][2], NA)
      return(data.frame(locus = locus, start = start, end = end, gene = gene, protein = protein))
    }))
  }
  
  gene_info <- extract_gene_info(gff_data)
  write.table(gene_info, file = file.path(diff_path, "feature_info.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

# Function to process operon data and perform KEGG enrichment
process_operon_and_enrich <- function(diff_path, ref) {
  operon_file <- list.files(diff_path, pattern = "list.*", full.names = TRUE)
  
  df <- fread(operon_file, header = TRUE, fill = TRUE)
  df$Operon <- na.locf(df$Operon)
  df <- subset(df, IdGene != "")
  colnames(df)[2] <- "ID"
  
  operon_gff <- list.files(diff_path, pattern = "ORF.*", full.names = TRUE)
  operon_dat <- fread(operon_gff, fill = TRUE, header = FALSE) %>% as.data.frame()
  operon_dat$ID <- str_extract(operon_dat$V9, "(?<=ID=).*(?=;pro)")
  colnames(operon_dat)[c(1, 4, 5)] <- c("locus", "start", "end")
  
  gff <- list.files(diff_path, pattern = "feature_info.*", full.names = TRUE)
  gff_dat <- fread(gff, fill = TRUE, header = FALSE)
  colnames(gff_dat)[c(1, 2, 3)] <- c("locus", "start", "end")
  
  merge_gff <- merge(gff_dat, operon_dat[c(1, 4, 5, 10)], by = c("locus", "start", "end"))
  merge_gff <- merge(merge_gff, df[, 1:2], by = "ID")
  merge_gff <- merge_gff[!duplicated(merge_gff$ID), ]
  colnames(merge_gff)[5:6] <- c("gene", "symbol")
  
  kegg <- list.files(diff_path, pattern = "kegg", full.names = TRUE)
  kegg_id <- read.table(kegg, header = FALSE, fill = TRUE)
  colnames(kegg_id) <- c("symbol", "keggid")
  merge_gff <- merge(merge_gff, kegg_id, by = "symbol")
  
  results <- merge_gff %>%
    filter(keggid != "") %>%
    group_by(Operon) %>%
    summarize(keggid_list = list(unique(keggid))) %>%
    pull(keggid_list, Operon)
  
  genome_kegg_ids <- unlist(results)
  res_kegg <- lapply(names(results), function(x) {
    tryCatch({
      enrich_result <- enrichKEGG(
        gene = results[[x]],
        universe = genome_kegg_ids,
        organism = 'ko',
        keyType = 'kegg',
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      )
      res <- enrich_result@result
      data.frame(Operon = x, category = res$category, subcategory = res$subcategory, ID = res$ID,
                 Description = res$Description, pvalue = res$pvalue, p.adjust = res$p.adjust, qvalue = res$qvalue)
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NULL)
    })
  }) %>% do.call(rbind, .)
  
  kegg_dat_filter <- subset(res_kegg, !duplicated(Operon))
  
  sigs <- list.files(diff_path, pattern = "^diff", full.names = TRUE)
  sigs_dat <- fread(sigs[ref], header = FALSE, sep = "\t") %>% as.data.frame()
  colnames(sigs_dat)[4] <- "symbol"
  
  merge_gff_sigs <- merge(sigs_dat, merge_gff[, c("symbol", "Operon")], by = "symbol")
  merge_gff_sigs <- merge(merge_gff_sigs, kegg_dat_filter, by = "Operon")
  merge_gff_sigs <- merge_gff_sigs[!duplicated(merge_gff_sigs), ]
  
  return(merge_gff_sigs)
}

# Example usage:
base_path <- "./78578/HVKP1/"
diff_path = "./HVKP1"
extract_feature_info(base_path, diff_path)
#ref = 1是k2044参考基因组
results <- process_operon_and_enrich(diff_path, ref = 1)

#kegg_operon
#k2044
lapply(folders, function(x){
  base_path <- file.path("k2044", x)
  diff_path = x
  extract_feature_info(base_path, diff_path)
  #ref = 1是k2044参考基因组
  results <- process_operon_and_enrich(diff_path, ref = 1)
  write.table(results, file = file.path(x, "operon_k2044.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
})

#78578
lapply(folders, function(x){
  base_path <- file.path("78578", x)
  diff_path = x
  extract_feature_info(base_path, diff_path)
  #ref = 1是k2044参考基因组
  results <- process_operon_and_enrich(diff_path, ref = 2)
  write.table(results, file = file.path(x, "operon_78578.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
})

debugonce("process_operon_and_enrich")
results <- process_operon_and_enrich(diff_path = "HVKP1", ref = 2)

