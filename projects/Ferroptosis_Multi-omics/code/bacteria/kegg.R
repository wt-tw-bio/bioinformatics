feature_info <- function(base_path) {
  library(data.table)
  library(dplyr)
  
  # Read genome differences
  genome_diff <- fread(file.path(base_path, "genome_differences.txt"), fill = TRUE, skip = 4)
  colnames(genome_diff) <- c("locus", "type", "start", "end", "len", "len_ref", "diff")
  
  # Find GFF file
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
      return(NULL)
    }))
  }
  
  gene_info <- extract_gene_info(gff_data)
  write.table(gene_info, file = file.path(base_path, "feature_info.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}


feature_info_operon <- function(){
  # Find file
  base_path <- "./k2044/CKP28"
  operon_file <- list.files(base_path, pattern = "list.*", full.names = TRUE)
  # 加载必要的包
  install.packages("zoo")
  library(zoo)
  
  # 读取数据
  df <- data.table::fread(operon_file, header = TRUE, fill = TRUE)
  # 填充空白的 'Operon' 列
  df$Operon <- zoo::na.locf(df$Operon)
  df <- subset(df, IdGene != "")
  colnames(df)[2] <- "ID"
  operon_gff <- list.files(base_path, pattern = "ORF.*", full.names = TRUE)
  operon_dat <- data.table::fread(operon_gff, fill = TRUE, header = FALSE) %>% as.data.frame()
  operon_dat$ID <- stringr::str_extract(operon_dat$V9, "PBMFBDJK_\\d+")
  colnames(operon_dat)[c(1, 4, 5)] <- c("locus", "start", "end")
  gff <- list.files(base_path, pattern = "feature_info.*", full.names = TRUE)
  gff_dat <- data.table::fread(gff, fill = TRUE, header = FALSE)
  colnames(gff_dat)[c(1, 2, 3)] <- c("locus", "start", "end")
  merge_gff <- merge(gff_dat, operon_dat[c(1,4,5,10)], by = c("locus", "start", "end"))
  merge_gff <- merge(merge_gff, df[, 1:2], by = "ID")
  merge_gff <- merge_gff[!duplicated(merge_gff$ID), ]
  colnames(merge_gff)[5:6] <- c("gene", "symbol")
  kegg <- list.files(base_path, pattern = "kegg", full.names = TRUE)
  kegg_id <- read.table(kegg, header = FALSE, fill = TRUE)
  colnames(kegg_id) <- c("symbol", "keggid")
  merge_gff <- merge(merge_gff, kegg_id, by = "symbol")
}

results <- merge_gff %>%
  filter(keggid != "") %>%
  group_by(Operon) %>%
  summarize(keggid_list = list(unique(keggid))) %>%
  pull(keggid_list, Operon)

library(clusterProfiler)
genome_kegg_ids <- unlist(results)
enrich_result <- enrichKEGG(
  gene = results[3][[1]],
  universe = genome_kegg_ids,
  organism = 'ko',
  keyType = 'kegg',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

debugonce(lapply)
res_kegg <- lapply(results, function(x){
    enrich_result <- enrichKEGG(
      gene = x,
      universe = genome_kegg_ids,
      organism = 'ko',
      keyType = 'kegg',
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    return(enrich_result)
})

kegg_dat <- lapply(names(res_kegg), function(x){
  tryCatch(expr = {
    res <- res_kegg[[x]]@result
    dat <- data.frame(Operon = x, category = res$category, subcategory = res$subcategory, ID = res$ID,
                      Description = res$Description, pvalue = res$pvalue, p.adjust = res$p.adjust, 
                      qvalue = res$qvalue)
    return(dat)
  },
  error = function(e) {
    # 处理错误
    message("An error occurred: ", e$message)
    return(NULL)  # 返回一个默认值
  })
}) %>% do.call(rbind, .)

kegg_dat_filter <- subset(kegg_dat, !duplicated(Operon))

sigs <- list.files("./CKP28/", pattern = "^diff", full.names = TRUE)

sigs_dat <- read.table(sigs[1], header = FALSE, sep = "\t")
colnames(sigs_dat)[4] <- "symbol"
merge_gff_sigs <- merge(sigs_dat, merge_gff[, c("symbol", "Operon")], by = "symbol")
merge_gff_sigs <- merge(merge_gff_sigs, kegg_dat_filter, by = "Operon")
merge_gff_sigs <- merge_gff_sigs[!duplicated(merge_gff_sigs), ]
feature_info("./k2044/CKP28")


kegg_id <- read.table("CKP28/CKP28_kegg.txt", header = FALSE, fill = TRUE)
colnames(kegg_id) <- c("symbol", "keggid")
diff_gene <- read.table("CKP28/feature_diff_gene.txt", header = FALSE, fill = TRUE)
kegg_id <- kegg_id[kegg_id$symbol %in% diff_gene$V5, ]
kegg_id_use <- kegg_id[kegg_id$keggid != "", ]
kegg_id_use <- unique(kegg_id_use$keggid)
enrich_result <- enrichKEGG(
  gene = kegg_id_use,
  organism = 'ko',
  keyType = 'kegg',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
res_78578 <- enrich_result@result


kegg_id <- read.table("CKP28/CKP28_kegg.txt", header = FALSE, fill = TRUE)
colnames(kegg_id) <- c("symbol", "keggid")
diff_gene <- read.table("CKP28/feature_diff_gene_k2044.txt", header = FALSE, fill = TRUE)
kegg_id <- kegg_id[kegg_id$symbol %in% diff_gene$V5, ]
kegg_id_use <- kegg_id[kegg_id$keggid != "", ]
kegg_id_use <- unique(kegg_id_use$keggid)
enrich_result <- enrichKEGG(
  gene = kegg_id_use,
  organism = 'ko',
  keyType = 'kegg',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
res_k2044 <- enrich_result@result





library(clusterProfiler)

performEnrichment <- function(directory) {
  # Helper function to perform enrichment analysis
  enrichKEGGAnalysis <- function(diff_gene_file) {
    # Read data
    kegg_id <- read.table(file.path(directory, paste0(basename(directory), "_kegg.txt")), header = FALSE, fill = TRUE)
    colnames(kegg_id) <- c("symbol", "keggid")
    diff_gene <- read.table(file.path(directory, diff_gene_file), header = FALSE, fill = TRUE)
    
    # Filter and process data
    kegg_id <- kegg_id[kegg_id$symbol %in% diff_gene$V5, ]
    kegg_id_use <- kegg_id[kegg_id$keggid != "", ]
    kegg_id_use <- unique(kegg_id_use$keggid)
    
    # Perform enrichment analysis
    enrich_result <- enrichKEGG(
      gene = kegg_id_use,
      organism = 'ko',
      keyType = 'kegg',
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )
    
    return(enrich_result@result)
  }
  
  # Perform enrichment for both datasets
  res_78578 <- enrichKEGGAnalysis("feature_diff_gene.txt")
  res_k2044 <- enrichKEGGAnalysis("feature_diff_gene_k2044.txt")
  
  # Return results as a list
  return(list(
    "78578" = res_78578,
    "k2044" = res_k2044
  ))
}

# Example usage
result <- performEnrichment("CKP28")
folders <- dir("./", pattern = "CKP|HVKP", include.dirs = TRUE)
folders <- folders[file.info(folders)$isdir]
res_kegg_diff <- lapply(folders, performEnrichment)
names(res_kegg_diff) <- folders
library(KEGGREST)

# 获取所有 KEGG 通路
kegg_pathways <- keggList("pathway", "ko")

# 查看前几个通路
head(kegg_pathways)

kegg_pathways[grepl("ferroptosis", kegg_pathways, ignore.case = TRUE)]

lapply(names(res_kegg_diff), function(x) {
  if("KO04216" %in% res_kegg_diff[[x]][["78578"]][["ID"]]) {
    print(paste0(x, " has KO04216 in 78578"))
  }
  if("KO04216" %in% res_kegg_diff[[x]][["k2044"]][["ID"]]) {
    print(paste0(x, " has KO04216 in k2044"))
  }
})


Sys.setenv(http_proxy = "http://127.0.0.1:10809")
Sys.setenv(https_proxy = "http://127.0.0.1:10809")


# 原始字符串
string <- "Amino sugar and nucleotide sugar metabolism，Citrate cycle (TCA cycle)，Biosynthesis of nucleotide sugars，Alanine, aspartate and glutamate metabolism，Pyruvate metabolism"

# 使用 strsplit 分隔字符串
elements <- strsplit(string, "，")[[1]]

meta_keggs <- names(kegg_pathways[kegg_pathways %in% elements])


# 定义一个函数来处理每个元素
process_element <- function(x) {
  check_ids <- function(ids, label) {
    ids_s <- ids[ids %in% meta_keggs]
    if (length(ids_s) > 0) {
      message <- paste0(x, " has ", paste0(ids_s, collapse = "--"), " in ", label)
    } else {
      message <- paste0(x, " has no meta kegg in ", label)
    }
    print(message)
    return(ids_s)
  }
  
  res_78578 <- check_ids(res_kegg_diff[[x]][["78578"]][["ID"]], "78578")
  res_k2044 <- check_ids(res_kegg_diff[[x]][["k2044"]][["ID"]], "k2044")
  return(list(res_78578 = res_78578, res_k2044 = res_k2044))
}

# 使用 lapply 处理每个元素
a <- lapply(names(res_kegg_diff), process_element)
names(a) <- names(res_kegg_diff)
node_meta <- readRDS("node_meta.rds")
idx_78578 <- match(a[["CKP28"]][["res_78578"]], res_kegg_diff[["CKP28"]][["78578"]][["ID"]])
node_78578 <- data.frame(pathway = res_kegg_diff[["CKP28"]][["78578"]][["Description"]][idx_78578],
                         keggid = res_kegg_diff[["CKP28"]][["78578"]][["geneID"]][idx_78578])

node_78578 <- node_78578 %>%
  separate_rows(keggid, sep = "/") %>%
  as.data.frame()

# Read data
kegg_id <- read.table(file.path("CKP28", paste0("CKP28", "_kegg.txt")), header = FALSE, fill = TRUE)
colnames(kegg_id) <- c("symbol", "keggid")
diff_gene <- read.table(file.path("CKP28", "feature_diff_gene.txt"), header = FALSE, fill = TRUE)
colnames(diff_gene)[5] <- "symbol"
diff_gene <- merge(diff_gene, kegg_id, by = "symbol")
colnames(diff_gene)[5] <- "gene"
node_78578 <- merge(node_78578, diff_gene[, c("gene", "keggid")], by = "keggid")

dir.create("cytoscape")
process_kegg_data <- function(prefix) {
  # 构建文件路径
  kegg_file <- file.path(prefix, paste0(prefix, "_kegg.txt"))
  diff_gene_file <- file.path(prefix, "feature_diff_gene.txt")
  
  # 读取数据
  kegg_id <- read.table(kegg_file, header = FALSE, fill = TRUE)
  colnames(kegg_id) <- c("symbol", "keggid")
  diff_gene <- read.table(diff_gene_file, header = FALSE, fill = TRUE)
  colnames(diff_gene)[5] <- "symbol"
  
  # 合并数据
  diff_gene <- merge(diff_gene, kegg_id, by = "symbol")
  colnames(diff_gene)[5] <- "gene"
  
  # 匹配和分割
  idx_78578 <- match(a[[prefix]][["res_78578"]], res_kegg_diff[[prefix]][["78578"]][["ID"]])
  node_78578 <- data.frame(
    pathway = res_kegg_diff[[prefix]][["78578"]][["Description"]][idx_78578],
    keggid = res_kegg_diff[[prefix]][["78578"]][["geneID"]][idx_78578]
  )
  
  node_78578 <- node_78578 %>%
    separate_rows(keggid, sep = "/") %>%
    as.data.frame()
  
  # 最后合并结果
  node_78578 <- merge(node_78578, diff_gene[, c("gene", "keggid")], by = "keggid")
  
  node_78578$group <- "78578"
  
  diff_gene_file_k2044 <- file.path(prefix, "feature_diff_gene_k2044.txt")
  
  diff_gene_k2044 <- read.table(diff_gene_file_k2044, header = FALSE, fill = TRUE)
  colnames(diff_gene_k2044)[5] <- "symbol"
  
  # 合并数据
  diff_gene_k2044 <- merge(diff_gene_k2044, kegg_id, by = "symbol")
  colnames(diff_gene_k2044)[5] <- "gene"
  
  idx_k2044 <- match(a[[prefix]][["res_k2044"]], res_kegg_diff[[prefix]][["k2044"]][["ID"]])
  node_k2044 <- data.frame(
    pathway = res_kegg_diff[[prefix]][["k2044"]][["Description"]][idx_k2044],
    keggid = res_kegg_diff[[prefix]][["k2044"]][["geneID"]][idx_k2044]
  )
  
  node_k2044 <- node_k2044 %>%
    separate_rows(keggid, sep = "/") %>%
    as.data.frame()
  
  node_k2044 <- merge(node_k2044, diff_gene_k2044[, c("gene", "keggid")], by = "keggid")
  
  node_k2044$group <- "k2044"
  node_k2044$gene <- gsub("_\\d+$", "", node_k2044$gene)
  node_78578$gene <- gsub("_\\d+$", "", node_78578$gene)
  dir.create(file.path("cytoscape", prefix), showWarnings = FALSE)
  write.table(node_k2044[, c("pathway", "gene")], file = file.path("cytoscape", prefix, "k2044.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(node_78578[, c("pathway", "gene")], file = file.path("cytoscape", prefix, "78578.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(rbind(node_78578, node_k2044))
}

# 使用函数
result <- process_kegg_data("CKP28")

res <- lapply(folders, process_kegg_data)
names(res) <- folders


write.table(node_meta[, -1], file = file.path("cytoscape", "node_meta.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


