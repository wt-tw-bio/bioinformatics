process_genome_data <- function(base_path) {
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
      if (grepl("\\bgene\\b", line)) {
        fields <- strsplit(line, "\t")[[1]]
        locus <- fields[1]
        start <- as.numeric(fields[4])
        end <- as.numeric(fields[5])
        gene_match <- regmatches(line, regexec("gene=([^;]+)", line))
        gene <- ifelse(length(gene_match[[1]]) > 1, gene_match[[1]][2], NA)
        protein_match <- regmatches(line, regexec("locus_tag=([^;]+)", line))
        protein <- ifelse(length(protein_match[[1]]) > 1, protein_match[[1]][2], NA)
        return(data.frame(locus = locus, start = start, end = end, gene = gene, protein = protein))
      }
      return(NULL)
    }))
  }
  
  gene_info <- extract_gene_info(gff_data)
  
  # Find GBK file
  gbk_file <- list.files(base_path, pattern = "\\.gbk$", full.names = TRUE)
  if (length(gbk_file) == 0) stop("No GBK file found.")
  gbk_lines <- readLines(gbk_file[1])
  
  locus_lines <- gbk_lines[grep("LOCUS", gbk_lines)]
  locus_info <- list()
  for (i in 1:length(locus_lines)) {
    matches <- regmatches(locus_lines[i], regexec("LOCUS\\s+(\\S+)\\s+(\\d+)", locus_lines[i]))
    locus <- matches[[1]][2]
    length <- matches[[1]][3]
    locus_info[[i]] <- data.frame(locus = locus, length = length)
  }
  
  cumulative_lengths <- do.call(rbind, locus_info)
  cumulative_lengths$length <- c(0, head(cumulative_lengths$length, -1)) %>% as.numeric()
  cumulative_lengths$length <- cumsum(cumulative_lengths$length)
  merged_gene_info <- merge(gene_info, cumulative_lengths, by = "locus")
  merged_gene_info$start <- merged_gene_info$start + merged_gene_info$length
  merged_gene_info$end <- merged_gene_info$end + merged_gene_info$length
  
  genome_diff$locus <- as.character(genome_diff$locus)
  merged_diff_info <- merge(genome_diff, cumulative_lengths, by = "locus")
  merged_diff_info$start <- merged_diff_info$start + merged_diff_info$length
  merged_diff_info$end <- merged_diff_info$end + merged_diff_info$length
  merged_diff_info <- within(merged_diff_info, {
    swap_needed <- end < start
    temp <- ifelse(swap_needed, start, end)
    start <- ifelse(swap_needed, end, start)
    end <- temp
  })
  
  # Remove temporary columns
  merged_diff_info$swap_needed <- NULL
  merged_diff_info$temp <- NULL
  merged_diff_info <- as.data.frame(merged_diff_info)
  feature_diff <- merged_diff_info[c("start", "end", "type")]
  colnames(feature_diff) <- c("#start", "end", "label")
  feature_diff <- cbind(chr = "chr1", feature_diff)
  write.table(feature_diff, file = file.path(base_path, "feature_diff.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Merge data
  result <- merged_diff_info %>%
    rowwise() %>%
    mutate(
      type = case_when(
        type == "GAP" & diff > 0 ~ "GAP-insertion",
        type == "GAP" & diff < 0 ~ "GAP-deletion",
        TRUE ~ type
      ),
      gene = paste(
        merged_gene_info$gene[merged_gene_info$start <= end & merged_gene_info$end >= start],
        collapse = ","
      )
    )
  
  # Construct list
  gene_list <- result %>%
    ungroup() %>%
    group_by(type) %>%
    summarise(genes = paste(unique(gene[gene != ""]), collapse = ",")) %>%
    { setNames(as.list(.$genes), .$type) }
  gene_list <- lapply(gene_list, function(x) strsplit(x, ",")[[1]])
  
  feature_diff_gene <- merged_gene_info[merged_gene_info$gene %in% unlist(gene_list), ]
  
  # Read and filter VFDB results
  #vfdb_file <- file.path(base_path, "vfdb_results.txt")
  #vfdb_data <- fread(vfdb_file)
  #siderophore_subset <- subset(vfdb_data, grepl("siderophore", vfdb_data$V10))
  #feature_diff_gene <- feature_diff_gene[feature_diff_gene$protein %in% siderophore_subset$V1, ]
  feature_diff_gene2 <- feature_diff_gene[c("start", "end", "gene", "protein")]
  colnames(feature_diff_gene2)[1:3] <- c("#start", "end", "label")
  feature_diff_gene2 <- cbind(chr = "chr1", feature_diff_gene2)
  write.table(feature_diff_gene2, file = file.path(base_path, "feature_diff_gene.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}


process_genome_data("./k2044/CKP28")
dirs <- list.dirs("./k2044")
dirs <- dirs[grepl("./C|./H", dirs)]
lapply(dirs, process_genome_data)

dirs <- list.dirs("./78578/")
dirs <- dirs[grepl("./C|./H", dirs)]
lapply(dirs, process_genome_data)



process_genome_data("./k2044/CKP28")
process_genome_data("./k2044/")


a <- lapply(basename(list.dirs("./78578/"))[-1], function(x){
  res <- data.table::fread(file.path("./78578", x, "feature_diff_gene.txt")) %>% as.data.frame()
  return(res[, 3])
})


file.copy(list.files("./78578/", pattern = "fna", recursive = TRUE, full.names = TRUE), "./78578/query/")
dir.create("./k2044/query/")
file.copy(list.files("./k2044/", pattern = "fna", recursive = TRUE, full.names = TRUE), "./k2044/query/")
