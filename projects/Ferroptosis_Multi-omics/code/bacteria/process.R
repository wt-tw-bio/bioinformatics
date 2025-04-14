# 加载必要的库
library(dplyr)
library(stringr)
# 读取 GFF 文件
gff_file <- "../outdir/Prokka/CKP28/CKP28.gff"
gff_data <- read.table(gff_file, sep="\t", header=FALSE, comment.char="#")

gff_file <- list.files("../outdir/Prokka/CKP28/", pattern = "\\.gff$", full.names = TRUE)
gff_lines <- readLines(gff_file[1])
header_indices <- grep("^##", gff_lines)
gene_indices <- tail(header_indices, 2)
gff_data <- gff_lines[(gene_indices[1] + 1):(gene_indices[2] - 1)]

gff_df <- lapply(gff_data, function(x) {
  split_data <- strsplit(x, "\t")[[1]]
  description <- split_data[9]
  res <- data.frame(locus = split_data[1], type = split_data[3], start = split_data[4], end = split_data[5], strand = split_data[7], 
                    id = str_extract(description, "(?<=ID=).+?(?=;|$)"),
                    gene = str_extract(description, "(?<=gene=).+?(?=;|$)"), 
                    cog = str_extract(description, "(?<=COG:).+?(?=;|$)"),
                    description = str_extract(description, "(?<=product=).+?(?=;|$)"))
  return(res)
}) %>% do.call(rbind, .) %>% as.data.frame()

gbk_file <- list.files("../outdir/Prokka/CKP28/", pattern = "\\.gbk$", full.names = TRUE)
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
merged_gene_info <- merge(gff_df, cumulative_lengths, by = "locus")
merged_gene_info$start <- as.numeric(merged_gene_info$start) + merged_gene_info$length
merged_gene_info$end <- as.numeric(merged_gene_info$end) + merged_gene_info$length


f <- data.table::fread("./vfdb_results.txt")
f <- subset(f, grepl("Klebsiella pneumoniae ", V10))
f <- f[!duplicated(f$V1), ]

# 正则表达式模式
pattern <- "\\(([^)]+)\\) \\(([^)]+)\\)\\s*([^\\[]+) \\[([^\\(]+)\\(([^\\)]+)\\) - ([^\\(]+)\\(([^\\)]+)\\)\\] \\[([^\\]]+)\\]"

# 使用str_match提取信息
matches <- str_match(f$V10[2], pattern)

df_f <- lapply(seq_len(nrow(f)), function(x){
  matches <- str_match(f$V10[x], pattern)
  data.frame(id_vfdb = matches[5], category_vfdb = matches[7], des_vfdb = matches[4])
}) %>% do.call(rbind, .) %>% as.data.frame()

df_f$id <- f$V1

merged_gene_info <- merge(merged_gene_info, df_f, by = "id", all.x = TRUE)


merged_gene_info <- merge(merged_gene_info, m_cog_select[, c("id", "cogid")], all.x = TRUE)

merged_gene_info <- merge(merged_gene_info, cog_map[, c("cogid", "Gene", "Func")], by = "cogid", all.x = TRUE)

capitalize_first_lower <- function(x) {
  paste0(tolower(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}


merged_gene_info$Gene <- capitalize_first_lower(merged_gene_info$Gene)
merged_gene_info <- merged_gene_info[!duplicated(merged_gene_info$id), ]
colnames(merged_gene_info)[13] <- c("category")
merged_gene_info <- merge(merged_gene_info, cog_category_des, by = "category", all.x = TRUE)







#给我一个具有25个颜色的色盘
circos_colors <- c(
  "red", "green", "blue", "yellow", "purple", "cyan",
  "orange", "brown", "pink", "grey", "violet", "magenta", "lime",
  "olive", "navy", "teal", "maroon", "aqua", "fuchsia", "silver",
  "gold", "lavender", "coral", "indigo","black"
)

merged_gene_info <- merge(merged_gene_info, data.frame(color = circos_colors, category = cog_category_des$category), by = "category", all.x = TRUE)



forward <- data.table::fread("./forward.txt")
forward <- as.data.frame(forward)
reverse <- data.table::fread("./reverse.txt")
colnames(forward)[4] <- "id"
forward <- merge(forward, merged_gene_info[, c("id", "color.y")], by = "id", all.x = TRUE)
forward$color.y[is.na(forward$color.y)] <- "white"
colnames(reverse)[4] <- "id"
reverse <- merge(reverse, merged_gene_info[, c("id", "color")], by = "id", all.x = TRUE)
reverse$color[is.na(reverse$color)] <- "#ffffff"
write.table(forward[,-1], file = "./forward.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(reverse[,-1], file = "./reverse.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

k2044 <- data.table::fread("./k2044/CKP28/feature_diff.txt")
k2044 <- cbind(v1 = "chr1", k2044[, 1:2], V2 = "red")
m78578 <- data.table::fread("./78578/CKP28/feature_diff.txt")
m78578 <- cbind(v1 = "chr1", m78578[, 1:2], V2 = "blue")
write.table(k2044, file = "./feature_diff_k2044.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(m78578, file = "./feature_diff_m78578.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
#白色和黑色的颜色编码
"#ffffff"
"#000000"

library(clusterProfiler)

# 加载 KEGGREST 包
library(KEGGREST)

# 定义物种缩写
species <- "kpn"

# 获取所有基因信息
genes <- keggList("genes", "kpn")

# 提取基因 ID
gene_ids <- names(genes)

# 打印前几个基因 ID
print(head(gene_ids))


gene <- gsub("_\\d+", "", a$gene) %>% unique()
enrichKEGG(
  gene,
  organism = "kpn",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)



library(httr)
library(jsonlite)

# KEGG REST API 基本 URL
base_url <- "https://rest.kegg.jp"

# 假设你知道物种缩写
species <- "kpn"

# 获取基因列表
response <- GET(paste0(base_url, "/list/genes/", species))
genes <- content(response, "text")

# 解析基因列表
gene_lines <- strsplit(genes, "\n")[[1]]
gene_ids <- sapply(gene_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) > 1) {
    return(parts[1])
  }
  return(NA)
})

# 过滤掉 NA 的值
gene_ids <- gene_ids[!is.na(gene_ids)]

# 打印前几个基因 ID
print(head(gene_ids))

