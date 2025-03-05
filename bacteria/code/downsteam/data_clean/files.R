# 设置源文件夹路径
source_folder <- "ref"

# 获取所有子文件夹的名称
subfolders <- list.dirs(source_folder, recursive = FALSE, full.names = FALSE)

# 遍历每个子文件夹
for (subfolder in subfolders) {
  # 构建源文件路径
  source_file <- file.path(source_folder, subfolder, "feature_diff_gene.txt")
  
  # 创建目标子文件夹
  dir.create(subfolder, showWarnings = FALSE)
  
  # 构建目标文件路径
  target_file <- file.path(subfolder, "feature_diff_gene_ref.txt")
  
  # 复制并重命名文件
  if (file.exists(source_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  } else {
    warning(paste("File not found:", source_file))
  }
}

subfolders

dirs <- head(subfolders, -1)

ls_vfdb_ref <- lapply(dirs, function(x) {
  f <- data.table::fread(file.path(x, "vfdb_results.txt")) %>% as.data.frame()
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
  rownames(df_f) <- df_f$id
  diff_gene <- data.table::fread(file.path(x, "feature_diff_gene_ref.txt")) %>% as.data.frame()
  diff_gene$V4 <- gsub("_\\d+$", "", diff_gene$V4)
  rownames(diff_gene) <- diff_gene$V5
  inter <- intersect(diff_gene$V5, df_f$id)
  df_f <- cbind(df_f[inter,], diff_gene[inter, "V4"])
  colnames(df_f)[5] <- ("gene")
  write.table(df_f, file = file.path(x, "diff_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  return(df_f)
})

names(ls_vfdb_ref) <- dirs


library(ggplot2)

lapply(dirs, function(x){
  
  categories <- ls_vfdb_ref[[x]][["category_vfdb"]]
  
  # 创建数据框，计算每个类别的数量
  category_counts <- as.data.frame(table(categories))
  
  # 绘制柱形图
  ggplot(category_counts, aes(x = categories, y = Freq, fill = categories)) +
    geom_bar(stat = "identity") +
    labs(x = "Category", y = "Count", title = "ref") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set3")  # 使用预设调色板
  
  ggsave(file = file.path(x, "category_counts_ref.pdf"), width = 8, height = 6)
})

x <- function(i){
  dg <- data.table::fread(file.path(i, "diff_genes.txt"))
  cat("The differentially expressed genes in ", i, " is:\n ", paste0(unique(dg[[5]]), collapse = "，"), "\n")
}


