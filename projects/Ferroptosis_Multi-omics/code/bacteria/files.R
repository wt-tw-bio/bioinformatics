# 设置源文件夹路径
source_folder <- "k2044"

# 获取所有子文件夹的名称
subfolders <- list.dirs(source_folder, recursive = FALSE, full.names = FALSE)

# 遍历每个子文件夹
for (subfolder in subfolders) {
  # 构建源文件路径
  source_file <- file.path(source_folder, subfolder, "feature_diff_gene.txt")
  
  # 创建目标子文件夹
  dir.create(subfolder, showWarnings = FALSE)
  
  # 构建目标文件路径
  target_file <- file.path(subfolder, "feature_diff_gene_k2044.txt")
  
  # 复制并重命名文件
  if (file.exists(source_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  } else {
    warning(paste("File not found:", source_file))
  }
}


# 设置源文件夹路径
source_folder <- "78578"

# 获取所有子文件夹的名称
subfolders <- list.dirs(source_folder, recursive = FALSE, full.names = FALSE)

# 遍历每个子文件夹
for (subfolder in subfolders) {
  # 构建源文件路径
  source_file <- file.path(source_folder, subfolder, "feature_diff_gene.txt")
  
  
  # 构建目标文件路径
  target_file <- file.path(subfolder, "feature_diff_gene.txt")
  
  # 复制并重命名文件
  if (file.exists(source_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  } else {
    warning(paste("File not found:", source_file))
  }
}


# 设置源文件夹路径
source_folder <- "78578"

# 获取所有子文件夹的名称
subfolders <- list.dirs(source_folder, recursive = FALSE, full.names = FALSE)

# 遍历每个子文件夹
for (subfolder in subfolders) {
  # 构建源文件路径
  source_file <- file.path(source_folder, subfolder, "vfdb_results.txt")
  
  
  # 构建目标文件路径
  target_file <- file.path(subfolder, "vfdb_results.txt")
  
  # 复制并重命名文件
  if (file.exists(source_file)) {
    file.copy(source_file, target_file, overwrite = TRUE)
  } else {
    warning(paste("File not found:", source_file))
  }
}

subfolders

dirs <- head(subfolders, -1)

ls_vfdb_k2044 <- lapply(dirs, function(x) {
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
  diff_gene <- data.table::fread(file.path(x, "feature_diff_gene_k2044.txt")) %>% as.data.frame()
  diff_gene$V4 <- gsub("_\\d+$", "", diff_gene$V4)
  rownames(diff_gene) <- diff_gene$V5
  inter <- intersect(diff_gene$V5, df_f$id)
  df_f <- cbind(df_f[inter,], diff_gene[inter, "V4"])
  colnames(df_f)[5] <- ("gene")
  write.table(df_f, file = file.path(x, "diff_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  return(df_f)
})

names(ls_vfdb_k2044) <- dirs

ls_vfdb_78578 <- lapply(dirs, function(x) {
  f <- data.table::fread(file.path(x, "vfdb_results.txt"))
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
  diff_gene <- data.table::fread(file.path(x, "feature_diff_gene.txt")) %>% as.data.frame()
  diff_gene$V4 <- gsub("_\\d+$", "", diff_gene$V4)
  rownames(diff_gene) <- diff_gene$V5
  inter <- intersect(diff_gene$V5, df_f$id)
  df_f <- cbind(df_f[inter,], diff_gene[inter, "V4"])
  colnames(df_f)[5] <- ("gene")
  write.table(df_f, file = file.path(x, "diff_genes_78578.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  return(df_f)
})

names(ls_vfdb_78578) <- dirs

library(ggplot2)

lapply(dirs, function(x){
  
  categories <- ls_vfdb_78578[[x]][["category_vfdb"]]
  
  # 创建数据框，计算每个类别的数量
  category_counts <- as.data.frame(table(categories))
  
  # 绘制柱形图
  ggplot(category_counts, aes(x = categories, y = Freq, fill = categories)) +
    geom_bar(stat = "identity") +
    labs(x = "Category", y = "Count", title = "MGH 78578") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set3")  # 使用预设调色板
  
  ggsave(file = file.path(x, "category_counts.pdf"), width = 8, height = 6)
})


lapply(dirs, function(x){
  
  categories <- ls_vfdb_k2044[[x]][["category_vfdb"]]
  
  # 创建数据框，计算每个类别的数量
  category_counts <- as.data.frame(table(categories))
  
  # 绘制柱形图
  ggplot(category_counts, aes(x = categories, y = Freq, fill = categories)) +
    geom_bar(stat = "identity") +
    labs(x = "Category", y = "Count", title = "NUTH-K2044") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    scale_fill_brewer(palette = "Set3")  # 使用预设调色板
  
  ggsave(file = file.path(x, "category_counts_k2044.pdf"), width = 8, height = 6)
})

x <- function(i){
  dg <- data.table::fread(file.path(i, "diff_genes.txt"))
  cat("The differentially expressed genes in ", i, " is:\n ", paste0(unique(dg[[5]]), collapse = "，"), "\n")
  dg_78578 <- data.table::fread(file.path(i, "diff_genes_78578.txt"))
  cat("The differentially expressed genes of 78578 in ", i, " is:\n ", paste0(unique(dg_78578[[5]]), collapse = "，"), "\n")
}
x("HVKP5")
x("HVKP6")
x("HVKP7")
x("HVKP8")
x("HVKP17")
x("HVKP18")
x("HVKP19")
x("HVKP20")
x("HVKP21")
x("HVKP22")
x("HVKP23")

# 加载必要的库
library(utils)

# 定义文件名和文件夹
files <- list.files("./", pattern = "tar.gz")
folders <- dir("./", pattern = "CKP|HVKP", include.dirs = TRUE)
folders <- folders[file.info(folders)$isdir]
# 批量解压缩文件
for (i in seq_along(files)) {
  # 解压到指定文件夹
  untar(files[i], exdir = folders[i])
  
  # 搜索并重命名文件
  files_in_folder <- list.files(folders[i], full.names = TRUE, recursive = TRUE)
  
  # 找到并重命名包含"list_of_operons"的文件
  operons_file <- grep("list_of_operons", files_in_folder, value = TRUE)
  if (length(operons_file) > 0) {
    file.rename(operons_file, file.path(dirname(operons_file), "list_of_operons.txt"))
    # 移动到上层文件夹
    file.copy(file.path(dirname(operons_file), "list_of_operons.txt"), file.path(folders[i], "list_of_operons.txt"))
  }
  
  # 找到并重命名包含"ORFs_coordinates"的文件
  orfs_file <- grep("ORFs_coordinates", files_in_folder, value = TRUE)
  if (length(orfs_file) > 0) {
    file.rename(orfs_file, file.path(dirname(orfs_file), "ORFs_coordinates.txt"))
    # 移动到上层文件夹
    file.copy(file.path(dirname(orfs_file), "ORFs_coordinates.txt"), file.path(folders[i], "ORFs_coordinates.txt"))
  }
}


# 解压kegg_res.zip到kegg_res文件夹
unzip("kegg_res.zip")

# 定义目标文件和目标文件夹
files_to_copy <- list.files("kegg_res", full.names = TRUE)

# 批量复制文件到指定文件夹
for (i in seq_along(files_to_copy)) {
  # 构建源文件路径
  source_file <- files_to_copy[i]
  
  # 复制文件到目标文件夹
  file.copy(source_file, file.path(stringr::str_extract(basename(source_file), ".*(?=_kegg)"), basename(source_file)))
}



