library(edgeR)
library(limma)
group <- gsub("\\d", "", colnames(bulk_count))
#DEGs
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(bulk_count)
# 创建 DGEList 对象
dge <- DGEList(counts = bulk_count)
# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)
# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE)
bulk_norm <- v$E
meta_dat <- expers[, c(12:27, 6:11)]
colnames(meta_dat) <- gsub("-", "", colnames(meta_dat))
meta_group <- data.frame(sample = colnames(meta_dat), group = rep(c("treat", "normal"), c(16, 6)))
meta_norm <- dem(meta_dat, meta_group, dir_path = "normalization_all", plot = TRUE, control = "normal", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter")
meta_norm_dat <- meta_norm$dat_norm
colnames(meta_norm_dat)
colnames(bulk_norm)
bulk_norm <- bulk_norm[, c(17:22, 1:16)]
meta_norm_sigs <- meta_norm_dat[sigs, ]
# 计算相关性和 p 值
Meta_Rna_cor <- WGCNA::cor(t(meta_norm_sigs), t(bulk_norm))
n_samples <- ncol(meta_norm_sigs)
Meta_Rna_cor.p <- WGCNA::corPvalueStudent(Meta_Rna_cor, n_samples)


# 重塑数据
Meta_Rna_melt <- reshape2::melt(Meta_Rna_cor)
Meta_Rna_melt.p <- reshape2::melt(Meta_Rna_cor.p)
Meta_Rna_cor_all <- Meta_Rna_melt
Meta_Rna_cor_all$p.value <- Meta_Rna_melt.p$value
colnames(Meta_Rna_cor_all)[1:3] <- c("meta", "gene", "cor")
res_filter <- subset(Meta_Rna_cor_all, cor >= 0.8 & p.value < 0.05) %>% na.omit()
node <- data.frame(node1 = res_filter$meta, node2 = res_filter$gene)
attr <- rbind(data.frame(node = unique(c(node$node2)), type = "gene"), data.frame(node = unique(c(node$node1)), type = "meta"))
write.table(node, "correlation.txt", row.names = F, quote = F, sep = "\t")
write.table(attr, "attr.txt", row.names = F, quote = F, sep = "\t")
