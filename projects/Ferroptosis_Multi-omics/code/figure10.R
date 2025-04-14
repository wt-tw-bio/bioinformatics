library(edgeR)
library(limma)
#DEGs
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(result)
# 创建 DGEList 对象
dge <- DGEList(counts = result)
# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)
# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE)
# 使用线性模型进行拟合
fit <- lmFit(v, design)
# 和上面两个包一样，需要说明是谁比谁
# 创建对比矩阵
con <- "CKP-PBS"
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
deg_cp <- topTable(fit2,coef = 1,number = Inf)

deg_cp$threshold=factor(ifelse(deg_cp$adj.P.Val < 0.05 & abs(deg_cp$logFC) >= 1, 
                               ifelse(deg_cp$logFC>= 1 ,'Up','Down'),'Not sig'),levels=c('Up','Down','Not sig'))
deg_cp$logP <- -log10(deg_cp$adj.P.Val)
deg_cp_sig <- subset(deg_cp,threshold %in% c("Up","Down"))
write.csv(deg_cp, file = "plot/bulk/ckp-pbs/deg.csv", row.names = FALSE)

mycol <- scales::hue_pal()(4)[c(1,4)]
library(ggthemes)
# 火山图
p_ckp_pbs <- ggplot(data = deg_cp,aes(x = logFC,y = logP,color = threshold))+
  xlab("log2 Fold Change")+
  ylab("-log10(p-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(size=2,alpha=0.6) +
  scale_color_manual(values = c("Up" = mycol[1], "Not sig"= 'grey', 'Down' = mycol[2]))+
  geom_hline(yintercept=1.3,linetype="longdash",col="grey")+
  geom_vline(xintercept=c(-1,1),linetype="longdash",col="grey")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.key = element_rect(colour = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", fill=NA))
#theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5),legend.position = 'none')
label <- deg_cp[intersect(rownames(deg_cp_sig), c(FRGs$driver$MGI.symbol, FRGs$suppressor$MGI.symbol, FRGs$marker$MGI.symbol)),]
volcano_plot <- p_ckp_pbs +
  geom_point(size = 5, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color="black",label.size = 0.5,
                            cex = 4,max.overlaps = 30)
ggsave("plot/bulk/ckp-pbs/vlcano_ckp-pbs.pdf",height = 6,width = 8, plot = volcano_plot)



library(clusterProfiler)
library(enrichplot)
GSEA_input <- deg_cp$logFC
names(GSEA_input) <- rownames(deg_cp)
GSEA_input <- sort(GSEA_input, decreasing = TRUE)
res_gsea <- GSEA(GSEA_input, #排序后的gene
                 TERM2GENE = df, #基因集
                 pvalueCutoff = 0.05, #P值阈值
                 minGSSize = 2, #最小基因数量
                 maxGSSize = 500, #最大基因数量
                 eps = 0, #P值边界
                 pAdjustMethod = "BH") #校正P值的计算方法

result_gsea <- data.frame(res_gsea)


gene_symbol <- bitr(rownames(deg_cp_sig),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
KEGG1<-enrichKEGG(gene_symbol$ENTREZID,#KEGG富集分析
                  organism  = "mmu",pAdjustMethod="none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1)

library('org.Mm.eg.db')

GO1<-enrichGO(gene_symbol$ENTREZID,#GO富集分析
              OrgDb ='org.Mm.eg.db',pAdjustMethod="none",
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 1,#设定q值阈值
              readable = T)

GO_res <- GO1@result
Kegg_res <- KEGG1@result


write.csv(result_gsea, file = "plot/bulk/ckp-pbs/gsea_all.csv", row.names = FALSE)
write.csv(GO_res, file = "plot/bulk/ckp-pbs/GO.csv", row.names = FALSE)
write.csv(Kegg_res, file = "plot/bulk/ckp-pbs/Kegg.csv", row.names = FALSE)

res_gsea_fer <- GSEA(GSEA_input, #排序后的gene
                     TERM2GENE = gsea_FRGs, #基因集
                     pvalueCutoff = 0.05, #P值阈值
                     minGSSize = 2, #最小基因数量
                     maxGSSize = 500, #最大基因数量
                     eps = 0, #P值边界
                     pAdjustMethod = "BH") #校正P值的计算方法

result_gsea_fer <- data.frame(res_gsea_fer)
p_gsea <- gsea_plot(res_gsea_fer, "suppressor")
ggsave("plot/bulk/ckp-pbs/gsea.pdf", width = 7, height = 6, plot = p_gsea)

debugonce(gsea_plot)

# 热图
library(pheatmap)
group_heat <- group[grep("PBS|CKP", group)]
anno <- data.frame("annotation" = as.factor(group_heat))
rownames(anno) <- colnames(result)[grep("PBS|CKP", colnames(result))]
annotation_col <- list(annotation = c("PBS" = mycol[1], "CKP" = mycol[2]))
label <- arrange(label,desc(logFC))
dat_heat <- v[["E"]][rownames(label),rownames(anno)]
p <- pheatmap(dat_heat,cluster_col=F,cluster_rows = F, 
              scale = "row", 
              cluster_cols = F,
              show_rownames = T,show_colnames = F,
              fontsize_col = 2, fontsize =10,fontsize_row = 12,
              annotation_col=anno,annotation_colors = annotation_col,
              color = colorRampPalette(c(mycol[2],'white',mycol[1]))(100),
              border_color = NA)
pdf("plot/bulk/ckp-pbs/heatmap_ckp-pbs.pdf",height = 4,width = 6)
p
dev.off()
