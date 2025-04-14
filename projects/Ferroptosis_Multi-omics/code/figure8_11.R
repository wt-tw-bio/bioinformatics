ls_flies <- list.files("./test/count/", full.names = TRUE)
ls_count <- lapply(ls_flies, function(x){
  data.table::fread(x) %>% as.data.frame()
})
dat_count <- data.frame(cbind(ls_count[[1]][-c(2:6)], ls_count[[2]][-c(1:6)], ls_count[[3]][-c(1:6)], ls_count[[4]][-c(1:6)]))
#all(dat_count$Geneid == dat_count$Geneid.3)
colnames(dat_count) <- gsub(".markdup.sorted.bam", "", colnames(dat_count))
anno <- data.table::fread("./test/outdir_hvkp/star_salmon/salmon.merged.gene_counts.tsv")
all(anno$gene_id %in% dat_count$Geneid)
table(duplicated(anno$gene_id))
rownames(dat_count) <- dat_count$Geneid
dat_count <- dat_count[anno$gene_id, ]
all(dat_count$Geneid == anno$gene_id)
table(duplicated(anno$gene_name))
dat_count$Geneid <- anno$gene_name
colnames(dat_count)[1] <- "gene"
dat_count_filter <- dat_count[rowSums(dat_count[-1]) > 0, ]
table(grepl("ENSMUSG", dat_count_filter$gene))
dat_count_filter <- dat_count_filter[!grepl("ENSMUSG", dat_count_filter$gene), ]

# 为每行计算除'a'列外的行和
dat_count_filter$row_sum <- rowSums(dat_count_filter[, -1])

# 按'a'列分组，选择每组中行和最大的行
result <- dat_count_filter[with(dat_count_filter, ave(row_sum, gene, FUN = max) == row_sum), ]

# 去掉辅助列
result$row_sum <- NULL

result <- result %>% tibble::remove_rownames()  %>% tibble::column_to_rownames("gene")

group <- gsub("\\d", "", colnames(result))

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
con <- "HVKP-PBS"
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
deg_hp <- topTable(fit2,coef = 1,number = Inf)

deg_hp$threshold=factor(ifelse(deg_hp$adj.P.Val < 0.05 & abs(deg_hp$logFC) >= 1, 
                                 ifelse(deg_hp$logFC>= 1 ,'Up','Down'),'Not sig'),levels=c('Up','Down','Not sig'))
deg_hp$logP <- -log10(deg_hp$adj.P.Val)
deg_hp_sig <- subset(deg_hp,threshold %in% c("Up","Down"))
write.csv(deg_hp, file = "plot/bulk/hvpk-pbs/deg.csv", row.names = FALSE)

mycol <- scales::hue_pal()(4)[c(1,4)]
library(ggthemes)
# 火山图
p_hvpk_pbs <- ggplot(data = deg_hp,aes(x = logFC,y = logP,color = threshold))+
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
label <- deg_hp[intersect(rownames(deg_hp_sig), c(FRGs$driver$MGI.symbol, FRGs$suppressor$MGI.symbol, FRGs$marker$MGI.symbol)),]
volcano_plot <- p_hvpk_pbs +
  geom_point(size = 5, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color="black",label.size = 0.5,
                            cex = 4,max.overlaps = 30)
ggsave("plot/bulk/vlcano_hvkp-pbs.pdf",height = 6,width = 8, plot = volcano_plot)



library(clusterProfiler)
library(enrichplot)
GSEA_input <- deg_hp$logFC
names(GSEA_input) <- rownames(deg_hp)
GSEA_input <- sort(GSEA_input, decreasing = TRUE)
df <- read.gmt("gsea/m2.cp.v2024.1.Mm.symbols.gmt")
res_gsea <- GSEA(GSEA_input, #排序后的gene
                 TERM2GENE = df, #基因集
                 pvalueCutoff = 0.05, #P值阈值
                 minGSSize = 2, #最小基因数量
                 maxGSSize = 500, #最大基因数量
                 eps = 0, #P值边界
                 pAdjustMethod = "BH") #校正P值的计算方法

result_gsea <- data.frame(res_gsea)


gene_symbol <- bitr(rownames(deg_hp_sig),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
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


write.csv(result_gsea, file = "plot/bulk/hvpk-pbs/gsea.csv", row.names = FALSE)
write.csv(GO_res, file = "plot/bulk/hvpk-pbs/GO.csv", row.names = FALSE)
write.csv(Kegg_res, file = "plot/bulk/hvpk-pbs/Kegg.csv", row.names = FALSE)

library(biomaRt)
listMarts()
# 人类数据库
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org/")
# 小鼠数据库
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl", 
                 host = "https://dec2021.archive.ensembl.org/")

FRGs <- lapply(list.files(pattern = "ferroptosis", full.names = TRUE), function(x){
  a <- data.table::fread(x)
  b <- a$symbol
  c <- getLDS(values = unique(b),                  # 要转换的基因集
              attributes = "hgnc_symbol",    # 需要查询的来源物种的attributes信息
              filters = "hgnc_symbol",       # 过滤参数
              mart = human,                 # 来源种属来源数据库
              martL = mouse,                # 目标种属数据库
              attributesL = "mgi_symbol",  #  需要查询的目标物种的attributes信息
              uniqueRows = TRUE)
  return(c)
}
)
names(FRGs) <- c("driver", "marker", "suppressor")
ls_FRGs <- lapply(FRGs, function(x){
  return(x[["MGI.symbol"]])
})
gsea_FRGs <- do.call(cbind, ls_FRGs) %>% reshape2::melt()
gsea_FRGs <- gsea_FRGs[-1]
colnames(gsea_FRGs) <- c("term", "gene")

res_gsea <- GSEA(GSEA_input, #排序后的gene
                 TERM2GENE = gsea_FRGs, #基因集
                 pvalueCutoff = 0.05, #P值阈值
                 minGSSize = 2, #最小基因数量
                 maxGSSize = 500, #最大基因数量
                 eps = 0, #P值边界
                 pAdjustMethod = "BH") #校正P值的计算方法

result_gsea <- data.frame(res_gsea)

p_gsea <- gsea_plot(res_gsea, "suppressor")
ggsave("plot/bulk/gsea.pdf", width = 7, height = 6, plot = p_gsea)

# 热图
library(pheatmap)
group_heat <- group[grep("HVKP|PBS", group)]
anno <- data.frame("annotation" = as.factor(group_heat))
rownames(anno) <- colnames(result)[grep("HVKP|PBS", colnames(result))]
annotation_col <- list(annotation = c("HVKP" = mycol[1], "PBS" = mycol[2]))
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
pdf("plot/bulk/heatmap_hvkp-pbs.pdf",height = 4,width = 6)
p
dev.off()


library(BayesPrism)
bk.dat <- t(result)
sc.dat <- GetAssayData(sc_data_py_sub[["RNA"]], layer = "counts") %>% as.matrix() %>% t()
dim(sc.dat)
cell.type.labels <- sc_data_py_sub$cell_type

dir.create("plot/py_plot/bulk/deconvolution", recursive = TRUE)
pdf("plot/py_plot/bulk/deconvolution/cor.pdf", height = 6,width = 10)
plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)
dev.off()


bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)


sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="mm", 
                                  #c("other_Rb","chrM","chrX","chrY","Rb","Mrp","act","hb","MALAT1")
                                  gene.group=c("chrX","chrY","Rb","Mrp","act","hb"),
                                  exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = NULL,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)


bp.res <- run.prism(prism = myPrism, n.cores=100)

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")



theta <- as.data.frame(theta)
theta$group <- gsub("\\d", "", rownames(theta))
theta_sub <- subset(theta, group %in% c("CKP", "HVKP"))
theta_box <- reshape2::melt(theta_sub,id.vars = "group")
library(ggpubr)
ggboxplot(theta_box, x="variable", y="value", fill  = "group", 
          ylab="score",
          xlab="level",
          #palette = colorplate[1:2],
          width=0.6, add = "none")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  rotate_x_text(90)+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
ggsave("figure9/theta_box.pdf",height = 5,width = 10)

plot(1:5)

dev.off()

library(ggalluvial)
theta$sample <- rownames(theta)
dat_ratio <- reshape2::melt(theta[-11], id.var = "sample")
p_ratio <- ggplot(dat_ratio, aes(x =sample, y= value, fill = variable,
                                        stratum=variable, alluvium=variable)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='',y = 'Ratio')+
  theme(axis.text = element_text(size = 12, color = "black"))
ggsave("plot/bulk/cell_ratio_bulk.pdf", height = 5, width =14, plot = p_ratio)


library(GSVA)
inter_score <- data.frame("gene" = myoverlap)
ARG.score <- gsva(expr = as.matrix(mdat),
                  gset.idx.list = inter_score,
                  method = "ssgsea", 
                  kcdf = 'Gaussian', 
                  abs.ranking = TRUE)


ssgsea <- ssgseaParam(
  v[["E"]],
  purrr::map(FRGs, list(2)))

res_gsva <- gsva(
  ssgsea
)

box_data <- as.data.frame(t(res_gsva))
box_data$group <- gsub("\\d$", "", rownames(box_data))
box_data <- reshape2::melt(box_data,id.vars = "group")
library(ggplot2)
library(ggpubr)
ggboxplot(box_data, x="variable", y="value", fill  = "group", 
          ylab="Expression",
          xlab="Gene",
          #palette =  colorplatte,
          width=0.6, add = "none")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  rotate_x_text(90)+
  stat_compare_means(aes(group=group),
                     method="kruskal.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")
