library("dittoSeq")

cell_percent <- dittoBarPlot(sc_data_sub, "group", group.by = "cell_type",main = '') +
  xlab('') +
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(family = 'serif',size = 12, colour = "black"),
        text = element_text(family = 'serif',size = 12),
        legend.position = 'top') 
ggsave("plot/py_plot/sc_percent.pdf", height = 6, width = 12, plot = cell_percent)

kup_data <- subset(sc_data_sub, cell_type == "Kupffer cells")
Idents(kup_data) <- kup_data$group
deg_kup <- FindAllMarkers(kup_data)

res_deg_gk <- deg_kup %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  as.data.frame()

frgs <- purrr::map(FRGs, list(2)) %>% unlist
inter_frgs <- intersect(frgs, res_deg_gk$gene)


sc_volcano <- scRNAtoolVis::jjVolcano(diffData = res_deg_gk, aesCol = rev(scales::hue_pal()(2)))+
  xlab("")+
  theme(legend.position = "right")
ggsave("plot/py_plot/sc_volcano.pdf", height = 6, width = 8, plot = sc_volcano)  

#pseudo_data <- AggregateExpression(sc_data_sub, assays = "RNA", return.seurat = T, group.by = c("cell_type", "sample"))
#pseudo_data$group <- gsub("\\d$", "", pseudo_data$sample)
#pseudo_data$celltype.group <- paste(pseudo_data$cell_type, pseudo_data$group, sep = "_")
#Idents(pseudo_data) <- "celltype.group"
#bulk.kup.de <- FindMarkers(object = pseudo_data, 
#                           ident.1 = "Kupffer cells_HVKP", 
#                           ident.2 = "Kupffer cells_CKP",
#                           test.use = "DESeq2")
#head(bulk.kup.de, n = 15)
#deg_bulk_up <- subset(bulk.kup.de, avg_log2FC > 0 & p_val_adj<0.05)
#deg_bulk_down <- subset(bulk.kup.de, avg_log2FC < 0 & p_val_adj<0.05)

deg_kup_top10 <- deg_kup %>%
  filter(avg_log2FC>1) %>%
  group_by(cluster) %>%
  slice_head(n = 10)


sc_heat <- DoHeatmap(kup_data, features = unique(deg_kup_top10$gene))
ggsave("plot/py_plot/sc_heat.pdf", height = 6, width = 8, plot = sc_heat)
sc_heat_avg <- scRNAtoolVis::AverageHeatmap(object = kup_data, markerGene = unique(deg_kup_top10$gene), gene.order = unique(deg_kup_top10$gene), row_title = "")
pdf("plot/py_plot/sc_heat_average.pdf", height = 6, width = 5)
sc_heat_avg
dev.off()

library(clusterProfiler)
library(org.Mm.eg.db)
gene_entre <- bitr(res_deg_gk$gene, "SYMBOL", "ENTREZID", "org.Mm.eg.db")
res_deg_gk <- merge(res_deg_gk, gene_entre, by.x = "gene", by.y = "SYMBOL")
head(res_deg_gk)
sc_kegg <- compareCluster(ENTREZID~cluster, data = res_deg_gk, fun = enrichKEGG, organism="mmu", pvalueCutoff=0.05)
sc_go <- compareCluster(ENTREZID~cluster, data = res_deg_gk, OrgDb='org.Mm.eg.db', fun = enrichGO, ont = "ALL", pvalueCutoff=0.05)
res_kegg <- sc_kegg@compareClusterResult

res_kegg_sub <- res_kegg %>%
  group_by(Description) %>%
  filter(n() == 1) %>%
  ungroup()

res_kegg_sub <- subset(res_kegg, !grepl("Infectious disease|Translation|Development and regeneration|Immune disease|Digestive system|Signaling molecules and interaction|Cancer|Neurodegenerative disease|Folding, sorting and degradation|Endocrine and metabolic disease|Excretory system|Circulatory system|Endocrine system|Cellular community - eukaryotes", subcategory))
res_kegg_sub <- subset(res_kegg_sub, !grepl("Neutrophil extracellular trap formation|Platelet activation|C-type lectin receptor signaling pathway|Leukocyte transendothelial migration|Fc gamma R-mediated phagocytosis|Antigen processing and presentation|Hematopoietic cell lineage|Neutrophil extracellular trap formation|Thermogenesis|Regulation of actin cytoskeleton|Viral myocarditis", Description))
res_kegg_sub$Description <- stringr::str_extract(res_kegg_sub$Description, pattern = ".*(?= - Mus)")
sc_kegg_sub <- sc_kegg
sc_kegg_sub@compareClusterResult <- res_kegg_sub
library(enrichplot)
p_dot_kegg <- dotplot(sc_kegg_sub, showCategory = 20, label_format = 70)+ xlab("")

ggsave("plot/py_plot/sc_filter_kegg.pdf", height = 8, width = 8, plot = p_dot_kegg)

res_go <- sc_go@compareClusterResult
res_go_sub <- res_go %>%
  group_by(Description) %>%
  filter(n() == 1) %>%
  ungroup()
sc_go_sub <- sc_go
sc_go_sub@compareClusterResult <- res_go_sub
p_dot_go <- dotplot(sc_go_sub, showCategory = 10, label_format = 70)
ggsave("plot/py_plot/sc_filter_go.pdf", height = 8, width = 8, plot = p_dot_go)



gene_list <-list(ferroptosis = frgs)
gene_list2 <- purrr::map(FRGs, list(2))
#AUCell
library(AUCell)
#AUcell是分析感兴趣的基因集在所有细胞是否存在富集
#细胞内基因排序
data_auc <- GetAssayData(sc_data_sub[["RNA"]], layer = "data")
cells_rankings <- AUCell_buildRankings(data_auc, splitByBlocks=TRUE)  #基因排序，使用10个核，加速计算
#计算每个细胞特定基因集的富集分数
cells_AUC <- AUCell_calcAUC(gene_list, cells_rankings, nCores = 10)
#提取AUC值
sc_data_sub$AUC <- as.numeric(getAUC(cells_AUC))
cells_AUC <- AUCell_calcAUC(gene_list2, cells_rankings, nCores = 10)
#提取AUC值
sc_data_sub@meta.data <- cbind(sc_data_sub@meta.data, as.data.frame(t(getAUC(cells_AUC)))) %>% as.data.frame()
#绘制细胞AUC值的umap分布图
plot.df<- data.frame(sc_data_sub@meta.data, sc_data_sub@reductions$umap@cell.embeddings)

class_avg <- plot.df %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    UMAP_1 = median(umap_1),
    UMAP_2 = median(umap_2)
  )

aucell_p1 <- ggplot() + 
  geom_point(data=plot.df, aes(x=umap_1,y=umap_2,colour=AUC), size = 1.5) +
  viridis::scale_color_viridis(option="F") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)+
  ggrepel::geom_label_repel(aes(label = cell_type),
                            data = class_avg,
                            size = 4,
                            x = class_avg$UMAP_1,
                            y = class_avg$UMAP_2,
                            label.size = 0,
                            segment.color = NA
  )
ggsave("AUCell.pdf", plot = p1, scale = 1, width = 12, height =10, units =c("cm"))

# 假设你的数据框为data，分组变量为group
groups <- unique(plot.df$cell_type)
# 生成所有可能的组合
all_pairs <- combn(groups, 2)
index <- apply(all_pairs, 2, function(x){grepl("Kupffer cells",x) %>% any()})
ls_pairs <- all_pairs[, index] %<>% as.data.frame() %>% as.list()

#细胞的AUC箱线图
library(ggpubr)
aucell_p1 <- ggboxplot(plot.df,x = 'cell_type',y = 'AUC',color = 'cell_type',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))



aucell_p2 <- ggboxplot(plot.df,x = 'group',y = 'AUC',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))



aucell_p3 <- ggboxplot(plot.df,x = 'group',y = 'driver',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))

aucell_p4 <- ggboxplot(plot.df,x = 'group',y = 'suppressor',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))

aucell_p5 <- ggboxplot(plot.df,x = 'group',y = 'marker',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))


plot.df.sub <- subset(plot.df, cell_type == "Kupffer cells")


aucell_p6 <- ggboxplot(plot.df.sub,x = 'group',y = 'AUC',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  labs(title = "Kupffer cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5))



aucell_p7 <- ggboxplot(plot.df.sub,x = 'group',y = 'driver',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  labs(title = "Kupffer cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5))


aucell_p8 <- ggboxplot(plot.df.sub,x = 'group',y = 'suppressor',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  labs(title = "Kupffer cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5))


aucell_p9 <- ggboxplot(plot.df.sub,x = 'group',y = 'marker',color = 'group',palette = 'jco') +
  stat_compare_means(method = "anova")+
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none')+
  rotate_x_text(45)+
  labs(title = "Kupffer cells")+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5))




dir.create("plot/py_plot/AUC")
dir.create("plot/py_plot/AUC/all")
dir.create("plot/py_plot/AUC/kupffer")
ggsave("plot/py_plot/AUC/all/AUCell.pdf", plot = aucell_p1, scale = 1, width = 16, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/all/AUCell_sample.pdf", plot = aucell_p2, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/all/AUCell_driver.pdf", plot = aucell_p3, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/all/AUCell_suppressor.pdf", plot = aucell_p4, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/all/AUCell_marker.pdf", plot = aucell_p5, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/kupffer/AUCell_sample.pdf", plot = aucell_p6, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/kupffer/AUCell_driver.pdf", plot = aucell_p7, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/kupffer/AUCell_suppressor.pdf", plot = aucell_p8, scale = 1, width = 12, height =10, units =c("cm"))
ggsave("plot/py_plot/AUC/kupffer/AUCell_marker.pdf", plot = aucell_p9, scale = 1, width = 12, height =10, units =c("cm"))
