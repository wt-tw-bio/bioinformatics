metanr_packages <- function(){
  
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

library(MetaboAnalystR)
library(ggplot2)
devtools::install_github("jaspershen/MetNormalizer")

devtools::install_github("jaspershen/demoData")

library(demoData)
library(MetNormalizer)
path <- system.file("MetNormalizer", package = "demoData")
file.copy(from = path, to = ".", overwrite = TRUE, recursive = TRUE)
new.path <- file.path("./MetNormalizer")

metNor(
  ms1.data.name = "data.csv",
  sample.info.name = "sample.info.csv",
  minfrac.qc = 0,
  minfrac.sample = 0,
  optimization = TRUE,
  multiple = 5,
  threads = 4,
  path = "./"
)

which(duplicated(expers$name))
expers <- readxl::read_xlsx("./neg.xlsx")
expers$MS1_name[expers$MS1_name == "--"] <- ""
expers$MS2_name[expers$MS2_name == "--"] <- ""
expers$MS2_name <- paste0(expers$MS2_name, expers$MS1_name)
#table(duplicated(expers$MS2_name))
#expers <- subset(expers, !duplicated(expers$MS2_name))
#expers <- expers[c(2,6,5,7:31)]
#names(expers)[1:3] <- c("name", "mz", "rt")
#sample.info <- data.frame(sample.name = colnames(expers)[4:ncol(expers)],
                          #injection.order = c(2:7, 9:16, 18:25, 1, 8, 17),
                          #injection.order = c(1:25),
                          #class = c(rep("Subject", 22), rep("QC", 3)))
#write.csv(sample.info, "sample.info.csv", row.names = FALSE)
#write.csv(expers, "data.csv", row.names = FALSE)

expers <- subset(expers, !duplicated(expers$MS2_name))
expers <- tibble::remove_rownames(expers) %>% tibble::column_to_rownames("MS2_name")

expers_hp <- expers[c(6:11, 20:27)]
colnames(expers_hp) <- gsub("-", "_", colnames(expers_hp))
group <- data.frame(sample = colnames(expers_hp), group = c(rep("PBS", 6), rep("HVKP", 8)))
res_dem <- dem(expers_hp, group, dir_path = "normalization_hp", plot = TRUE, control = "PBS", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter")

dat_dem <- res_dem$res_dem

dat_dem$group <- ifelse(dat_dem$vip >= 1 & dat_dem$p_value < 0.05 & dat_dem$logFC>=1, "Sig", "Unsig")
table(dat_dem$group)



expers_hc <- expers[c(12:19, 20:27)]
colnames(expers_hc) <- gsub("-", "_", colnames(expers_hc))
group <- data.frame(sample = colnames(expers_hc), group = c(rep("CKP", 8), rep("HVKP", 8)))
res_dem_hc <- dem(expers_hc, group, dir_path = "normalization3", plot = TRUE, control = "CKP", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter")
dat_dem_hc <- res_dem_hc$res_dem

dat_dem_hc$group <- ifelse(dat_dem_hc$vip >= 1 & dat_dem_hc$p_value < 0.05 & dat_dem$logFC>=1, "Sig", "Unsig")
table(dat_dem_hc$group)

sigs <- intersect(rownames(dat_dem)[dat_dem$group == "Sig"], rownames(dat_dem_hc)[dat_dem_hc$group == "Sig"])

total <- kegg_ref(org = "mmu")
inputs <- meta.map(sigs)
inputs <- as.data.frame(inputs)
a <- na.omit(inputs$KEGG[inputs$KEGG != "NA"])
df <- kegg_meta(a)
p_kegg_meta <- ggplot(data = df, aes(GeneRatio, Description))+
  geom_point(aes(size=Count, color=pvalue, fill=pvalue), pch=21)+
  scale_color_gradientn(colours = (rev(RColorBrewer::brewer.pal(11,"RdBu"))))+
  scale_fill_gradientn(colours =(rev(RColorBrewer::brewer.pal(11,"RdBu"))))+
  guides(size=guide_legend(title="Count"))+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(color="black", size = 12),
        axis.text.y=element_text(color="black", size = 12),
        panel.background = element_rect(fill = NA,color = NA),
        panel.border = element_rect(fill=NA,color="black",size=1,linetype="solid"),
        legend.key=element_blank(),
        legend.title = element_text(color="black",size=14),
        legend.text = element_text(color="black",size=12),
        #legend.spacing.x=unit(0.1,'cm'),
        #legend.key.width=unit(0.5,'cm'),
        #legend.key.height=unit(0.5,'cm'),
        legend.background=element_blank()
  )+
  scale_y_discrete(labels = function(y) stringr::str_wrap(y,width=30))
ggsave("meta_kegg.pdf", plot = p_kegg_meta, width = 6, height = 4)

library(ggthemes)
# 火山图
dat_dem$logP <- -log10(dat_dem$p_value)
p_vip <- ggplot(data = dat_dem,aes(x = vip,y = logP, color = group))+
  xlab("VIP score")+
  ylab("-log10(p-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(size=1,alpha=0.6) +
  scale_color_manual(values = c("Sig" = scales::hue_pal()(2)[1],"Unsig" = 'grey'))+
  geom_hline(yintercept=1.3,linetype="longdash",col="grey")+
  geom_vline(xintercept=c(1),linetype="longdash",col="grey")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.key = element_rect(colour = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", fill=NA))

ggsave("vip_hp.pdf",height = 4,width = 5, plot = p_vip)

mycol <- scales::hue_pal()(4)[c(1,4)]
# 热图
dat_norm <- res_dem$dat_norm
group <- group[match(colnames(dat_norm), group$sample), ]
all(colnames(dat_norm) == group$sample)
library(pheatmap)
library(scales)
anno <- data.frame("annotation" = as.factor(group$group))
rownames(anno) <- group$sample
annotation_col <- list(annotation = c("HVKP" = mycol[1], "PBS" = mycol[2]))
label <- arrange(subset(dat_dem, group == "Sig"), desc(vip))
label <- subset(label, logFC >= 1)
dat_heat <- dat_norm[rownames(label)[1:20], ]
p_heatmap_meta <- pheatmap(as.matrix(dat_heat),cluster_col=F,cluster_rows = F, 
                           scale = "row", 
                           cluster_cols = F,
                           show_rownames = T,show_colnames = F,
                           fontsize_col = 2, fontsize =10, fontsize_row = 12,
                           annotation_col=anno, annotation_colors = annotation_col,
                           color = colorRampPalette(c(mycol[2],'white', mycol[1]))(100),
                           border_color = NA)
pdf("heatmap_meta_hp.pdf", height = 8, width = 16)
p_heatmap_meta
dev.off()

# 火山图
dat_dem_hc$logP <- -log10(dat_dem_hc$p_value)
p_vip_hc <- ggplot(data = dat_dem_hc,aes(x = vip,y = logP, color = group))+
  xlab("VIP score")+
  ylab("-log10(p-value)")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(size=1,alpha=0.6) +
  scale_color_manual(values = c("Sig" = scales::hue_pal()(2)[1],"Unsig" = 'grey'))+
  geom_hline(yintercept=1.3,linetype="longdash",col="grey")+
  geom_vline(xintercept=c(1),linetype="longdash",col="grey")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.key = element_rect(colour = NA),
        panel.grid = element_blank(),
        panel.background = element_rect(colour = "black", fill=NA))

ggsave("vip_hc.pdf",height = 4,width = 5, plot = p_vip)


# 热图
dat_norm_hc <- res_dem_hc$dat_norm
group <- group[match(colnames(dat_norm_hc), group$sample), ]
group <- group[order(group$sample, decreasing = TRUE),]
dat_norm_hc <- dat_norm_hc[,group$sample]
all(colnames(dat_norm_hc) == group$sample)
library(pheatmap)
library(scales)

anno <- data.frame("annotation" = as.factor(group$group))
rownames(anno) <- group$sample
annotation_col <- list(annotation = c("HVKP" = mycol[1], "CKP" = mycol[2]))
label <- arrange(subset(dat_dem_hc, group == "Sig"), desc(vip))
label <- subset(label, logFC >= 1)
dat_heat <- dat_norm_hc[rownames(label)[1:20], ]
p_heatmap_hc <- pheatmap(as.matrix(dat_heat),cluster_col=F,cluster_rows = F, 
                           scale = "row", 
                           cluster_cols = F,
                           show_rownames = T,show_colnames = F,
                           fontsize_col = 2, fontsize =10, fontsize_row = 12,
                           annotation_col=anno, annotation_colors = annotation_col,
                           color = colorRampPalette(c(mycol[2],'white', mycol[1]))(100),
                           border_color = NA)
pdf("heatmap_meta_hc.pdf", height = 8, width = 16)
p_heatmap_hc
dev.off()


