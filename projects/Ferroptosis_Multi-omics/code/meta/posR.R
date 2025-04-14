expers_pos <- readxl::read_xlsx("./pos.xlsx")
expers_pos$MS1_name[expers_pos$MS1_name == "--"] <- ""
expers_pos$MS2_name[expers_pos$MS2_name == "--"] <- ""
expers_pos$MS2_name <- paste0(expers_pos$MS2_name, expers_pos$MS1_name)
expers_pos <- subset(expers_pos, !duplicated(expers_pos$MS2_name))
expers_pos <- tibble::remove_rownames(expers_pos) %>% tibble::column_to_rownames("MS2_name")

expers_pos_hp <- expers_pos[c(6:11, 20:27)]
colnames(expers_pos_hp) <- gsub("-", "_", colnames(expers_pos_hp))
group <- data.frame(sample = colnames(expers_pos_hp), group = c(rep("PBS", 6), rep("HVKP", 8)))
res_dem_pos <- dem(expers_pos_hp, group, dir_path = "normalization_hp_pos", plot = TRUE, control = "PBS", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter")

dat_dem_pos <- res_dem_pos$res_dem

dat_dem_pos$group <- ifelse(dat_dem_pos$vip >= 1 & dat_dem_pos$p_value < 0.05 & dat_dem_pos$logFC>=1, "Sig", "Unsig")
table(dat_dem_pos$group)

expers_pos_hc <- expers_pos[c(12:19, 20:27)]
colnames(expers_pos_hc) <- gsub("-", "_", colnames(expers_pos_hc))
group <- data.frame(sample = colnames(expers_pos_hc), group = c(rep("CKP", 8), rep("HVKP", 8)))
res_dem_pos_hc <- dem(expers_pos_hc, group, dir_path = "normalization_hc_pos", plot = TRUE, control = "CKP", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter")
dat_dem_hc_pos <- res_dem_pos_hc$res_dem

dat_dem_hc_pos$group <- ifelse(dat_dem_hc_pos$vip >= 1 & dat_dem_hc_pos$p_value < 0.05 & dat_dem_hc_pos$logFC>=1, "Sig", "Unsig")
table(dat_dem_hc_pos$group)

sigs_pos <- intersect(rownames(dat_dem_pos)[dat_dem_pos$group == "Sig"], rownames(dat_dem_hc_pos)[dat_dem_hc_pos$group == "Sig"])


inputs <- meta.map(c(sigs_pos, sigs))
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


# 热图
dat_norm <- res_dem$dat_norm
group <- group[match(colnames(dat_norm), group$sample), ]
all(colnames(dat_norm) == group$sample)
library(pheatmap)
library(scales)
anno <- data.frame("annotation" = as.factor(group$group))
rownames(anno) <- group$sample
annotation_col <- list(annotation = c("HVKP" = "red", "PBS" = "blue"))
label <- arrange(subset(dat_dem, group == "Sig"), desc(vip))
label <- subset(label, logFC >= 1)
dat_heat <- dat_norm[rownames(label)[1:20], ]
p_heatmap_meta <- pheatmap(as.matrix(dat_heat),cluster_col=F,cluster_rows = F, 
                           scale = "row", 
                           cluster_cols = F,
                           show_rownames = T,show_colnames = F,
                           fontsize_col = 2, fontsize =10, fontsize_row = 12,
                           annotation_col=anno, annotation_colors = annotation_col,
                           color = colorRampPalette(c(hue_pal()(2)[2],'white', hue_pal()(2)[1]))(100),
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
annotation_col <- list(annotation = c("HVKP" = "red", "CKP" = "blue"))
label <- arrange(subset(dat_dem_hc, group == "Sig"), desc(vip))
label <- subset(label, logFC >= 1)
dat_heat <- dat_norm_hc[rownames(label)[1:20], ]
p_heatmap_hc <- pheatmap(as.matrix(dat_heat),cluster_col=F,cluster_rows = F, 
                         scale = "row", 
                         cluster_cols = F,
                         show_rownames = T,show_colnames = F,
                         fontsize_col = 2, fontsize =10, fontsize_row = 12,
                         annotation_col=anno, annotation_colors = annotation_col,
                         color = colorRampPalette(c(hue_pal()(2)[2],'white', hue_pal()(2)[1]))(100),
                         border_color = NA)
pdf("heatmap_meta_hc.pdf", height = 8, width = 16)
p_heatmap_hc
dev.off()