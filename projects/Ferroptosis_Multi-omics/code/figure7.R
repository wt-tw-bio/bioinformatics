library(CellChat)
library(Seurat)
library(magrittr)
library(ggplot2)

#导出normalized data和细胞类型注释信息
data.input <- GetAssayData(sc_data_sub, assay = "RNA", layer = "data") # normalized data matrix
meta <- data.frame(row.names = rownames(sc_data_sub@meta.data),CellType = sc_data_sub$cell_type)
#创建一个Cell Chat对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
#设置cellType为默认的细胞分类方式
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "CellType")
levels(cellchat@idents)
#导入配体受体数据库
CellChatDB <- CellChatDB.mouse
#使用整个配体库来识别相互作用
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
#提取细胞通讯相关基因
cellchat <- subsetData(cellchat)
future::plan("default")
#鉴定高表达的配体和受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#储存到cellchat的cellchat@LR$LRsig中，将配体受体相关基因表达值投射到互作网络上并进行校正，时间较久
cellchat <- projectData(cellchat, PPI.mouse)
#推断细胞通讯网络，为每个通讯网络计算一个概率值和对应的p值
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) #raw.use表示使用未校正的数据
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
#统计细胞通讯数量和强度（概率）
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Number of interactions")
pdf("cellchat_circle_weight.pdf", height = 13/2.6, width = 13/2.6)
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Weight of interactions")
dev.off()

pdf("plot/py_plot/cellchat_heatmap.pdf", height = 10/2.6, width = 12/2.6)
netVisual_heatmap(cellchat, color.heatmap = "Reds")
dev.off()

pdf("plot/py_plot/cellchat_buubble.pdf", height = 18/2.6, width = 24/2.6)
netVisual_bubble(cellchat, 
                 targets.use = c(5, 6), 
                 remove.isolate = TRUE)
dev.off()

pdf("plot/py_plot/cellchat_buubble2.pdf", height = 18/2.6, width = 24/2.6)
netVisual_bubble(cellchat, 
                 sources.use = c(5, 6), 
                 remove.isolate = TRUE)
dev.off()


levels(cellchat@idents)


library(liana)
DefaultAssay(sc_immune)
Idents(sc_immune)
liana_res <- liana_wrap(sc_immune)
liana_res <- liana_aggregate(liana_res)
dplyr::glimpse(liana_res)
p_chat <- liana_res %>%
  liana_dotplot(ntop = 30, source_groups = c("Monocyte_Macrophages")) 
p_chat +
  theme(axis.text.x = element_text(angle = 90, size = 18, vjust = 0.5),
        axis.text.y = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 15))
ggsave("cell_chat/cell_chat.pdf", height = 10, width = 18)


liana_trunc <- liana_res %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

pdf("cell_chat/cell_chat_heatmap.pdf", height = 6, width = 10)
heat_freq(liana_trunc)
dev.off()


library(nichenetr)
organism = "human"

#https://lishensuo.github.io/posts/bioinfo/029%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--nichenet%E7%BB%86%E8%83%9E%E9%80%9A%E8%AE%AF%E5%88%86%E6%9E%90/

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}



nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = sc_gbm, 
  expression_pct = 0.10,
  #Group
  condition_colname = "Source", 
  condition_oi = "Tumor", condition_reference = "Tumor", 
  # receiver
  receiver = c("MES-like Malignant", "Mono/Macro"),
  #geneset = "DE", 
  #lfc_cutoff = 0.25,
  # sender
  sender = "all", 
  #top_n_ligands = 20,
  #top_n_targets = 200, 
  #cutoff_visualization = 0.33,
  # refer data
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
)


