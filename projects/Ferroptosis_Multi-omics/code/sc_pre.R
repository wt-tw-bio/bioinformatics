library(Seurat)
library(reticulate)
sc <- import("scanpy", convert = TRUE)
adata = sc$read_h5ad("/work/run/projects/wtao/czh/P24052701/projects/results/CZH/anndata_cluster_harmony.h5ad")
a = adata$raw$X
# Get the expression matrix
exprs <- Matrix::t(adata$raw$X)
colnames(exprs) <- adata$raw$obs_names$to_list()
rownames(exprs) <- adata$raw$var_names$to_list()
# Create the Seurat object
sc_data <- CreateSeuratObject(exprs, meta.data = adata$obs)
# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
sc_data@reductions[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")
DimPlot(sc_data, group.by = c("leiden_2", "celltypist"), label = TRUE)
DimPlot(sc_data, group.by = c("group"), label = TRUE)
DotPlot(sc_data, features = c("Clec4f"), group.by = "leiden_0_25")
FeaturePlot(sc_data, features = c("Clec4f"))
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
all_gene <- rownames(sc_data)
# z-score:基因的均值为0, 方差为1
sc_data <- ScaleData(sc_data, features = all_gene)
# singleR自动注释
library(SingleR)
igd <- celldex::ImmGenData()
mrd <- celldex::MouseRNAseqData()
testdata <- GetAssayData(sc_data[["RNA"]], layer = "data")
#可选不同分辨率的cluster
clusters <- sc_data$leiden_0_25

cellpred <- SingleR(test = testdata,
                    ref = igd,
                    labels = igd$label.main,
                    clusters = clusters, 
                    assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred),
                      cell=cellpred$labels,
                      stringsAsFactors = F)

#sc_data$celltype <- plyr::mapvalues(x = as.integer(sc_data$cluster), from = clusterids, to = celltype)
sc_data$cell_type <- factor(sc_data$leiden_0_25, labels = as.character(celltype$cell))
DimPlot(sc_data, group.by = "cell_type", label = TRUE)

DotPlot(sc_data, features = c("Csf3r", "S100a8",
                                     "Alb",
                                     "Nkg7","Ccl5", "Il2rb",
                                     "Igfbp7", "Ptprb", "Clec4g", 
                                     "C1qb","C1qa", "C1qc", "Clec4f",
                                     "Cd68","Mpeg1", 
                                     "Cd3d", "Cd3e", "Cd3g", 
                                     "Cd19", "Cd22", "Cd79a",
                                     "Cd63", "Gata2", "Ccl3"), group.by = "leiden_2")

DimPlot(sc_data, group.by = "leiden_0_25", label = TRUE)

sc_data_sub <- subset(sc_data, `leiden_0_25` != 9)


sc_data_sub$cell_type <- "1"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(3)] <- "Hepatocytes"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(4)] <- "NK cells"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(1)] <- "NKT"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(2)] <- "T cells"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(0)] <- "B cells"
sc_data_sub$cell_type[sc_data_sub$leiden_2 %in% c(15)] <- "Kupffer cells"
sc_data_sub$cell_type[sc_data_sub$leiden_2 %in% c(20)] <- "Macrophages"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(5)] <- "Neutrophils"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(7)] <- "Basophils"
sc_data_sub$cell_type[sc_data_sub$leiden_0_25 %in% c(8)] <- "Endothelial cells"

sc_data_sub <- subset(sc_data_sub, `cell_type` != "1")

DimPlot(sc_data_sub, group.by = "cell_type", label = TRUE)
DimPlot(sc_data_sub, group.by = "leiden_0_25", label = TRUE)


p_dot_cell_type <- DotPlot(sc_data_sub, features = c("Csf3r", "S100a8",
                                                  "Alb",
                                                  "Nkg7","Ccl5", "Il2rb",
                                                  "Igfbp7", "Ptprb", "Clec4g", 
                                                  "C1qb","C1qa", "C1qc", "Clec4f",
                                                  "Cd68","Mpeg1", 
                                                  "Cd3d", "Cd3e", "Cd3g", 
                                                  "Cd19", "Cd22", "Cd79a",
                                                  "Cd63", "Gata2", "Ccl3"), group.by = "cell_type")
p_dim_cell_type <- DimPlot(sc_data_sub, group.by = "cell_type")

ggsave("plot/cell_type_dot.pdf", height = 5, width = 16, plot =p_dot_cell_type)
ggsave("plot/cell_type_dim.pdf", height = 6, width = 6, plot =p_dim_cell_type)