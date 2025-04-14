# 加载必要的包
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)

# 创建模拟数据
set.seed(123)

# 创建类别、通路和基因的数据
categories <- c("类别1", "类别2", "类别3", "类别4")

# 通路富集数据
pathway_data <- data.frame(
  Category = categories,
  Pathway = c("细胞周期", "代谢通路", "信号转导", "免疫反应"),
  EnrichmentScore = c(3.2, 2.8, 3.5, 2.5),
  PValue = c(0.001, 0.005, 0.0008, 0.01)
)

# 为计算和显示添加-log10(p-value)
pathway_data$NegLogPValue <- -log10(pathway_data$PValue)

# 基因数据
gene_list <- list(
  "类别1" = c("a", "b", "c", "d", "e"),
  "类别2" = c("f", "g", "h", "i"),
  "类别3" = c("j", "k", "l", "m", "n", "o"),
  "类别4" = c("p", "q", "r", "s")
)

# 特殊基因
special_genes <- c("a", "h", "l", "p")

# 创建完整的基因数据框
gene_data <- data.frame(
  Category = rep(names(gene_list), sapply(gene_list, length)),
  Gene = unlist(gene_list),
  Expression = runif(sum(sapply(gene_list, length)), 1, 10)
)

# 标记特殊基因
gene_data$Special <- gene_data$Gene %in% special_genes

# 为类别分配颜色
category_colors <- brewer.pal(n = length(categories), name = "Set2")
names(category_colors) <- categories

# 1. 创建KEGG富集柱状图
create_enrichment_plot <- function(pathway_data, category_colors) {
  # 确保类别顺序一致
  pathway_data$Category <- factor(pathway_data$Category, levels = categories)
  
  # 创建富集柱状图
  p <- ggplot(pathway_data, aes(x = Category, y = NegLogPValue, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Pathway), vjust = -0.5, size = 3.5) +
    geom_text(aes(label = sprintf("p = %.4f", PValue)), vjust = 1.5, size = 3, color = "white") +
    scale_fill_manual(values = category_colors) +
    labs(
      title = "各类别KEGG通路富集分析",
      x = NULL,
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

# 2. 创建基因表格函数 - 每个类别一个表格
create_gene_tables <- function(gene_data, categories, special_genes) {
  # 创建一个空列表存储每个类别的表格
  table_list <- list()
  
  # 对每个类别创建表格
  for(cat in categories) {
    # 过滤该类别的基因
    cat_genes <- gene_data %>%
      filter(Category == cat)
    
    # 准备表格数据
    table_data <- cat_genes %>%
      select(Gene)
    
    # 创建表格
    # 获取哪些行是特殊基因
    special_rows <- which(cat_genes$Gene %in% special_genes)
    
    # 创建背景颜色向量
    bg_colors <- rep("white", nrow(table_data))
    bg_colors[special_rows] <- "mistyrose"
    
    # 创建表格
    tbl <- tableGrob(
      table_data, rows = NULL,
      theme = ttheme_minimal(
        core = list(
          fg_params = list(hjust = 0, x = 0.5, fontsize = 10),
          bg_params = list(fill = bg_colors)
        ),
        colhead = list(
          fg_params = list(hjust = 0.5, fontsize = 11, fontface = "bold"),
          bg_params = list(fill = category_colors[cat], alpha = 0.5)
        )
      ),
      widths = unit(2, "cm")
    )
    table_list[[cat]] <- tbl
  }
  
  return(table_list)
}

# 3. 创建说明文本
create_legend <- function() {
  textGrob("* 浅红色背景标记的基因为特殊关注基因", 
           gp = gpar(fontsize = 10, col = "black"), 
           just = "left", x = 0.05)
}

# 4. 创建最终可视化 - 将富集图和表格对齐
create_aligned_visualization <- function() {
  # 创建富集图
  enrichment_plot <- create_enrichment_plot(pathway_data, category_colors)
  
  # 创建基因表格
  gene_tables <- create_gene_tables(gene_data, categories, special_genes)
  
  # 创建说明
  legend <- create_legend()
  
  # 计算每个表格应该占据的宽度比例
  genes_per_category <- sapply(gene_list, length)
  width_ratios <- genes_per_category / sum(genes_per_category)
  
  # 创建图形设备并保存当前参数
  current_dev <- dev.cur()
  
  # 绘制富集图
  enrichment_plot_grob <- ggplotGrob(enrichment_plot)
  
  # 创建表格区布局
  tables_layout <- grid.layout(
    nrow = 1, 
    ncol = length(categories),
    widths = unit(width_ratios, "null")
  )
  
  # 创建表格容器
  tables_panel <- frameGrob(layout = tables_layout)
  
  # 将表格添加到容器中
  for(i in 1:length(categories)) {
    cat <- categories[i]
    tables_panel <- placeGrob(
      tables_panel,
      grob = gene_tables[[cat]],
      row = 1,
      col = i
    )
  }
  
  # 组合所有元素
  final_viz <- grid.arrange(
    enrichment_plot_grob,
    tables_panel,
    legend,
    heights = c(3, 4, 0.5),
    ncol = 1
  )
  
  # 恢复原始图形设备
  if(current_dev > 1) dev.set(current_dev)
  
  return(final_viz)
}

# 执行可视化
final_visualization <- create_aligned_visualization()

# 保存图表
# ggsave("aligned_category_pathway_gene.pdf", final_visualization, width = 12, height = 10)