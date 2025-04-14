library(ggplot2)
library(gridExtra)
library(grid)
library(data.table)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(viridis)


# Define the main function to execute the analysis and visualization
perform_analysis <- function(diff_path, ref) {
  # Load and process operon data
  operon_file <- list.files(diff_path, pattern = "^operon.*", full.names = TRUE)
  all_data <- read.table(operon_file[ref], header = FALSE, fill = TRUE, sep = "\t")
  
  pathway_data <- all_data[, c(1, 10, 11)]
  pathway_data <- pathway_data[!duplicated(pathway_data), ]
  colnames(pathway_data) <- c("Category", "Pathway", "PValue")
  pathway_data$NegLogPValue <- -log10(pathway_data$PValue)
  
  # Load and process gene data
  operon_file <- list.files(diff_path, pattern = "list.*", full.names = TRUE)
  df <- fread(operon_file, header = TRUE, fill = TRUE)
  df$Operon <- zoo::na.locf(df$Operon)
  df <- subset(df, IdGene != "")
  colnames(df)[2] <- "ID"
  
  operon_gff <- list.files(diff_path, pattern = "ORF.*", full.names = TRUE)
  operon_dat <- fread(operon_gff, fill = TRUE, header = FALSE) %>% as.data.frame()
  operon_dat$ID <- str_extract(operon_dat$V9, "(?<=ID=).*(?=;pro)")
  colnames(operon_dat)[c(1, 4, 5)] <- c("locus", "start", "end")
  
  gff <- list.files(diff_path, pattern = "feature_info.*", full.names = TRUE)
  gff_dat <- fread(gff, fill = TRUE, header = FALSE)
  colnames(gff_dat)[c(1, 2, 3)] <- c("locus", "start", "end")
  
  merge_gff <- merge(gff_dat, operon_dat[c(1, 4, 5, 10)], by = c("locus", "start", "end"))
  merge_gff <- merge(merge_gff, df[, 1:2], by = "ID")
  merge_gff <- merge_gff[!duplicated(merge_gff$ID), ]
  colnames(merge_gff)[5:6] <- c("gene", "symbol")
  
  kegg <- list.files(diff_path, pattern = "kegg", full.names = TRUE)
  kegg_id <- read.table(kegg, header = FALSE, fill = TRUE)
  colnames(kegg_id) <- c("symbol", "keggid")
  merge_gff <- merge(merge_gff, kegg_id, by = "symbol")
  gene_list <- subset(merge_gff, !is.na(gene))
  gene_list$gene <- gsub("_\\d+", "", gene_list$gene)
  gene_list <- gene_list[, c("gene", "Operon")]
  gene_list <- gene_list[!duplicated(gene_list), ]
  gene_list <- split(gene_list$gene, gene_list$Operon)
  gene_list <- gene_list[as.character(pathway_data$Category)]
  
  # Special genes
  special_genes <- unique(all_data$V6)
  
  # Create gene data frame
  gene_data <- data.frame(
    Category = rep(names(gene_list), sapply(gene_list, length)),
    Gene = unlist(gene_list)
  )
  
  # Mark special genes
  gene_data$Special <- gene_data$Gene %in% special_genes
  
  # Assign colors to categories
  categories <- unique(pathway_data$Category)
  category_colors <- viridis::viridis(length(categories))
  names(category_colors) <- categories
  
  # Create enrichment plot
  create_enrichment_plot <- function(pathway_data, category_colors) {
    pathway_data$Category <- factor(pathway_data$Category, levels = categories)
    
    p <- ggplot(pathway_data, aes(x = Category, y = NegLogPValue, fill = Category)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = Category), vjust = -0.5, size = 3.5) +
      geom_text(aes(label = sprintf("p = %.4f", PValue)), vjust = 1.5, size = 3, color = "white") +
      scale_fill_manual(values = category_colors) +
      scale_x_discrete(labels = pathway_data$Pathway) +  # 使用Pathway的值作为x轴标签
      labs(
        title = "KEGG Pathway",
        x = NULL,
        y = "-log10(p-value)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 1, colour = "black"),
        legend.position = "none"
      )
    
    return(p)
  }
  
  # Create gene tables
  create_gene_tables <- function(gene_data, categories, special_genes, category_colors) {
    max_rows <- max(sapply(gene_list, length))
    table_list <- list()
    
    for(i in 1:length(categories)) {
      cat <- categories[i]
      cat_genes <- gene_data %>% filter(Category == cat)
      table_data <- cat_genes %>% select(Gene)
      
      if(nrow(table_data) < max_rows) {
        empty_rows <- data.frame(Gene = rep("", max_rows - nrow(table_data)))
        table_data <- rbind(table_data, empty_rows)
      }
      
      special_rows <- which(table_data$Gene %in% special_genes)
      bg_colors <- rep("white", nrow(table_data))
      bg_colors[special_rows] <- "mistyrose"
      
      tbl <- tableGrob(
        table_data, 
        rows = NULL, 
        theme = ttheme_default(
          core = list(
            fg_params = list(hjust = 0.5, x = 0.5, fontsize = 10),
            bg_params = list(fill = bg_colors)
          ),
          colhead = list(
            fg_params = list(hjust = 0.5, fontsize = 11, fontface = "bold"),
            bg_params = list(fill = category_colors[as.character(cat)], alpha = 0.5)
          )
        )
      )
      
      tbl <- gtable::gtable_add_grob(
        tbl,
        grobs = rectGrob(
          gp = gpar(fill = NA, col = "black", lwd = 1)
        ),
        t = 1, b = nrow(tbl), l = 1, r = ncol(tbl)
      )
      
      table_list[[as.character(cat)]] <- tbl
    }
    
    return(table_list)
  }
  
  # Create legend
  create_legend <- function() {
    textGrob("* The genes marked in red with a light background are differentially expressed genes", 
             gp = gpar(fontsize = 10, col = "black"), 
             just = "left", x = 0.05)
  }
  
  # Create aligned visualization
  create_aligned_visualization <- function() {
    enrichment_plot <- create_enrichment_plot(pathway_data, category_colors)
    gene_tables <- create_gene_tables(gene_data, categories, special_genes, category_colors)
    width_ratios <- rep(1, length(categories))
    
    tables_row <- do.call(
      grid.arrange,
      c(gene_tables, list(ncol = length(categories), widths = width_ratios))
    )
    
    legend <- create_legend()
    
    #final_viz <- grid.arrange(
    #  enrichment_plot,
    #  tables_row,
    #  legend,
    #  heights = c(3, 4, 0.5),
    #  ncol = 1
    #)
    final_viz <- enrichment_plot / 
      tables_row / 
      legend + 
      plot_layout(heights = c(3, 4, 0.5)) &
      theme(plot.margin = margin(0, 0, 0, 0))
    
    return(final_viz)
  }
  
  # Execute visualization
  final_visualization <- create_aligned_visualization()
  return(final_visualization)
}
# Example usage:
final_plot <- perform_analysis("./HVKP1")
ggsave("Operon_kegg_plot.png", final_plot, width = 12, height = 10, units = "in", dpi = 300)


dir.create("Operon/K2044", recursive = TRUE)
for(i in folders){
  final_plot <- perform_analysis(i, ref = 2)
  ggsave(paste0("Operon/K2044/", i, ".pdf"), final_plot, width = 16, height = 12)
}


dir.create("Operon/78578", recursive = TRUE)
for(i in folders){
  final_plot <- perform_analysis(i, ref = 1)
  ggsave(paste0("Operon/78578/", i, ".pdf"), final_plot, width = 16, height = 12)
}



debugonce("perform_analysis")
