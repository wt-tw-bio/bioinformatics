library(KEGGREST)
library(xml2)
library(XML)
library(RCurl)
library(MetaboAnalystR)
library(magrittr)
library(dplyr)
library(tibble)

dem <- function(exprs, group,  dir_path, plot = FALSE, control = "normal", rowNorm = "SumNorm", transNorm = "LogNorm", scaleNorm = "MeanCenter") {
  
  # 检查列名是否对应
  if (!all(colnames(exprs) %in% group$sample)) {
    stop("exprs 列名和 group 的样本名称不对应")
  }
  
  # 确保列名顺序一致
  exprs <- exprs[, group$sample]
  
  # 将所有列转换为数值类型并去除缺失值
  exprs <- exprs %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  
  # 检查并创建目录
  if (!dir.exists(dir_path)) {
    dir.create(dir_path)
  } else {
    message("目录 'normalization' 已存在。创建 'normalization_meta'。")
    dir_path <- "normalization_meta"
    dir.create(dir_path)
  }
  
  # 添加分组信息并保存为 CSV 文件
  exprs_with_group <- cbind(group = group$group, t(exprs)) %>% as.data.frame() %>% rownames_to_column("sample")
  write.csv(exprs_with_group, file = file.path(dir_path, "norm_input.csv"), row.names = FALSE)
  
  # 初始化数据对象
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet <- Read.TextData(mSet, file.path(dir_path, "norm_input.csv"), "rowu", "disc")
  mSet <- SanityCheckData(mSet) %>% ReplaceMin()
  mSet <- FilterVariable(mSet, qc.filter = "F", rsd = 25, var.filter = "iqr", 
                         var.cutoff = -1, int.filter = "mean", int.cutoff = 0) %>% 
    PreparePrenormData()
  
  # 数据归一化
  mSet <- Normalization(mSet, rowNorm, transNorm, scaleNorm, ratio = FALSE, ratioNum = 20)
  if(control == levels(mSet$dataSet$cls)[1]){
    cmp.type = 1
  } else {
    cmp.type = 0
  }
  mSet <- FC.Anal(mSet, fc.thresh = 1, cmp.type = cmp.type, paired=FALSE)
  mSet <- Ttests.Anal(mSet, F, threshp = 0.05, paired = FALSE, equal.var = TRUE, pvalType = "fdr", all_results = TRUE)
  mSet <- OPLSR.Anal(mSet, reg=TRUE)
  dat_dem <- data.frame(vip = mSet$analSet$oplsda$vipVn, row.names = names(mSet$analSet$oplsda$vipVn)) %>%
    mutate(FC = mSet$analSet$fc$fc.all,
           logFC = mSet$analSet$fc$fc.log,
           p_value = mSet$analSet$tt$p.value,
           p_fdr = mSet$analSet$tt$fdr.p)
  dat_norm <- mSet$dataSet$norm
  
  # 可选绘图
  if (plot) {
    mSet <- PlotNormSummary(mSet, file.path(dir_path, "norm_summary_"), "pdf", dpi = "")
    mSet <- PlotSampleNormSummary(mSet, file.path(dir_path, "norm_sample_"), "pdf", dpi = "")
  }
  
  # 清理对象
  rm(mSet)
  
  # 返回归一化后的数据框
  return(list(dat_norm = as.data.frame(t(dat_norm)),
              res_dem = dat_dem))
}

kegg_ref <- function(org){
  #kegg_total
  kegg_total <- function(x){
    res <- keggGet(x, "kgml")
    xml_doc <- read_xml(res)
    entries <- xml_find_all(xml_doc, "//entry")
    graphics <- xml_find_all(xml_doc, "//graphics")
    # 获取每个 'entry' 元素的 id 和内容
    entry_names <- xml_attr(entries, "name")
    entry_types <- xml_attr(entries, "type")
    first_graphics_list <- lapply(entries, function(entry) {
      graphics_nodes <- xml_find_all(entry, "./graphics")
      if (length(graphics_nodes) > 0) {
        return(xml_attr(graphics_nodes[[1]], "name"))
      } else {
        return(xml_attr(graphics_nodes, "name"))
      }
    })
    graphics_names <- unlist(first_graphics_list)
    tryCatch({b <- data.frame(id = entry_names, type = entry_types, ID = graphics_names)
    c <- subset(b, type %in% c("compound", "gene"))
    c$ID <- gsub("\\.\\.\\.", "", c$ID)
    c$pathway <- xml_attr(xml_doc, "title")
    total <- c[c("pathway", "ID", "type")]
    return(total)
    },error = function(e) {
      message("An error occurred: ", x, e$message)
      total <- list(id = entry_names, type = entry_types, ID = graphics_names)
      return(total)
    })
  }
  
  org <- paste0("https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=", org)
  orgDoc <- htmlParse(getURL(org), encoding = "UTF-8")
  ul_nodes <- getNodeSet(orgDoc, "//ul")
  idx <- xpathSApply(ul_nodes[[1]], ".//a", function(a) {
    sub("/pathway/", "", xmlGetAttr(a, "href"))
  })
  ls_res <- lapply(idx, kegg_total)
  total <- do.call(rbind, ls_res)
  total <- subset(total, type == "compound")
  total <- total %>%
    dplyr::group_by(pathway) %>%
    dplyr::distinct(ID)
  
  return(total)
}

meta.map <- function(x){
  mSet<-InitDataObjects("conc", "msetora", FALSE)
  cmpd.vec<-x
  mSet<-Setup.MapData(mSet, cmpd.vec)
  mSet<-CrossReferencing(mSet, "name")
  mSet<-CreateMappingResultTable(mSet)
  res <- mSet$dataSet$map.table
  rm(mSet)
  return(res)
}


kegg_meta <- function(input){
  x <- clusterProfiler::enricher(gene = input, TERM2GENE = total, minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 1)
  df <- x@result[-1] %>%
    tidyr::separate(`GeneRatio`,into=c("A","B"),sep="/") %>%
    dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>%
    dplyr::mutate(GeneRatio=A/B) %>%
    dplyr::arrange(pvalue)
  #' 2.定义因子
  df <- subset(df, pvalue < 0.05)
  df$Description <- factor(df$Description,levels = c(df$Description %>% as.data.frame() %>% dplyr::pull()))
  return(df)
}


