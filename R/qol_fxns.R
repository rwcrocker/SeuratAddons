### Quality of Life functions for Seurat ###

require(ggplot2)
require(dplyr)
require(cowplot)
require(Seurat)
require(fgsea)
require(scales)


###______________________________________###
###   A more convenient AddModuleScore   ###

#IN:      Seurat Object       (Must be self-replacing, i.e.   obj <- fxn(obj)  )
#         Vector of module component genes
#         Name for your module
#OUT:     Seurat object w/ new module score in meta.data


add_Module = function(obj, geneset, label){

  label_exists = label %in% colnames(obj@meta.data)
  if(label_exists){print(paste0("Overwriting existing label: ", label))}

  obj = AddModuleScore(obj,
                       features = list(geneset),
                       name = label)
  auto_label = paste0(label, 1)
  obj@meta.data[[label]]= obj@meta.data[[auto_label]]
  removal = names(obj@meta.data) == auto_label
  obj@meta.data = obj@meta.data[,!removal]

  label_correct = label %in% colnames(obj@meta.data)
  if(!label_correct){stop("Label could not be generated in the object@meta.data...")}

  return(obj)
}

### Stacked Vln  ###
stacked_vln <- function(obj, feats, clusters, alt.condition.names=c()){
  #adjust x axis names if desired
  sample_names <- levels(obj@meta.data$Sample)
  if (length(alt.condition.names) != 0){
    sample_names= alt.condition.names
    print("Alternative condition names given: ")
    print(alt.condition.names)
  }

  #Check if given clusters match obj clusters
  if (all(clusters %in% obj@active.ident) == FALSE){
    stop("Given clusters do not match object clusters!")
  }

  #dummy df
  ddf <- data.frame(x=c(0,1), y=c(0,1))
  #make plots
  stat_line <-stat_summary(fun.y = mean, geom='crossbar', size = .4, colour = "red")
  plots <- list()
  for (feature_number in 1:length(feats)){
    labelname <- paste0(feats[feature_number], "_Label")
    plots[[labelname]] <- ggplot(ddf)+annotate("text", family="Arial", x = 1, y = .5, size=8, label = feats[feature_number], hjust=1)+xlim(0,1)+ylim(0,1)+theme_void()
    for (cluster_number in 1:length(clusters)){
      plotname <- paste0(feats[feature_number], "_", clusters[cluster_number])
      plots[[plotname]] <- VlnPlot(obj, features = feats[feature_number], group.by = "sample", idents = c(clusters[cluster_number]), pt.size = 0.01)+
        stat_line+theme_bare+
        theme(plot.title = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Arial"), legend.position = 'none')+
        scale_x_discrete(labels=sample_names)
    }
  }

  cluster_labels <- list()
  cluster_labels[["Empty"]] <- ggplot(ddf)+theme_void()
  for (num in 1:length(clusters)){
    cluster_labels[[clusters[num]]] <-  ggplot(ddf)+annotate("text", x = .5, y = .25, size=8, label = clusters[num], hjust=0.5)+xlim(0,1)+ylim(0,1)+theme_void()
  }

  plots <- c(cluster_labels, plots)
  print(names(plots))
  plot_grid(plotlist = plots, ncol = length(clusters)+1)
}

###__________________###
###  ggplot2 themes  ###

ppg_family <- "Arial"
ppg_face <- "bold"

theme_dotplot <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(family = "Arial", face = "bold"),
                       panel.grid.major = element_line(color = "#e0e0e0"))

umap_theme <- theme(legend.text=element_text(size=20, family = ppg_family),
                    axis.text=element_text(size=18, family = ppg_family),
                    axis.title=element_text(size=20, family = ppg_family),
                    plot.title = element_text(size=25, hjust = 0.5, family = ppg_family))
theme_bare <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(family = ppg_family, face = ppg_face))


big_color <- c(hue_pal()(10), "navy", "black", "grey", "red")

###  Get genes from Go term downloads  ###

read_GO = function(path){
  require(readr)
  genename_synonyms = c("Symbol", "symbol", "SYMBOL", "Gene_Symbol", "Gene_symbol", "gene_symbol", "GENE_SYMBOL")


  files = dir(path)
  if(length(files) != 0){
    files = paste0(path, files)
    }
  else{
    files = c(path)
  }
  files_exist = file.exists(files)
  if(!all(files_exist)){
    stop(paste0("One or more files do not exist among:", paste(files,collapse = "\t")))
    }

  GO_list = list()
  for(file in files){
    print(paste0("Reading file: ", paste(file, collapse = "\n")))
    df = read_tsv(file)
    genename_column = genename_synonyms[which(genename_synonyms %in% colnames(df))]
    genes = as.vector(df[[genename_column]])
    GO_list[[file]] = genes
  }
  return(GO_list)
}


# From a directory of tsv formated go term file, read all and add to list
read_GO_directory = function(directory){
  print(directory)
  go_list = list()
  for (file in dir(directory)){
    name = tools::file_path_sans_ext(file)
    go_list[[name]] = read_GO(file = paste0(directory, file))
  }
  return(go_list)
}


###  GSEA on Seurat Marker list  ###

#  use signed -log p value to rank
run_GSEA = function(total_markers, genelist, use_signed_pvalue=TRUE){
  total_markers = total_markers %>%
    mutate(sign = sign(avg_log2FC)) %>%
    mutate(signed_neg_log_padj = avg_log2FC * (-log10(p_val_adj)))%>%
    arrange(desc(signed_neg_log_padj))

  total_markers$signed_neg_log_padj[which(is.infinite(total_markers$signed_neg_log_padj))] = 10000000 # this is jank

  ranks = as.vector(total_markers$signed_neg_log_padj)
  names(ranks) = rownames(total_markers)
  output = list()
  output[["table"]]= fgsea(stats = ranks,
                           pathways = genelist)
  output[["plots"]] = list()
  for (num in seq_along(genelist)){
    output[["plots"]][[names(genelist)[num]]] = plotEnrichment(stats = ranks, pathway = genelist[[num]])
  }
  return(output)
}

### Get nonzero percent  ###
get_Nonzero = function(obj, genes, select_cells=NULL){
  cells = colnames(obj)
  if(!(is.null(select_cells))){
    cells = selected_cells
  }
  detected_genes = genes[which(genes %in% rownames(obj))]
  data = FetchData(neuts, detected_genes, cells = cells)
  data = data.frame(percent_expressed = colMeans(data  > 0))*100
  return(data)
}



# Find gene correlations
find_Correlations = function(obj, querry, gene_subset=NULL){
  genes = c()
  if(!(is.null(gene_subset))){
    genes = gene_subset
  }
  else{ genes = rownames(obj@assays$RNA@data)}

  correlations = lapply(genes, FUN = function(gene){
    cor(obj@assays$RNA@data[querry,], obj@assays$RNA@data[gene,])
  })
  names(correlations)=genes
  return(correlations)
}

#
add_Annotation = function(obj, col_label, convert_vector){
  obj@meta.data[[col_label]] = NA
  for (cluster_number in levels(obj@meta.data$seurat_clusters)){
    index_num = as.numeric(cluster_number)+1
    obj@meta.data[[col_label]][which(str_detect(obj@meta.data$seurat_clusters, cluster_number))] = convert_vector[index_num]
  }
  return(obj)
}



