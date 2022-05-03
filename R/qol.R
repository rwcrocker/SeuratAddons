### Quality of Life functions for Seurat ###

#' Generate ggplot2 themes for common Seurat plots
#'
#' This function provides custom themes for several
#' common Seurat functions, making them presentation-
#' ready.
#'
#' @import ggplot2
#' @param plot_type One of "base", "dot", "dim", "vln", or "comp"
#' @return ggplot2 theme
#' @export

theme_seurat = function(plot_type="base"){

  family = "sans"
  face = "bold"

  plot_themes = list()
  plot_themes[["base"]] = theme(axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.text = element_text(family = family, face = face),
                               legend.text = element_text(family = family, face = face),
                               legend.key.size = unit(1, 'cm'))
  plot_themes[["dot"]] = plot_themes[["base"]]+
    theme(panel.grid.major = element_line(color = "#e0e0e0"),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = face)
          )
  plot_themes[["dim"]] = plot_themes[["base"]]+
    theme(legend.text=element_text(size=20),
          axis.text=element_text(size=18),
          axis.title=element_text(size=20),
          plot.title = element_text(size=25, hjust = 0.5)
    )
  plot_themes[["vln"]] = plot_themes[["base"]]+
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 20)
    )
  plot_themes[["comp"]] = plot_themes[["base"]]+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(family = family, face = face, size = 12),
          legend.title = element_blank()
    )

  plot_types = names(plot_themes)
  is_legit_type = plot_type %in% plot_types

  if(!is_legit_type){
    stop("Given plot_type doesn't match available themes!")
  }

  theme_out = plot_themes[[plot_type]]
  return(theme_out)
}

#' Add a Gene Module Score without numeric tail
#'
#' This function adds a gene module score to a Seurat object.
#' Built off the Seurat-native AddModuleScore() function. Unlike AddModuleScore(),
#' this function generates the new column with the exact string given.
#'
#' @import Seurat
#' @param obj Seurat object to be scored
#' @param geneset Vector of gene symbols
#' @return Seurat object with calculated module score as meta.data column
#' @export

add_Module = function(obj, genelist){
  labels = names(genelist)
  labels_preexist = labels %in% colnames(obj@meta.data)
  if(any(labels_preexist)){
    warning(paste0("Overwriting existing meta.data label(s): ",
                   paste(labels[labels_preexist], collapse = ", ")))
  }

  holder = AddModuleScore(obj,
                          features = genelist,
                          name = labels)

  for (label in labels){
    index = which(labels == label)
    colnames(obj@meta.data) = gsub(colnames(obj@meta.data),
                                   pattern = paste0(label, index),
                                   replacement = label
    )
  }

  labels_exist = labels %in% colnames(obj@meta.data)
  if(!all(labels_exist)){stop("Label could not be generated in the object@meta.data...")}
  else
    return(obj)
}


#' Get genes from a downloadable Gene Ontology text file (MGI)
#'
#' This function reads a MGI formated GO file and extracts the genes that
#' make up the gene ontology term.
#'
#' @import readr tools
#' @param path File or directory path
#' @return List of gene vectors
#' @export

read_GO = function(path){

  if(!file.exists(path)){stop("File path does not exist.")}

  pathtype = (length(dir(path)) > 1) + 1
  pathtype = switch(pathtype, "file", "directory")
  print(paste0("Given path is a: ", pathtype))

  files = NA
  if (pathtype == "directory"){
    files = paste0(path, dir(path))
  }
  else{
    files = c(path)
  }

  genename_synonyms = c("Symbol", "symbol", "SYMBOL", "Gene_Symbol", "Gene_symbol", "gene_symbol", "GENE_SYMBOL")
  GO_list = list()
  for(file in files){
    print(paste0("Reading file: ", paste(file, collapse = "\n")))
    df = read_tsv(file)
    genename_column = genename_synonyms[which(genename_synonyms %in% colnames(df))]
    genes = as.vector(df[[genename_column]])
    GO_list[[file_path_sans_ext(basename(file))]] = genes
  }
  return(GO_list)
}

###  GSEA on Seurat Marker list  ###

#  use signed -log p value to rank
# run_GSEA = function(total_markers, genelist, use_signed_pvalue=TRUE){
#   total_markers = total_markers %>%
#     mutate(sign = sign(avg_log2FC)) %>%
#     mutate(signed_neg_log_padj = avg_log2FC * (-log10(p_val_adj)))%>%
#     arrange(desc(signed_neg_log_padj))
#
#   total_markers$signed_neg_log_padj[which(is.infinite(total_markers$signed_neg_log_padj))] = 10000000 # this is jank
#
#   ranks = as.vector(total_markers$signed_neg_log_padj)
#   names(ranks) = rownames(total_markers)
#   output = list()
#   output[["table"]]= fgsea(stats = ranks,
#                            pathways = genelist)
#   output[["plots"]] = list()
#   for (num in seq_along(genelist)){
#     output[["plots"]][[names(genelist)[num]]] = plotEnrichment(stats = ranks, pathway = genelist[[num]])
#   }
#   return(output)
# }

#' Get percent nonzero expressing cells in Seurat obj
#'
#' This function fetches expression data for a given geneset and
#' calculates the percentage of cells in the whole object with non-zero
#' expression levels for each gene.
#'
#' @import Seurat
#' @param obj Seurat object
#' @param geneset Vector of genes

percent_Nonzero = function(object, genes){
  genes_legit = genes %in% rownames(object)

  if(any(!genes_legit)){
    warning(
      paste0("Some genes not found in object: ",
             paste(genes[!genes_legit], collapse = ", "))
    )
  }

  genes = genes[genes_legit]
  counts = FetchData(object,
                     vars = genes,
                     slot = "counts")
  percents = as.vector(colMeans(counts > 0))*100
  names(percents) = genes
  return(percents)
}

#' Find correlations between a query and a geneset
#'
#' This function calculates the Pearson correlation between a query gene and
#' all detected genes in the Seurat object, returning a dataframe.
#' Specify a gene subset to improve efficiency.
#'
#' @import dplyr
#' @param obj Seurat object
#' @param query Vector of genes to assess
#' @param subset Subset of total genes in obj to correlate against
#' @export

find_Correlations = function(obj, query, subset=NULL){
  genes = rownames(obj@assays$RNA@data)
  if(!(is.null(subset))){
    genes = genes[which(genes %in% subset)]
  }

  cors = sapply(genes, FUN = function(gene){
    cor(obj@assays$RNA@data[query,], obj@assays$RNA@data[gene,])
  })

  df = data.frame(r = unlist(cors))
  df = df %>%
    mutate(gene = rownames(df))%>%
    arrange(r)
  return(df)
}


