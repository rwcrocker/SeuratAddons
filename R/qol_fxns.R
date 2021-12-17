### Quality of Life functions for Seurat ###

# ggplot themes for consistency
FAMILY = "Arial"
FACE = "bold"

theme_dotplot = theme_bare+
  theme(panel.grid.major = element_line(color = "#e0e0e0"))

theme_umap <- theme(legend.text=element_text(size=20, family = FAMILY),
                    axis.text=element_text(size=18, family = FAMILY),
                    axis.title=element_text(size=20, family = FAMILY),
                    plot.title = element_text(size=25, hjust = 0.5, family = FAMILY))

theme_bare <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                    axis.text = element_text(family = FAMILY, face = FACE))


#' Add a Gene Module Score
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

add_Module = function(obj, genelist, labels){

  labels_preexist = labels %in% colnames(obj@meta.data)
  preexisting = labels[labels_preexist]
  if(any(labels_preexist)){
    print(paste0("Overwriting existing label(s): ", paste(preexisting, collapse = ", ")))
  }

  obj = AddModuleScore(obj,
                       features = genelist,
                       name = labels)

  for (num in seq_along(labels)){
    label = labels[num]
    generated_label = paste0(label, num)
    obj@meta.data[[label]] = obj@meta.data[[generated_label]]
    generated_column_number = which(colnames(obj@meta.data) == generated_label)
    col2remove = colnames(obj@meta.data) == generated_label
    obj@meta.data = obj@meta.data[,!col2remove]
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
#' @param path file or directory path
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
#' @import
#' @param obj Seurat object
#' @param geneset Vector of genes

percent_Nonzero = function(obj, geneset, select_cells=NULL){
  cells = colnames(obj)
  if(!(is.null(select_cells))){
    cells = selected_cells
  }
  detected_genes = geneset[which(geneset %in% rownames(obj))]
  data = FetchData(neuts, detected_genes, cells = cells)
  data = data.frame(percent_expressed = colMeans(data  > 0))*100
  return(data)
}

#' Find correlations between a query and a geneset
#'
#' This function calculates the Pearson correlation between a query gene and
#' all detected genes in the Seurat object, returning a dataframe.
#' Specify a gene subset to improve efficiency.
#'
#' @import dplyr
#' @param obj
#' @param query
#' @param subset
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

#' Add meta-data annotation
#'
#' This function uses cluster meta-data and a matching conversion vector to add
#' a new annotation to a Seurat object's meta-data.
#'
#' @import stringr
#' @param obj Seurat object
#' @param col_label Name for new meta.data column
#' @param convert Vector with values corresoning to reference index
#' @param reference Meta.data column-name to use for cluster identity
#' @export

add_Annotation = function(obj, col_label, convert, reference = "seurat_clusters"){
  obj@meta.data[[col_label]] = NA
  levels_match = length(convert) == length(unique(obj@meta.data[[reference]]))
  if(!levels_match){stop("Unequal number of cluster and conversion levels.")}

  for (cluster_number in levels(obj@meta.data[[reference]])){
    index_num = as.numeric(cluster_number)+1
    obj@meta.data[[col_label]][which(str_detect(obj@meta.data[[reference]], cluster_number))] = convert[index_num]
  }
  return(obj)
}
