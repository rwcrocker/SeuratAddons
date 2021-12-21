# Plotting Functions


# ggplot themes
FAMILY = "Arial"
FACE = "bold"
theme_bare <- theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                    axis.text = element_text(family = FAMILY, face = FACE))
theme_dotplot <- theme_bare+
  theme(panel.grid.major = element_line(color = "#e0e0e0"))
umap_theme <- theme_bare +
  theme(legend.text=element_text(size=20, family = FAMILY),
        plot.title = element_text(size=25, hjust = 0.5, family = FAMILY))


#' Plot QC metrics for Seurat Object
#'
#' This function plots common QC metrics (nUMI, nGene,nGene:nUMI, and mitoRatio)
#'  for a Seurat object with corresponding meta.data columns.
#'
#' @import scales ggplot2 cowplot
#' @param obj Seurat object
#' @param annotation
#' @param annotation_colors
#' @param annotation_levels
#' @param expected_cells
#' @param nUMI_cuofft
#' @param nUMI_name
#' @param nGene_cutoff
#' @param nGene_name
#' @param complexity_cutoff
#' @param complexity_name
#' @param mitoRatio_cutoff
#' @param mitoRatio_name
#' @param axis_label_adjusts Length 4 vector for axis title adjustment
#' @export

plot_QC <- function(obj,
                    annotation = "sample", annotation_colors = NULL, annotation_levels=NULL,
                    expected_cells = NULL,
                    nUMI_cutoff=500, nUMI_name = "nUMI",
                    nGene_cutoff=250, nGene_name = "nGene",
                    complexity_cutoff=0.8, complexity_name = "log10GenesPerUMI",
                    mitoRatio_cutoff=0.2, mitoRatio_name = "mitoRatio",
                    axis_label_adjusts=c(15, 20, 22, 16)
){
  meta = obj@meta.data

  # Checking for condition name and corresponding column
  if(!(annotation %in% colnames(meta))){
    stop("Provided condition name must be a column in the object meta.data")
  }

  #apply levels to factor if given
  if(!is.null(annotation_levels) & length(annotation_levels) == length(unique(meta[[annotation]]))){
    meta[[annotation]] = factor(meta[[annotation]], levels = annotation_levels)
  }


  annotation_column = meta[[paste(annotation)]]


  # Setting and checking colors
  num_annotations = length(unique(annotation_column))
  if (is.null(annotation_colors)){
    annotation_colors = hue_pal()(num_annotations)
  }
  # Plot theme
  theme_qc = theme(text = element_text(face = "bold", colour = "black", size = 20),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   plot.title = element_text(hjust=0.5),
                   axis.title.y = element_text(size = 15),
                   axis.title.x = element_text(size = 15),
                   legend.position = 'none')
  # Plotting
  plots = list()
  # nCell
  if (!(is.null(annotation))){
    plots[["Total Cell Number"]] = ggplot(meta, aes(x=annotation_column, fill=annotation_column))+
      geom_bar()+
      geom_hline(yintercept = expected_cells, linetype="dashed")+
      theme_classic()+
      ggtitle("Cell Count")+
      theme_qc+
      theme(axis.title.x = element_blank())+
      scale_fill_manual(values = annotation_colors)
  }

  #nUMI
  if (!(is.null(nUMI_name))){
    nUMI_column = meta[[paste(nUMI_name)]]
    plots[["nUMI"]] = ggplot(meta, aes(x=nUMI_column, fill=annotation_column, color=annotation_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = nUMI_cutoff, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("nUMI")+
      theme_qc+
      scale_fill_manual(values = annotation_colors)+
      scale_color_manual(values = annotation_colors)+
      xlab("UMI/Cell")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[1])))
  }
  #nGene
  if (!(is.null(nGene_name))){
    nGene_column = meta[[paste(nGene_name)]]
    plots[["nGene"]] = ggplot(meta, aes(x=nGene_column, fill=annotation_column, color=annotation_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = nGene_cutoff, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("nGene")+
      theme_qc+
      scale_fill_manual(values = annotation_colors)+
      scale_color_manual(values = annotation_colors)+
      xlab("Genes/Cell")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[2])))
  }
  #Complexity
  if (!(is.null(complexity_name))){
    complexity_column = meta[[paste(complexity_name)]]
    plots[["Complexity"]] = ggplot(meta, aes(x=complexity_column, fill=annotation_column, color=annotation_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = complexity_cutoff, linetype="dashed")+
      #scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("Txn Complexity")+
      theme_qc+
      scale_fill_manual(values = annotation_colors)+
      scale_color_manual(values = annotation_colors)+
      xlab(label = "log10(nGene/nUMI)")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[3])))
  }
  #mitoRatio
  if (!(is.null(mitoRatio_name))){
    mitoRatio_column = meta[[paste(mitoRatio_name)]]
    plots[["mitoRatio"]] = ggplot(meta, aes(x=mitoRatio_column, fill=annotation_column, color=annotation_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = mitoRatio_cutoff, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("Mitochondrial Ratio")+
      theme_qc+
      scale_fill_manual(values = annotation_colors)+
      scale_color_manual(values = annotation_colors)+
      xlab(label = "mitoUMI / nUMI")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[4])))
  }
  # Add legend to last
  plots[["Legend"]] = get_legend(plots[[1]]+theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(1.5, units = "cm")))
  # Combine
  plot_grid(plotlist = plots, nrow = 1, axis = "b", align = "h")
}


##_________________________________________###
###   Generate Cell Distribution Metrics    ###

#IN:      Seurat obj, meta.data column for cell annotation and for sample id (condition)
#OPTIONS: List of vectors length(2) comparing conditions, releveling options, and whether to assign the raw dfs to temp_
#OUT:     Combined plot of cell composition of seurat dataset

#' Plot differences in cell composition
#'
#' This function plots the distribution of cell populations in a Seurat object,
#' grouping by a meta.data factor column.
#'
#' @import dplyr
#' @param obj Seurat object
#' @param group_by Column name of meta.data, i.e.cell-type
#' @param split_by Column name of meta.data, i.e condition
#' @param comparisons List of pairwise comparison vectors
#' @param comparison_names Vector of comparison names
#' @param split_levels Vector of levels for split_by
#' @param group_levels Vector of levels for group_by
#' @param show_legend Logical value to show legend
#' @export

plot_composition = function(obj,
                            group_by, split_by,
                            comparisons=NULL, comparison_names = NULL,
                            split_levels = NULL, group_levels=NULL,
                            temp_out=FALSE, show_legend=TRUE){

  splits = unique(obj@meta.data[[split_by]])
  groups = unique(obj@meta.data[[group_by]])

  # method for getting ncells from meta.data
  get_ncells = function(obj, group=NULL, split=NULL){
    meta = as.data.frame(obj@meta.data)
    if(!is.null(group)){
      meta = meta %>% filter(get(group_by) == group)
    }
    if(!is.null(split)){
      meta = meta %>% filter(get(split_by) == split)
    }
    ncells = dim(meta)[1]
    return(ncells)
  }

  # construct df of cell counts
  df = data.frame(
    df_splits = as.vector(sapply(splits, function(x){rep(x, length(groups))})),
    df_groups = rep(groups, length(splits))
    )

  df = df %>%
    mutate(df_ncell = apply(df, 1, FUN = function(x){
      get_ncells(obj, group = x[2], split = x[1])
    })) %>%
    mutate(df_totalsplit = apply(df, 1, FUN = function(x){
      get_ncells(obj, split = x[1])
    })) %>%
    mutate(df_percent = (df_ncell / df_totalsplit ))

  # Construct df of percent change per condition
  comparisons_given = is.list(comparisons)
  comparisons_legit = unlist(lapply(comparisons, function(vector){length(vector) == 2 & all(vector %in% splits)}))

  if (comparisons_given & all(comparisons_legit)){
    change_df = data.frame(change_groups = groups)
    for (comp in comparisons){
      reference = comp[1]
      experimental = comp[2]
      comp_name = paste0(reference, "_vs_", experimental)

      change_df[[comp_name]] = sapply(change_df$change_groups, FUN = function(change_group){
        exp = df$df_percent[which(df$df_splits == experimental & df$df_groups == change_group)]
        ref = df$df_percent[which(df$df_splits == reference & df$df_groups == change_group)]
        return(exp - ref)
      })
    }
  }

  # Refactor groups if provided
  if(!is.null(group_levels) & all(group_levels %in% groups)){
    df$df_groups = factor(df$df_groups, levels = group_levels)
    change_df$change_groups = factor(change_df$change_groups, levels = group_levels)
  }
  # Refactor splits if provided
  if(!is.null(split_levels) & all(split_levels %in% splits)){
    df$df_splits = factor(df$df_splits, levels = split_levels)
    #no splits present in change_df
  }

  ###  PLOTTING  ###
  group_colors = hue_pal()(length(groups))
  combined_plot = list()

  # Ncell Bar
  combined_plot[["ncell"]] <- ggplot(data = df, aes(x=df_splits, y=df_ncell, fill=df_groups)) +
    geom_col(position="dodge") +
    ggtitle("Cells Count")+
    ylab("Number of Cells\n")+
    xlab("")+
    theme_classic()+
    theme(text = element_text(size = 20, face = "bold", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust =1, vjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size=22, hjust=0.5, face="bold"),
          legend.position = "none",)+
    scale_fill_manual(values=c(group_colors))

  # Pies
  subpies = list()
  for (cond in splits){
    sub_df <- df %>% filter(df_splits == cond)
    subpies[[cond]] <- ggplot(data = sub_df,
                              aes(x="", y=df_percent, fill=df_groups))+
      geom_col(stat="identity", width = 1)+
      ggtitle(cond)+
      coord_polar("y", start = 0)+
      theme_void()+
      theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            legend.position = 'none')+
      scale_fill_manual(values=c(group_colors))
  }
  combined_plot[["pies"]] = plot_grid(plotlist = subpies, ncol = 1)

  # Percent change bar
  if(comparisons_given & all(comparisons_legit)){
    comparison_plots = lapply(comparisons, FUN = function(comp){

    comp_name = paste(comp, collapse = "_vs_")
    xmax = max(abs(change_df[[comp_name]])) *100

    given_comp_names = names(comparisons)

    ggplot(change_df, aes(x = .data[[comp_name]] * 100,
                            y = change_groups,
                          fill = change_groups))+
      geom_col()+
      geom_vline(xintercept = 0)+
      ggtitle(comp_name)+
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            text = element_text(face = "bold", size = 16, colour = "black"),
            plot.title = element_text(size = 18, hjust = 0.5))+
      xlim(-xmax, xmax)+
      scale_y_discrete(limits = unique(
        change_df$change_groups)[rev(order(unique(change_df$change_groups)))]
        )
    })

    # Change comparison titles if provided
    if(!is.null(comparison_names) & length(comparison_names) == length(comparisons)){
      for(num in seq_along(comparison_plots)){
        comparison_plots[[num]] = comparison_plots[[num]] + ggtitle(comparison_names[num])
      }
    }
    combined_plot[["comparison"]] = plot_grid(plotlist = comparison_plots, nrow = 1)
  }

  # Add legend depending on option
  if(show_legend){
    plot_to_steal_legend = combined_plot[[1]]
    combined_plot[["legend"]] = get_legend(plot_to_steal_legend+
                                           theme(legend.position = "right",
                                                 legend.title = element_blank(),
                                                 legend.key.size = unit(1.5, units = "cm")))
  }
  # Final combined plot
  combined_plot = plot_grid(plotlist = combined_plot, nrow = 1, rel_widths = c(2, 1, (1 + length(comparisons)), 2))
  return(combined_plot)
}


#ENRICHR PLOTTING
##take a DF and make a bar plot
enrichment_bar <- function(df, title=""){
  df$Term <- factor(df$Term, levels = rev(df$Term))
  df <- df[which(df$P.value <=0.05),]
  anno <- ""
  if (length(df$P.value) > 20){
    anno <- paste0( (length(df$P.value) - 20), " more\nwith padj < 0.05")
    df <- df[1:20,]
  }

  plot <- ggplot(df, aes(x=-log10(P.value), y=Term))+
    geom_col()+geom_vline(xintercept = -log10(0.05), color="red", size=1.25)+
    ggtitle(title)+
    annotate("text", x = (-log10(df$P.value[1])*.75), y=1, label=anno)+
    theme_bw()+
    theme(text = element_text(family = "Arial", face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 10),
          plot.title.position = "plot", plot.title = element_text(hjust = 0.5))
  return(plot)
}

#take a list of dfs and make plots for each
enrichment_bar_list <- function(df_list, title_prefix=""){
  plot_list <- list()
  for (df_num in 1:length(df_list)){
    df_name <- names(df_list)[df_num]
    plot_list[[df_name]] <- enrichment_bar(df_list[[df_name]], title = paste0(title_prefix, " ", df_name))
  }
  return(plot_list)
}

### VlnPlot with associated nonzero percent
vln_nonzero_plot = function(obj, gene, plot_ymax = NULL){


  # Assign ylim if given
  if (!(is.null(plot_ymax))){
    plot = VlnPlot(obj, features = gene, flip = T, y.max = plot_ymax)+
      theme(legend.position = "none")+theme_bare
  }
  # Grab pre-assigned ylim if not given
  else{
    plot = VlnPlot(obj, features = gene, flip = T)+
      theme(legend.position = "none")+theme_bare
    plot_ymax = ggplot_build(plot)$layout$panel_params[[1]]$y.range[2]
  }
  plot_ymax_buffered = plot_ymax * .95 #add buffer

  # get nonzero expression percentage
  nonzero = unlist(nonzero_percent(obj = obj, genes = gene))
  nonzero_decimal = nonzero / 100

  # Add nonzero plot
  plot = plot+
    geom_rect(aes(xmin=1.75, xmax=2.25, ymin=nonzero_decimal*plot_ymax_buffered, ymax=plot_ymax_buffered),
              color="black", fill="white")+
    geom_rect(aes(xmin=1.75, xmax=2.25, ymin=0, ymax=nonzero_decimal*plot_ymax_buffered),
              color="black", fill = "grey")+
    geom_text(aes(x=2, y=.5*plot_ymax_buffered,
                  label = paste0(round(nonzero, digits = 1), " %")))

  print(plot)
}


