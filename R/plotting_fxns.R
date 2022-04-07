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


