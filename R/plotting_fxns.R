##_________________________________________###
###   Generate Cell Distribution Metrics    ###

#' Plot differences in cell composition
#'
#' This function plots the distribution of cell populations in a Seurat object,
#' grouping by a meta.data factor column.
#'
#' @import dplyr ggplot2
#' @param obj Seurat object
#' @param group.by Meta.data column name to annotate
#' @param split.by Meta.data column name to split plots
#' @param format Format of output plot. Acceptable formats: "pie", "bar".
#' @param group.levels Vector of levels to order group.by with
#' @param split.levels Vector of levels to order split.by with
#' @export

plot_composition = function(obj, group.by, split.by, format = "pie", group.levels=NULL, split.levels=NULL){
  df = obj@meta.data
  vars_real = all(c(group.by, split.by) %in% colnames(df))
  if(!vars_real) stop("Group.by or split.by are not a meta.data column.")
  vars_legit = every(df[c(group.by, split.by)], ~is.character(.x)|is.factor(.x))
  if(!vars_legit) stop("Group.by or split.by are not character/vector columns.")

  calculate_prop = function(vec){prop.table(table(vec))}

  usplits = unique(df[[split.by]])
  ugroups = unique(df[[group.by]])

  props = map(usplits, function(usplit){
    sub_df = df %>%
      filter(get(split.by) == usplit)
    sub_prop = calculate_prop(sub_df[[group.by]])
    return(sub_prop)
  })
  props = unlist(props)

  prop_df = data.frame(
    split = rep(usplits, each=length(unique(ugroups))),
    group = names(props),
    prop = props
  )

  if(!is.null(split.levels)) prop_df$split = factor(prop_df$split, levels = split.levels)
  if(!is.null(group.levels)) prop_df$group = factor(prop_df$group, levels = group.levels)

  plot = ggplot(prop_df)+
    geom_col(aes(x="", y=prop, fill = group))+
    facet_grid(~split)+
    theme_seurat("comp")
  if(format == "pie") plot = plot + coord_polar(theta = "y")

  return(plot)
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


