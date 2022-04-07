# General Plotting Fxns

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

#' Plot percent of cells with any detectable expression
#' @import ggplot2
#' @param obj Seurat object
#' @param gene Single gene name
#' @param x.margin X margin size
#' @param y.margin Y margin size
#' @param title.margin Title margin from rectangle
#' @param title.size Title text size
#' @param color Color of percent fill
#' @export
#'

plot_Nonzero = function(obj, gene, x.margin=5, y.margin=0, title.margin = 5, title.size = 8, color = "blue"){
  df = data.frame(
    pct = percent_Nonzero(obj, genes = gene),
    gene = gene
  )
  df$label = paste0(round(df$pct), "%")

  plot = ggplot(df)+
    geom_rect(aes(xmin=-1, xmax=1, ymin=0, ymax=100), color="black", size = 1.5, fill="light grey")+
    geom_rect(aes(xmin=-1, xmax=1, ymin=0, ymax=pct), color="black", size = 1.5, fill=color)+
    geom_text(aes(x=0, y=pct/2, label = label), color = "white", size = 12)+
    geom_text(aes(x=0, y=100+title.margin, label = "Percent\nNon-Zero"), color = "black", size = title.size, vjust=0, fontface="bold")+
    ylim(0-y.margin,100+y.margin+title.margin)+
    xlim(-1*x.margin, x.margin)+
    theme_void()+theme(text = element_text(face = "bold"))
  return(plot)
}

