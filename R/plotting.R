# General Plotting Fxns

#' Plot differences in cell composition
#'
#' @import dplyr ggplot2
#' @param obj Seurat object
#' @param group.by Meta.data column name to annotate
#' @param split.by Meta.data column name to split plots
#' @param as.pie Output plot as pie chart, defaults to FALSE
#' @param group.levels Vector of levels to order group.by with
#' @param split.levels Vector of levels to order split.by with
#' @param label Labels to be plotted on plot instead of in legend, defaults to FALSE
#' @param label.pct Label to be plotted with proportion, defaults to FALSE
#' @export

plot_composition = function(obj, group.by, split.by=NULL, as.pie = FALSE, group.levels=NULL, split.levels=NULL, label = FALSE, label.pct = FALSE){
  df = obj@meta.data
  calculate_prop = function(vec){
    t = prop.table(table(vec))
    out = as.vector(t)
    names(out) = names(t)
    return(out)
  }

  # Check for erroneous inputs
  vars_real = all(c(group.by, split.by) %in% colnames(df))
  if(!vars_real) stop("Group.by or split.by are not a meta.data column.")
  vars_legit = every(df[c(group.by, split.by)], ~ is.character(.x) | is.factor(.x))
  if(!vars_legit) stop("Group.by or split.by are not character/vector columns.")

  # Dynamically set groups/splits
  groups = unique(df[[group.by]])

  splits = NULL
  should_split = !is.null(split.by)
  if(should_split){
    splits = unique(df[[split.by]])
    props = map(splits,
                function(split){
                  sub_df = df %>%
                    filter(get(split.by) == split)
                  sub_prop = calculate_prop(sub_df[[group.by]])
                  return(sub_prop)
                }
    )
    props = unlist(props)
  }

  else{
    props = calculate_prop(df[[group.by]])
  }

  # Generate proportion dataframe to plot
  prop_df = data.frame(
    group = names(props),
    prop = props
  )
  prop_df$label = paste0(round(prop_df$prop*100, 1), "% ", prop_df$group)

  if(!is.null(splits)){
    prop_df$split = rep(splits, each=length(unique(groups)))
  }

  # Organize levels by proportion (big on top)
  prop_df = prop_df %>%
    arrange(desc(prop))
  prop_df$group = factor(prop_df$group, levels = unique(prop_df$group))

  if(!is.null(split.levels)) prop_df$split = factor(prop_df$split, levels = split.levels)
  if(!is.null(group.levels)) prop_df$group = factor(prop_df$group, levels = group.levels)

  # Plotting with conditionals for provided split.by, format, and label
  plot = ggplot(prop_df, aes(x = "", y = prop, fill = group, label = label))+
    geom_col()+
    theme_seurat("comp")
  if(!is.null(splits)){plot = plot+facet_grid(~split)}
  if(as.pie){plot = plot+coord_polar(theta = "y")}
  if(label){
    if(label.pct){
      plot = plot+
        geom_text(aes(label = label), position = position_stack(vjust = 0.5))+
        theme(legend.position = "none")
    }
    else{
      plot = plot+
        geom_text(aes(label = group), position = position_stack(vjust = 0.5))+
        theme(legend.position = "none")
    }
  }

  return(plot)
}

#' Plot percent of cells with any detectable expression
#' @import ggplot2 Seurat
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
    xlim(( x.margin * -1 ), x.margin)+
    theme_void()+
    theme(text = element_text(face = "bold"))
  return(plot)
}

