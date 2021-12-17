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

###_________________________###
###   Plot QC Metrics   ###

#IN:      Seurat Obj
#Options: Cut-offs, alternate column names, custom colors
#OUT:     Combined plot of QC metrics

plot_QC <- function(obj,
                    condition_name = "sample", condition_colors = NULL, condition_levels=NULL,
                    expected_cells = NULL,
                    nUMI_cut=500, nUMI_name = "nUMI",
                    nGene_cut=250, nGene_name = "nGene",
                    complexity_cut=0.8, complexity_name = "log10GenesPerUMI",
                    mitoRatio_cut=0.2, mitoRatio_name = "mitoRatio",
                    axis_label_adjusts=c(15, 20, 28, 16)
){
  meta = obj@meta.data

  # Checking for condition name and corresponding column
  if(!(condition_name %in% colnames(meta))){
    stop("Provided condition name must be a column in the obj meta.data")
  }
  condition_column = meta[[paste(condition_name)]]

  # Setting and checking colors
  num_conditions = length(unique(condition_column))
  if (is.null(condition_colors)){
    condition_colors = hue_pal()(num_conditions)
  }

  ###  Plotting  ###

  # Plot theme
  theme_qc = theme(text = element_text(face = "bold", colour = "black", size = 20),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   plot.title = element_text(hjust=0.5),
                   axis.title.y = element_text(size = 15),
                   axis.title.x = element_text(size = 15),
                   legend.position = 'none')

  plots = list()

  # nCell
  if (!(is.null(condition_name))){
    plots[["Total Cell Number"]] = ggplot(meta, aes(x=condition_column, fill=condition_column))+
      geom_bar()+
      geom_hline(yintercept = expected_cells, linetype="dashed")+
      theme_classic()+
      ggtitle("Cell Count")+
      theme_qc+
      theme(axis.title.x = element_blank())+
      scale_fill_manual(values = condition_colors)
  }

  #nUMI
  if (!(is.null(nUMI_name))){
    nUMI_column = meta[[paste(nUMI_name)]]
    plots[["nUMI"]] = ggplot(meta, aes(x=nUMI_column, fill=condition_column, color=condition_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = nUMI_cut, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("nUMI")+
      theme_qc+
      scale_fill_manual(values = condition_colors)+
      scale_color_manual(values = condition_colors)+
      xlab("UMI/Cell")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[1])))
  }

  if (!(is.null(nGene_name))){
    nGene_column = meta[[paste(nGene_name)]]
    plots[["nGene"]] = ggplot(meta, aes(x=nGene_column, fill=condition_column, color=condition_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = nGene_cut, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("nGene")+
      theme_qc+
      scale_fill_manual(values = condition_colors)+
      scale_color_manual(values = condition_colors)+
      xlab("Genes/Cell")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[2])))
  }

  if (!(is.null(complexity_name))){
    complexity_column = meta[[paste(complexity_name)]]
    plots[["Complexity"]] = ggplot(meta, aes(x=complexity_column, fill=condition_column, color=condition_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = complexity_cut, linetype="dashed")+
      #scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("Txn Complexity")+
      theme_qc+
      scale_fill_manual(values = condition_colors)+
      scale_color_manual(values = condition_colors)+
      xlab(label = "log10(nGene/nUMI)")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[3])))
  }

  if (!(is.null(mitoRatio_name))){
    mitoRatio_column = meta[[paste(mitoRatio_name)]]
    plots[["mitoRatio"]] = ggplot(meta, aes(x=mitoRatio_column, fill=condition_column, color=condition_column))+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = mitoRatio_cut, linetype="dashed")+
      scale_x_log10()+
      theme_classic()+
      ylab("Cell density")+
      ggtitle("Mitochondrial Ratio")+
      theme_qc+
      scale_fill_manual(values = condition_colors)+
      scale_color_manual(values = condition_colors)+
      xlab(label = "mitoUMI / nUMI")+
      theme(axis.text.x = element_text(margin=margin(b=axis_label_adjusts[4])))
  }

  # add legend to last plot
  plots[["Legend"]] = get_legend(plots[[1]]+theme(legend.position = "right", legend.title = element_blank(), legend.key.size = unit(1.5, units = "cm")))

  #plot
  plot_grid(plotlist = plots, nrow = 1, axis = "b", align = "h")
}


###_________________________________________###
###   Generate Cell Distribution Metrics    ###

#IN:      Seurat obj, meta.data column for cell annotation and for sample id (condition)
#OPTIONS: List of vectors length(2) comparing conditions, releveling options, and whether to assign the raw dfs to temp_
#OUT:     Combined plot of cell composition of seurat dataset

plot_composition = function(obj,
                            annotation_name, condition_name,
                            comparison_list=NULL,
                            condition_levels = NULL, annotation_levels=NULL,
                            temp_out=FALSE, show_legend=TRUE){

  conditions = unique(obj@meta.data[[condition_name]])
  annotations = unique(obj@meta.data[[annotation_name]])

  # get ncells from seurat obj
  get_ncells = function(obj, annotation, condition){
    meta = as.data.frame(obj@meta.data) %>%
      filter(get(annotation_name) == annotation)%>%
      filter(get(condition_name) == condition)

    ncells = length(meta[,1])
    return(ncells)
  }

  # construct dataframe for ggplotting
  df_condition = unlist(lapply(conditions, function(x){rep(x, length(annotations))}))
  df_annotation = rep(annotations, length(conditions))
  df = data.frame(df_condition, df_annotation)
  df$df_ncell = 0
  for (num in seq_along(df$df_condition)){
    df$df_ncell[num]=get_ncells(obj, annotation = df$df_annotation[num], condition = df$df_condition[num])
  }
  df$df_condition_total = 0
  for (cond in conditions){
    df$df_condition_total[which(df$df_condition == cond)]= sum(df$df_ncell[which(df$df_condition == cond)])
  }

  df = df %>%
    mutate(df_percent = (df_ncell / df_condition_total))

  # Construct df of percent change per condition
  change_df = data.frame(annotation = annotations)
  for (cond in conditions){
    change_df[[cond]]=NA
  }
  for (anno in change_df$annotation){
    for(cond in conditions){
      change_df[[cond]][which(change_df$annotation == anno)] = df$df_percent[which(df$df_annotation == anno & df$df_condition == cond)]
    }
  }

  # Confirm all comparisons provided are plottable
  comparisons_provided = is.list(comparison_list) & length(comparison_list) > 0

  is_comparison_legit = function(comparison, table){
    is_vector = is.vector(comparison)
    is_two = length(comparison) == 2
    is_present = all(comparison %in% colnames(table))
    is_legit = all(is_vector, is_two, is_present)
    return(is_legit)
  }

  num_comparisons = 0
  if(comparisons_provided){
    num_comparisons = length(comparison_list) # to use later for scaling plot widths
    comparisons_legit = sapply(comparison_list, FUN = function(x){
      is_comparison_legit(x, table = change_df)
    })
    all_comparisons_legit = all(comparisons_legit)
    if (!all_comparisons_legit){
      stop("Provided comparison_list is broken")
    }
  }

  ###  RELEVELING  ###
  # TODO make this prettier and add a check to make sure all levels provided are present

  # is_leveling_legit = function(levels, table, ){
  #   is_vector = is.vector(levels)
  #   is_present = all(levels %in% colnames(table))
  #   is_legit = all(is_vector, is_two, is_present)
  #   return(is_legit)
  # }

  if(is.vector(condition_levels) & length(condition_levels) == length(conditions)){
    df$df_condition = factor(df$df_condition, levels = condition_levels)
    conditions = condition_levels
  }

  if(is.vector(annotation_levels) & length(annotation_levels) == length(annotations)){
    df$df_annotation = factor(df$df_annotation, levels = annotation_levels)
  }

  ###  PLOTTING  ###
  annotation_colors = hue_pal()(length(annotations))
  final_plots = list()

  # Ncell Bar
  final_plots[["ncell_plot"]] <- ggplot(data = df, aes(x=df_condition, y=df_ncell, fill=df_annotation)) +
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
    scale_fill_manual(values=c(annotation_colors))

  # Pies
  subpies = list()
  for (cond in conditions){
    sub_df <- df %>% filter(df_condition == cond)
    subpies[[cond]] <- ggplot(data = sub_df, aes(x="", y=df_percent, fill=df_annotation))+
      geom_col(stat="identity", width = 1)+
      ggtitle(cond)+
      coord_polar("y", start = 0)+
      theme_void()+
      theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            legend.position = 'none')+
      scale_fill_manual(values=c(annotation_colors))
  }
  final_plots[["pie_plots"]] = plot_grid(plotlist = subpies, ncol = 1)

  # Percent change bar
  if(comparisons_provided){
    comparison_plots = list()
    for (num in seq_along(comparison_list)){
      comparison_name = names(comparison_list)[num]
      comparison = comparison_list[[num]]
      change_df[[paste0(comparison_name)]]= ( change_df[[paste0(comparison[2])]] - change_df[[paste0(comparison[1])]] ) * 100
      xmax = max(abs(change_df[[paste0(comparison_name)]]))

      comparison_plots[[comparison_name]] = ggplot(change_df,
                                                   aes(x=.data[[paste0(comparison_name)]],
                                                       y = annotation,
                                                       fill = annotation))+
        geom_col()+
        geom_vline(xintercept = 0)+
        ggtitle(comparison_name)+
        theme_bw()+
        theme(legend.position = "none",
              axis.title = element_blank(),
              text = element_text(face = "bold", size = 16, colour = "black"),
              plot.title = element_text(size = 18))+
        xlim(-xmax, xmax)+
        scale_y_discrete(limits = unique(change_df$annotation)[rev(order(unique(change_df$annotation)))])
    }
    final_plots[["comparison_plots"]] = plot_grid(plotlist = comparison_plots, nrow = 1)
  }


  # Add legend depending on option
  if(show_legend){
    plot_to_steal_legend = final_plots[[1]]
    final_plots[["legend"]] = get_legend(plot_to_steal_legend+
                                           theme(legend.position = "right",
                                                 legend.title = element_blank(),
                                                 legend.key.size = unit(1.5, units = "cm")))
  }

  # Final merged plot
  final = plot_grid(plotlist = final_plots, nrow = 1, rel_widths = c(2, 1, (1 + num_comparisons), 2))
  return(final)

  # Make raw data frames available
  if (temp_out){
    assign("temp_composition_df", df, envir = .GlobalEnv)
    assign("temp_change_df", change_df, .GlobalEnv)
  }
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


