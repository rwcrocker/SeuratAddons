### Functions for enrichR


#TODO add these to import
library(enrichR)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggprism)


major_dbs = c("WikiPathways_2019_Mouse", "MSigDB_Hallmark_2020",
                    "KEGG_2019_Mouse", "GO_Molecular_Function_2021",
                    "GO_Cellular_Component_2021", "GO_Biological_Process_2021")

#' Clean the dirty data output of enrichr()
#'

clean_enrichR = function(df){
  enrichr_vars = c("Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
  is_enrichr = all(enrichr_vars %in% colnames(df))
  if(!is_enrichr){stop("Data.frame doesn't look like an enrichR output (Column Names)")}

  df = df %>%
    arrange(Adjusted.P.value) %>%
    filter(Adjusted.P.value < 0.05)
  df$Term = factor(df$Term, levels = df$Term)

  # Correct gene symbols and reconcatonate
  df$Genes = df$Genes %>%
    str_split(pattern = ";") %>%
    lapply(FUN = function(i){
      i = stringr::str_to_title(i)
      i = paste(i, collapse = ",")
      return(i)
    })

  # Parese/Evaluate the fraction to a numeric
  df$Overlap = sapply(df$Overlap, function(i){eval(parse(text = i))})

  df = df %>%
    select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes)

  return(df)
}

#' Plot the output of the enrichR dataframe
plot_enrichr = function(df, ntop, fill_by = "Overlap", title = NULL, term_char=100){
  fill_by_legit = fill_by %in% colnames(df)
  if(!fill_by_legit){stop("Provided 'fill_by' is not a column of provided data.frame.")}

  df = df %>%
    arrange(desc(Combined.Score)) %>%
    slice(1:ntop)

  df$Term = substr(df$Term, start = 1, stop = term_char)
  df$Term = factor(df$Term, levels = rev(df$Term))

  plot = ggplot(df)+
    geom_col(aes(x=Combined.Score, y=Term, fill = get(fill_by)))+
    labs(title = title, fill = fill_by)+
    theme_prism()+
    theme(legend.title = element_text(),
          plot.title = element_text())

  print(plot)
}
