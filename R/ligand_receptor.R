# Ligand-Receptor analysis functions

#'
#'
get_all_degs = function(obj, annotation, condition, reference, experimental, padj.cut=0.05, log2fc.cut=0.1){
  cols_legit = annotation %in% colnames(obj@meta.data) & condition %in% colnames(obj@meta.data)
  if(!cols_legit){stop("Annotation & Condition need to be columns of the object meta.data!")}

  obj = SetIdent(obj, value = annotation)

  idents = unique(obj@meta.data[[annotation]])
  deglist = map(idents,
                ~ FindMarkers(object = obj, subset.ident = .x,
                              group.by = condition,
                              ident.1 = experimental, ident.2 = reference) %>%
                  filter(p_val_adj < padj.cut) %>%
                  filter(abs(avg_log2FC) > log2fc.cut) %>%
                  mutate(cluster = .x) %>%
                  rownames_to_column(var = "gene")
  )
  merged = reduce(deglist, rbind)
  merged = merged %>%
    mutate(id = paste0(gene, "_", cluster))
  return(merged)
}


#'
#'
analyze_LR = function(signatures, db_path){
  db = read_tsv(db_path)

  # Subset signatures to Ls and Rs
  signatures = signatures %>%
    mutate(is_ligand = gene %in% db[["ligand_gene_symbol"]]) %>%
    mutate(is_receptor = gene %in% db[["receptor_gene_symbol"]]) %>%
    filter(is_ligand | is_receptor)

  # Add unique ID to each signature
  signatures$lr = "R"
  signatures$lr[signatures$is_ligand] = "L"
  signatures$id = paste0(signatures$gene, "_", signatures$cluster)

  # Find the possible pairmates of a given gene from the LR database
  find_pairmates = function(gene, db){
    pairmates = db %>%
      filter(ligand_gene_symbol == gene | receptor_gene_symbol == gene) %>%
      select(ligand_gene_symbol, receptor_gene_symbol)
    pairmates = unique(unlist(pairmates))
    pairmates = pairmates[pairmates != gene]
    return(pairmates)
  }

  # Logic so that ligand always comes first in the pair-name
  organize_pairs = function(sig, pm, is.lig){
    if(is.lig){paste0(sig, "-", pm)}
    else{paste0(pm, "-", sig)}
  }

  # Find all pairmates
  potential_pairmates = map(signatures[["gene"]], find_pairmates, db)
  real_pairmates = map(potential_pairmates, ~ .x[.x %in% signatures[["gene"]]])
  names(real_pairmates) = signatures$gene
  n_real_pairmates = map(real_pairmates, length)

  # Covert to IDs and find all pairs
  real_pairmates = map(real_pairmates, ~ signatures$id[signatures$gene %in% .x])
  real_pairs = pmap(list(sig = signatures$id, pm = real_pairmates, is.lig = signatures$is_ligand), organize_pairs)
  signatures$n_pairmates = n_real_pairmates
  signatures$pairs = map(real_pairs, ~ paste0(.x, collapse = ", "))

  # Generate LR table from all interactions, do some parsing
  lr = data.frame(interaction = unique(unlist(real_pairs[n_real_pairmates > 0])))
  lr = lr %>%
    mutate(ligand_id = str_split(interaction, pattern = "-", simplify = T)[,1]) %>%
    mutate(receptor_id = str_split(interaction, pattern = "-", simplify = T)[,2]) %>%
    mutate(ligand = str_split(ligand_id, pattern = "_", simplify = T)[,1]) %>%
    mutate(receptor = str_split(receptor_id, pattern = "_", simplify = T)[,1]) %>%
    mutate(ligand_origin = str_split(ligand_id, pattern = "_", simplify = T)[,2]) %>%
    mutate(receptor_origin = str_split(receptor_id, pattern = "_", simplify = T)[,2])
  return(lr)
}


#'
#'
crossreference_degs = function(lr_table, deg_table, db_path){
  db = read_tsv(db_path)
  lr_table = lr_table %>%
    mutate(from_signature = TRUE)

  # Look for LR pairs regardless of signature status
  deg_lr = analyze_LR(deg_table, db_path = db_path)
  deg_lr = deg_lr %>%
    mutate(both_de = TRUE)

  # Combine
  merged = merge(lr_table, deg_lr, all=TRUE) %>%
    mutate(reflexive = ligand_origin == receptor_origin) %>%
    replace_na(list(from_signature = FALSE, both_de = FALSE))

  # Get FC and p value for all DE'ed ids
  crossreference = function(id, column){
    if (id %in% deg_table$id){
      return(deg_table[[column]][deg_table[["id"]] == id])
    }
    else
      return(NA)
  }

  merged = merged %>%
    mutate(ligand_log2FC = map_dbl(merged$ligand_id, ~ crossreference(.x, "avg_log2FC"))) %>%
    mutate(ligand_padj = map_dbl(merged$ligand_id, ~ crossreference(.x, "p_val_adj"))) %>%
    mutate(receptor_log2FC = map_dbl(merged$receptor_id, ~ crossreference(.x, "avg_log2FC"))) %>%
    mutate(receptor_padj = map_dbl(merged$receptor_id, ~ crossreference(.x, "p_val_adj"))) %>%
    mutate(any_de = !(is.na(ligand_log2FC) & is.na(receptor_log2FC))) %>%
    relocate(any_de, .before = both_de)

  return(merged)
}

