
fgsea_run <- function(term_data, res, min.size = 3) {
  res <- res |>
    filter(!is.na(value) & !is.na(feature_id)) |> 
    mutate(first_id = str_remove(feature_id, ";.*$"))
  ranks <-  set_names(res$value, res$first_id)
  fgsea::fgsea(pathways = term_data$term2feature, stats = ranks, nproc = 6, minSize = min.size, eps = 0) |>
    as_tibble() |> 
    rename(term_id = pathway) |> 
    mutate(term_name = term_data$term2name[term_id]) |> 
    arrange(NES) |>
    select(term_id, term_name, pval, padj, NES, size, leading_edge = leadingEdge)
}

fgsea_groups <- function(d, term_data, feature_var, value_var, group_var) {
  d |>
    mutate(value = get(value_var), feature_id = get(feature_var)) |>
    group_split(!!sym(group_var)) |>
    map_dfr(function(w) {
      fgsea_run(term_data, w) |>
        mutate(!!group_var := dplyr::first(w[[group_var]]))
    })
}

fgsea_all_terms <- function(d, terms, feature_var = "gene_id", value_var = "logFC", group_var = "contrast") {
  ontologies <- names(terms)
  map(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    fgsea_groups(d, terms[[ont]], feature_var, value_var, group_var)
  }) |>
    set_names(ontologies)
}





plot_fgsea_enrichment <- function(term_id, res, terms, valvar = "logFC") {
  lst <- terms$term2feature[[term_id]]
  rnks <- set_names(res[[valvar]], res$gene_id)
  fgsea::plotEnrichment(lst, rnks)
}

split_genes_fgsea <- function(se, fg, groupvar = "contrast") {
  fg |> 
    filter(padj < 0.05) |> 
    group_split(term, !!sym(groupvar)) |> 
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se |> 
        filter(gene_id %in% genes & !!sym(groupvar) == gr) |> 
        mutate(term_id = term, .before = "gene_id")
    })
}


gsea_de <- function(gse, res, fdr_limit = 0.01) {
  sig_genes <- res |> 
    filter(sig) |> 
    pull(gene_id)
  ontologies <- names(gse)
  map_dfr(ontologies, function(ont) {
    gse[[ont]] |> 
      filter(padj < fdr_limit) |> 
      unnest(leading_edge) |> 
      filter(leading_edge %in% sig_genes) |> 
      rename(gene_id = leading_edge) |> 
      left_join(res, by = "gene_id") |> 
      select(term_id, term_name, NES, padj, gene_id, gene_symbol, logFC, logCPM, PValue, FDR) |> 
      add_column(ontology = ont, .before = 1)
  })
}


get_terms_str <- function(gso, query, fdr_limit = 0.05) {
  gso |> 
    filter(str_detect(term_name, query) & padj < fdr_limit) |>
    pull(term_id)
}