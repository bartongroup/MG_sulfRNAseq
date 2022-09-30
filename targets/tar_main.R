targets_main <- function() {
  
  # biomart annotations
  get_annotations <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(gns, biomart_fetch_genes(mart)),
    tar_target(bm_gene_lengths, biomart_fetch_gene_lengths(mart, gns$gene_id)),
    tar_target(bm_genes, gns |> left_join(select(bm_gene_lengths, -chr), by = "gene_id")),
    tar_target(terms, get_functional_terms(mart, bm_genes, "human")),
    tar_target(fterms, prepare_functional_terms(terms, star$genes$gene_id))
  )

  # directories and metadata
  setup_experiment <- list(
    tar_target(rnaseq_dirs, make_dirs(PATH)),
    tar_target(metadata, make_metadata(SAMPLE_RENAME))
  )

  # quality control
  qc <- list(
    tar_target(fscreen, parse_fscreens(rnaseq_dirs, metadata, suffix = "_screen.txt")),
    tar_target(qcs, parse_qcs(rnaseq_dirs, metadata, paired = FALSE)),
    tar_target(idxstats, parse_idxstats(rnaseq_dirs, metadata)),
    tar_target(tab_star_log, star$star_log),
    
    tar_target(fig_fscreen, plot_fscreen_map(fscreen)),
    tar_target(fig_fscreen_sample, plot_fscreen_sample(fscreen, EXAMPLE)),
    tar_target(fig_read_qual, plot_qualities(qcs)),
    tar_target(fig_read_qual_clust, plot_cluster_qualities(qcs)),
    # tar_target(fig_chrom_proportion, plot_chrom_proportion(idxstats)),
    
    tar_target(fig_star_log, plot_star_log(star$star_log, metadata)),
    tar_target(fig_star_log_map, plot_star_log_map(star$star_log, metadata)),
    tar_target(fig_map_count, plot_mapped_count(star)),

    tar_target(example_count_file, one_count_file(rnaseq_dirs, 1)),
    tar_target(fig_star_sense, plot_star_sense(example_count_file)),
    
    tar_target(house, read_houskeeping_genes("info"))
  )
    
  # read star counts
  read_data <- list(
    tar_target(star, read_and_process_star(rnaseq_dirs, metadata, bm_genes, min.count = 10))
  )
  
  selections <- list(
    tar_target(gene2name, set_names(star$genes$gene_symbol, star$genes$gene_id))
  )

  # read count properties
  count_properties <- list(
    tar_target(fig_sample_dist, plot_sample_distributions(star, x_breaks = c(0, 1, 2, 3), x_lim = c(-1, 3), ncol = 4, colour_var = "treatment")),
    tar_target(fig_mean_var, plot_mean_var(star)),
    tar_target(fig_distance_mat, plot_distance_matrix(star)),
    tar_target(fig_clustering, plot_clustering(star)),
    tar_target(fig_pca, plot_pca(star, colour_var = "treatment", shape_var = "replicate")),
    tar_target(fig_umap, plot_umap(star, n_neighbours = 5, min_dist = 0.1, colour_var = "treatment", shape_var = "replicate"))
  )
  
  # differential expression
  differential_expression <- list(
    # DE pairwise group
    tar_target(edger, edger_de(star, gns, fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT)),
    tar_target(deseq, deseq2_de(star, gns, fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT)),
    tar_target(de_cmp, edger_deseq_merge(edger, deseq)),
    
    tar_target(de_genes, edger |> filter(sig) |> pull(gene_id) |> unique()),
    tar_target(de_genes_up, edger |> filter(sig & logFC > 0) |> pull(gene_id) |> unique()),

    # upset
    tar_target(upset_list_edger_fc1, de_list(edger, "contrast", "FDR", "logFC")),
    tar_target(upset_list_deseq_fc1, de_list(deseq, "contrast", "FDR", "logFC")),
    tar_target(upset_list_edger_vs_deseq, de_list(de_cmp, "tool", "FDR", "logFC")),
    
    # DE figures
    tar_target(fig_updown, plot_up_down(edger)),
    tar_target(fig_volcano, plot_volcano(edger)),
    tar_target(fig_ma, plot_ma(edger))
  )
  
  set_enrichment <- list(
    tar_target(gse, fgsea_all_terms(edger, fterms)),
    tar_target(gse_all, map_dfr(gse, identity)),
    tar_target(gse_edger, gsea_de(gse, edger))

    # tar_target(fig_fg_example_go_0030476, plot_fgsea_enrichment("GO:0030476", edger_sel |> filter(contrast == "Tfe2_60-WT_60"), fterms$go)),
  )
  
  for_report <- list(
    tar_target(fig_volcano_house, plot_volcano_house(edger, house)),
    
    tar_target(genes_of_interest, bm_genes |> filter(gene_symbol %in% GENES_OF_INTERST) |> pull(gene_id)),
    tar_target(fig_genes, plot_gene_groups(star, genes_of_interest, ncol = 4)),
    
    tar_target(gse_tab, gse_edger |> group_by(term_id, term_name, NES) |> summarise(genes = str_c(gene_symbol, collapse = ", ")) |> arrange(NES) |> ungroup()),
    
    tar_target(terms_cytokine, get_terms_str(gse_all, "[Cc]ytokine")),
    tar_target(terms_interleukin, get_terms_str(gse_all, "[Ii]nterleukin")),
    tar_target(terms_b_cell, get_terms_str(gse_all, "B cell")),
    tar_target(terms_monocyte, get_terms_str(gse_all, "[Mm]onocyte")),
    tar_target(terms_oxidative_stress, get_terms_str(gse_all, "[Oo]xidative stress", fdr_limit = 1)), # none significant
    
    tar_target(terms_up, gse_all |> filter(NES > 0 & padj < 0.05) |> pull(term_id)),
    
    tar_target(png_cytokine, plot_volcano_term(edger, gse_all, terms_cytokine) |> gs("terms_cytokine", 14, 20)),
    tar_target(png_interleukin, plot_volcano_term(edger, gse_all, terms_interleukin) |> gs("terms_interleukin", 14, 64)),
    tar_target(png_b_cell, plot_volcano_term(edger, gse_all, terms_b_cell)  |> gs("terms_b_cell", 14, 8)),
    tar_target(png_monocyte, plot_volcano_term(edger, gse_all, terms_monocyte)  |> gs("terms_monocyte", 14, 4)),
    tar_target(png_oxi, plot_volcano_term(edger, gse_all, terms_oxidative_stress)  |> gs("terms_oxidative_stress", 14, 20)),
    tar_target(png_up, plot_volcano_term(edger, gse_all, terms_up)  |> gs("terms_up", 14, 28))
  )
  
  
  # weighted gene co-expression network analysis
  wgcna <- list(
    tar_target(co_tab, wgcna_prepare(star)),
    tar_target(co_power, wgcna_thresholds(co_tab)),
    tar_target(co_net, wgcna_net(co_tab, power = 12)),
    tar_target(co_enr, wgcna_colour_enrichment(co_net, fterms, gene2name)),
    tar_target(co_edges, wgcna_network(co_tab, co_net, star$genes)),
    tar_target(test_net, plot_network(co_edges, 1))
  )
  
  make_tables <- list(
    tar_target(sav_counts, save_count_data(star, "tab/normalised_counts.tsv")),
    tar_target(sav_de_p, save_de(edger, "tab/de.tsv")),
    tar_target(sav_gse, gse_edger |> mutate(across(where(is.numeric), ~signif(.x, 4))) |> write_tsv("tab/gse.tsv"))
  )
  
  info <- list(
    tar_target(edger_version, packageVersion("edgeR")),
    tar_target(fgsea_version, packageVersion("fgsea"))
  )
  
  c(
    get_annotations,
    setup_experiment,
    selections,
    qc,
    read_data,
    count_properties,
    differential_expression,
    set_enrichment,
    #wgcna,
    make_tables,
    for_report,
    info
  )
}
