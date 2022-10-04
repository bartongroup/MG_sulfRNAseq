make_lograt_donors <- function(set) {
  dat <- set$dat |> 
    filter(good) |> 
    left_join(set$metadata, by = "sample") |> 
    pivot_wider(id_cols = c(gene_id, replicate), names_from = group, values_from = rlog) |> 
    mutate(
      lograt = LSF - CTRL,
      sample = str_glue("D{replicate}")
    ) |> 
    select(gene_id, sample, lograt)
  meta <- tibble(
    sample = unique(dat$sample)
  )
  
  list(
    metadata = meta,
    dat = dat
  )
}

