# Global constants, design of the experiment

PATH <-  "rna_seq"
SAMPLE_FILE <- "rna_seq/config/samples.txt"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
ENSEMBL_VERSION <- "107"
KEGG_SPECIES <- "hsa"
GO_SPECIES <- "goa_human"
REACTOME_SPECIES <- "Homo sapiens"

EXAMPLE <- "WT_0_1"

LOGFC_LIMIT <- 1
FDR_LIMIT <- 0.01

CONTRAST_SELECTION <-
  c(
    "CTRL-LSF"
  )


SAMPLE_RENAME <- tibble(
  raw_sapmle = c("SRR12926844", "SRR12926848", "SRR12926846", "SRR12926842", "SRR12926843", "SRR12926845", "SRR12926847", "SRR12926849"),
  treatment = c("CTRL", "CTRL", "CTRL", "CTRL", "LSF", "LSF", "LSF", "LSF"),
  replicate = c(2, 4, 3, 1, 1, 2, 3, 4)
)


make_metadata <- function(sren) {
  sren |> 
    unite(sample, c(treatment, replicate), sep = "-", remove = FALSE) |> 
    mutate(group = treatment) |> 
    mutate(
        across(c(sample, replicate, group), as_factor)
    )
}

make_dirs <- function(top_dir) {
  list(
    starmap = file.path(top_dir, "starmap"),
    fscreen = file.path(top_dir, "fscreen"),
    readcount = file.path(top_dir, "readcount"),
    chrcount = file.path(top_dir, "chrcount"),
    bedgraph = file.path(top_dir, "bedgraph"),
    qc = file.path(top_dir, "qc")
  )
}

# Gene names in the W303 genome have a suffix "_W303" or "mRNA". Need to remove
# these before processing.
fix_gene_names <- function(v, suffixes = c("_W303", "_mRNA")) {
  expr <- str_c(suffixes, collapse = "|")
  v |> 
    str_remove(expr)
}
