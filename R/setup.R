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

CHROMOSOMES <- c(1:22, "X", "Y")


SAMPLE_RENAME <- tibble(
  raw_sample = c("SRR12926844", "SRR12926848", "SRR12926846", "SRR12926842", "SRR12926843", "SRR12926845", "SRR12926847", "SRR12926849"),
  treatment = c("CTRL", "CTRL", "CTRL", "CTRL", "LSF", "LSF", "LSF", "LSF"),
  replicate = c(2, 4, 3, 1, 1, 2, 3, 4)
)


GENES_OF_INTERST <- c("NFE2L2", "IL6", "IL1B", "TNF", "CXCL1", "CCL1", "CCL2", "MMP9", "MMP14", "SOCS3", "CD74", "GPX2", "TXN", "TXNRD1", "SRXN1")

GO_EXAMPLES <- c("GO:0005125", "GO:0006952", "GO:0006954", "GO:0009986", "GO:0030335", "GO:0071356")
KG_EXAMPLES <- c("hsa04668", "hsa04060", "hsa04662")


make_metadata <- function(sren) {
  sren |> 
    unite(sample, c(treatment, replicate), sep = "-", remove = FALSE) |> 
    mutate(group = treatment) |> 
    arrange(treatment, replicate) |> 
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
