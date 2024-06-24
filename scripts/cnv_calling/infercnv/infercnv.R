# infercnv.R

app <- "infercnv.R"

args <- commandArgs(T)
if (length(args) < 8) {
  msg <- paste0("Usage: ", app, " <sample id> <matrix dir> \\
                                  <cell anno file> <ref cell type> \\
                                  <gene anno file> <out dir>   \\
                                  <sequencing platform> <number of threads>")
  write(msg, file = stderr())
  quit("no", 1)
}

sid <- args[1]
matrix_dir <- args[2]
cell_anno_file <- args[3]
ref_cell_type <- args[4]
gene_anno_file <- args[5]
out_dir <- args[6]
platform <- args[7]
n_threads <- as.numeric(args[8])

library(Seurat)
library(infercnv)


# check args
if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = T)
}
setwd(out_dir)

cutoff <- NULL
if (platform == "10x") {
  cutoff <- 0.1
} else if (platform == "smartseq") {
  cutoff <- 1
} else {
  stop(sprintf("invalid platform '%s'.", platform))
}


# run infercnv
gex_mtx <- Seurat::Read10X(data.dir = matrix_dir)    # gene x cell

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = gex_mtx,
  annotations_file = cell_anno_file,
  delim = '\t',
  gene_order_file = gene_anno_file,
  ref_group_names = c(ref_cell_type)
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff,  
  out_dir = out_dir, 
  cluster_by_groups = T,   
  denoise = T,
  HMM = T,
  num_threads = n_threads
)

print(paste0("[", app, "] All Done!"))
