# numbat.rdr.R

app <- "numbat.rdr.R"

args <- commandArgs(T)
if (length(args) < 7) {
  msg <- paste0("Usage: ", app, " <count matrix dir>  \\
                                  <cell anno file> <ref cell type>  \\
                                  <output dir> <out prefix>   \\
                                  <genome version> <ncores>")
  write(msg, file = stderr())
  quit("no", 1)
}

count_mtx_dir <- args[1]
cell_anno_fn <- args[2]
ref_cell_type <- args[3]
out_dir <- args[4]
out_prefix <- args[5]
genome <- args[6]
ncores <- as.numeric(args[7])


library(Seurat)
library(numbat)

# check args
if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = T)
}
setwd(out_dir)

if (! genome %in% c("hg19", "hg38"))
  stop(sprintf("invalid genome version '%s'.", genome))

# filter ref cells from count matrix
cell_anno <- read.delim(cell_anno_fn, header = F, stringsAsFactors = F)
colnames(cell_anno) <- c("cell", "group")
ref_cells <- cell_anno[cell_anno$group %in% ref_cell_type, ]
str(ref_cells)

count_mtx <- Seurat::Read10X(data.dir = count_mtx_dir)   # gene x cell
str(count_mtx)

cnt_mtx1 <- count_mtx[, ! (colnames(count_mtx) %in% ref_cells$cell)]
str(cnt_mtx1)
saveRDS(cnt_mtx1, sprintf("%s/%s.ref_filtered.count.mtx.rds", out_dir,
                          out_prefix))

# calculate average expression of ref cells
ref_mean_expr <- numbat::aggregate_counts(count_mtx, ref_cells)
str(ref_mean_expr)
saveRDS(ref_mean_expr, sprintf("%s/%s.ref.gene_by_celltype.mtx.rds",
                               out_dir, out_prefix))

obj <- run_numbat_rdr(
  cnt_mtx1,          # gene x cell integer UMI count matrix 
  ref_mean_expr,     # reference expression profile, a gene x cell type normalized expression level matrix
  genome = genome,
  ncores = ncores,
  plot = TRUE,
  out_dir = out_dir
)

saveRDS(obj, sprintf("%s/%s.out.object.rds", out_dir, out_prefix))

print(paste0("[", app, "] All Done!"))
