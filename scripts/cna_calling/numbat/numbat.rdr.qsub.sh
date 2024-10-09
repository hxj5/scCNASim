#!/bin/bash
#PBS -N numbat_rdr
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o numbat_rdr.out
#PBS -e numbat_rdr.err

source ~/.bashrc
conda activate numbat

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux

repo_dir=~/projects/scCNASim

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
  work_dir=$PBS_O_WORKDIR
fi

out_dir=$work_dir/result

cp  $repo_dir/scripts/cna_calling/numbat/numbat.rdr.R  $work_dir

# settings
sid=GX109
genome=hg38
ncores=10

count_mtx_dir=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/matrix/helen_filtered_matrices

cell_anno_fn=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/anno/GX109-T1c_scRNA_annotation_2column.tsv
ref_cell_type="immune cells"


# preprocess args
prefix=${sid}.numbat_rdr


#Rscript $work_dir/numbat.rdr.R  \
#  <count matrix dir>  \
#  <cell anno file>   \
#  <ref cell type>    \
#  <out dir>         \
#  <out prefix>      \
#  <genome version>   \
#  <ncores>

/usr/bin/time -v Rscript $work_dir/numbat.rdr.R  \
    $count_mtx_dir      \
    $cell_anno_fn    \
    "$ref_cell_type"       \
    $out_dir    \
    $prefix         \
    $genome         \
    $ncores


set +ux
conda deactivate
echo All Done!