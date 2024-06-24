#!/bin/bash
#PBS -N infercnv
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o infercnv.out
#PBS -e infercnv.err

source ~/.bashrc
conda activate XCLBM

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux

repo_dir=~/projects/scCNVSim

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


out_dir=$work_dir/result
if [ ! -e "$out_dir" ]; then
    mkdir -p $out_dir
fi

cp  $repo_dir/scripts/cnv_calling/infercnv/infercnv.R  $work_dir

#Rscript $work_dir/infercnv.R \
#  <sample id>       \
#  <matrix dir>      \
#  <cell anno file>     \
#  <ref cell type>      \
#  <gene anno file>     \
#  <out dir>            \
#  <sequencing platform>: "10x" or "smartseq"   \
#  <number of threads>

/usr/bin/time -v Rscript $work_dir/infercnv.R \
    GX109  \
    /groups/cgsd/xianjie/data/dataset/GX109/scRNA/matrix/helen_filtered_matrices  \
    /groups/cgsd/xianjie/data/dataset/GX109/scRNA/anno/GX109-T1c_scRNA_annotation_2column.tsv  \
    "immune cells"     \
    /groups/cgsd/xianjie/data/refapp/infercnv/hg38_gene_note_noheader_unique.txt  \
    $out_dir        \
    "10x"           \
    10


set +ux
conda deactivate
echo All Done!