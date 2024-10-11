#!/bin/bash
#PBS -N infercnv
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=200g,walltime=100:00:00
#PBS -o infercnv.out
#PBS -e infercnv.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux

repo_dir=~/projects/scCNASim

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi
work_dir=`cd $work_dir && cd .. && pwd`
script_dir=$work_dir/scripts

out_dir=$work_dir/infercnv
if [ ! -e "$out_dir" ]; then
    mkdir -p $out_dir
fi

cp  $repo_dir/scripts/cna_calling/infercnv/infercnv.R  $script_dir

#Rscript $script_dir/infercnv.R \
#  <sample id>       \
#  <matrix dir>      \
#  <cell anno file>     \
#  <ref cell type>      \
#  <gene anno file>     \
#  <out dir>            \
#  <sequencing platform>: "10x" or "smartseq"   \
#  <number of threads>

/usr/bin/time -v Rscript $script_dir/infercnv.R \
    st_liver  \
    $work_dir/rdr  \
    $work_dir/simu/4_rs/rs.cell_anno.tsv  \
    "normal"     \
    /groups/cgsd/xianjie/data/refapp/infercnv/hg38_gene_note_noheader_unique.txt  \
    $out_dir        \
    "10x"           \
    10


set +ux
conda deactivate
echo All Done!
