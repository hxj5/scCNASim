#!/bin/bash
#PBS -N rdr
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o rdr.out
#PBS -e rdr.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
  work_dir=$PBS_O_WORKDIR
fi
work_dir=`cd $work_dir && cd .. && pwd`


# PUT YOUR CODE HERE
xcltk basefc     \
    --sam          $work_dir/simu/4_rs/bam/rs.possorted.bam  \
    --barcode      $work_dir/simu/4_rs/rs.samples.tsv   \
    --region       /groups/cgsd/xianjie/data/refapp/xclone/xcltk_data/annotate_genes_hg38_update_20230126.txt   \
    --outdir       $work_dir/rdr   \
    --ncores       10


# the "Seurat::Read10X()" expects "genes.tsv" or "barcodes.tsv.gz".
# "Seurat::Read10X()" is used by inferCNV and Numbat for loading matrix data.
cat $work_dir/rdr/features.tsv | \
  awk '{printf("%s\t%s\n", $4, $4)}' > $work_dir/rdr/genes.tsv
mv $work_dir/rdr/features.tsv $work_dir/rdr/xcltk.features.tsv

set +ux
conda deactivate
echo All Done!

