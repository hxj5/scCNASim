#!/bin/bash
#PBS -N numbat
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o numbat.out
#PBS -e numbat.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux

repo_dir=~/projects/scCNASim

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
  work_dir=$PBS_O_WORKDIR
fi

out_dir=$work_dir/result

cp  $repo_dir/scripts/cna_calling/numbat/pileup_and_phase.R  $work_dir
cp  $repo_dir/scripts/cna_calling/numbat/numbat.R  $work_dir

# settings
sid=GX109
platform=10x              # 10x, smartseq, or bulk
genome=hg38
ncores=10

bams=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/bam/raw_GX109-T1c/possorted_genome_bam.bam
barcodes=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/matrix/helen_filtered_matrices/barcodes.tsv
count_mtx_dir=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/matrix/helen_filtered_matrices

cell_anno_fn=/groups/cgsd/xianjie/data/dataset/GX109/scRNA/anno/GX109-T1c_scRNA_annotation_2column.tsv
ref_cell_type="immune cells"

gmap=/groups/cgsd/xianjie/data/refapp/eagle/Eagle_v2.4.1/tables/genetic_map_${genome}_withX.txt.gz
eagle=~/.anaconda3/envs/XCLBM/bin/eagle
snpvcf=/groups/cgsd/xianjie/data/refapp/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.${genome}.vcf.gz
paneldir=/groups/cgsd/xianjie/data/refapp/numbat/1000G_${genome}


# preprocess args
prefix=${sid}.numbat

platform_option=
if [ "$platform" = "10x" ]; then
    platform_option=""
elif [ "$platform" = "smartseq" ]; then
    platform_option="--smartseq"
elif [ "$platform" = "bulk" ]; then
    platform_option = "--bulk"
else
    echo "invalid platform '$platform'"
    exit 1
fi


#Run SNP pileup and phasing with 1000G
if [ ! -e "$out_dir/allele" ]; then
    mkdir -p $out_dir/allele
fi

#usage: pileup_and_phase.R [-h] --label LABEL --samples SAMPLES --bams BAMS
#                          [--barcodes BARCODES] --gmap GMAP [--eagle EAGLE]
#                          --snpvcf SNPVCF --paneldir PANELDIR --outdir OUTDIR
#                          --ncores NCORES [--UMItag UMITAG]
#                          [--cellTAG CELLTAG] [--smartseq] [--bulk]
#Arguments:
#  -h, --help           show this help message and exit
#  --label LABEL        Individual label
#  --samples SAMPLES    Sample names, comma delimited
#  --bams BAMS          BAM files, one per sample, comma delimited
#  --barcodes BARCODES  Cell barcode files, one per sample, comma delimited
#  --gmap GMAP          Path to genetic map provided by Eagle2 (e.g.
#                       Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)
#  --eagle EAGLE        Path to Eagle2 binary file
#  --snpvcf SNPVCF      SNP VCF for pileup
#  --paneldir PANELDIR  Directory to phasing reference panel (BCF files)
#  --outdir OUTDIR      Output directory
#  --ncores NCORES      Number of cores

/usr/bin/time -v Rscript $work_dir/pileup_and_phase.R   \
    --label  $sid    \
    --samples  $sid  \
    --bams  $bams   \
    --barcodes  $barcodes   \
    --gmap  $gmap   \
    --eagle  $eagle  \
    --snpvcf  $snpvcf  \
    --paneldir  $paneldir  \
    --outdir  $out_dir/allele  \
    --ncores  $ncores          \
    $platform_option


#Rscript $work_dir/numbat.R  \
#  <allele file>       \
#  <count matrix dir>  \
#  <cell anno file>   \
#  <ref cell type>    \
#  <out dir>         \
#  <out prefix>      \
#  <genome version>   \
#  <sequencing platform>  \
#  <ncores>

/usr/bin/time -v Rscript $work_dir/numbat.R  \
    $out_dir/allele/${sid}_allele_counts.tsv.gz    \
    $count_mtx_dir      \
    $cell_anno_fn    \
    "$ref_cell_type"       \
    $out_dir/cna    \
    $prefix         \
    $genome         \
    $platform       \
    $ncores


set +ux
conda deactivate
echo All Done!