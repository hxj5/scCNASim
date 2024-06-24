# From BAM file to a VCF file containing phased SNPs

The `pipeline.py` script is used for generating reference-phased SNPs from
sequence alignments.
Its input is a BAM file, from RNA-seq, DNA-seq, ATAC-seq, either bulk or 
single cells, where for single cells, it supports both droplet-based, e.g.,
10x Genomics, and well-based, e.g., SMART-seq.
Its output is a VCF file containing reference phased SNPs.


## Dependencies

### Softwares

To use the pipeline, please install
[python (tested on python 3.11)](https://www.python.org/) and 
[xcltk >= 0.3.1][xcltk repo], together with a few dependencies listed below.

- [bcftools][bcftools]
- [bgzip or htslib][htslib]
- [cellsnp-lite >= 1.2.0][cellsnp-lite]
- [eagle2][eagle2]

Please install these softwares and add them to the system search path (i.e.,
the system variable `PATH`).

#### Use conda environment (recommended)

To make the installation easier, we highly recommend to set up a conda env and
then install the softwares in the env.

```
# Notes:
# 1. `pandas` supports python 3.9-3.12 while `pysam` supports python 3.6-3.11,
#     therefore, we tested with python 3.11;
# 2. bgzip/htslib will be automatically installed when installing bcftools
#    via conda.
# 3. Eagle2 has to be manually installed since there is not corresponding
#    installation in conda.

conda create -n xcltk python=3.11
conda activate xcltk
conda install -c conda-forge -c bioconda bcftools cellsnp-lite
pip install 'xcltk>=0.3.1'
```

Importantly, [eagle2][eagle2] has to be manually installed since there is not 
corresponding package in conda.
You can download the
[Eagle software](https://alkesgroup.broadinstitute.org/Eagle/#x1-30002)
and then unzip the file, e.g.,

```
# Commands below download and unzip Eagle v2.4.1;
# Download the latest version at https://alkesgroup.broadinstitute.org/Eagle/

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar xzvf Eagle_v2.4.1.tar.gz
```

### Files

In addition, the pipeline relies on a common SNP VCF (for pileup), 
phasing reference panel and genetic map (for phasing with eagle2).
You can download the files using following links.

#### 1000G SNP VCF

```
# hg38
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

# hg19
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz
```

#### 1000G Reference Panel

The pre-compiled files below will be used by Eagle2 as reference panel for
SNP phasing.
Credits to the authors of [Numbat][Numbat].

```
# hg38
wget http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip

# hg19
wget http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip
```

#### Genetic Map

Use commands in above section to download zipped Eagle2 file.
After unzip, the genetic map files are in subfolder `tables`,
e.g., `Eagle_v2.4.1/tables`.


## Quick Start

Note that you may use reference cells, instead of all valid cells 
(which is the default), for genotyping, hopefully to reduce the number of 
false positive calls, such as somatic variants.
To do that, please only specify reference cells in `--barcode` (for 10x data);
or in `--samList` and `--sampleList` (for SMART-seq data).

### 10x scRNA-seq data

An example for **hg38** data.
Type `python pipeline.py --help` for full parameters.

```
# conda activate xcltk

python pipeline.py  \
    --label        {sample name}    \
    --sam          {BAM file}       \
    --barcode      {barcode file}   \
    --snpvcf       {genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz}  \
    --outdir       {output folder}          \
    --gmap         {Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz}  \
    --eagle        {Eagle_v2.4.1/eagle}      \
    --paneldir     {1000G_hg38}              \
    --ncores       10
```


### SMART-seq data

An example for **hg19** data.
Type `python pipeline.py --help` for full parameters.

```
# conda activate xcltk

python pipeline.py  \
    --label        {sample name}           \
    --samList      {BAM list file}         \
    --sampleList   {sample ID list file}   \
    --snpvcf       {genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz}  \
    --outdir       {output folder}          \
    --gmap         {Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz}  \
    --eagle        {Eagle_v2.4.1/eagle}      \
    --paneldir     {1000G_hg19}              \
    --cellTAG      None         \
    --UMItag       None         \
    --ncores       10
```


## Pipeline walkthrough

The phasing pipeline mainly includes two steps, while each step has been 
implemented in a specific function in the `xcltk` repo.

#### 1. Germline SNPs calling

The pipeline firstly calls germline heterozygous SNPs using a pseudo-bulk 
strategy, i.e., by aggregating UMIs/reads of all cells as one bulk sample.
By default, it only keeps SNPs with `--minCOUNT 20` and `--minMAF 0.1`.

This step is implemented in the function `pileup()` in module
`xcltk.baf.genotype`.

#### 2. Reference based SNP phasing

Then the pipeline performs SNP phasing with Eagle2 using 1000G reference panel.

This step is implemented in the function `ref_phasing()` in module
`xcltk.baf.genotype`.



[bcftools]: https://github.com/samtools/bcftools
[cellsnp-lite]: https://github.com/single-cell-genetics/cellsnp-lite
[eagle2]: https://alkesgroup.broadinstitute.org/Eagle/
[htslib]: https://github.com/samtools/htslib
[Numbat]: https://github.com/kharchenkolab/numbat/
[xcltk repo]: https://github.com/hxj5/xcltk