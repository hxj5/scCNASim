# pipeline.py - preprocess the input BAM file to generate phased SNPs in VCF format.

# Note:
# 1. - this file was copied and modifed from `xcltk/baf/pipeline.py` in `xcltk` repo.


import getopt
import os
import sys

from logging import error, info

# requires xcltk>=0.3.0
from xcltk.baf.genotype import pileup, ref_phasing, vcf_add_genotype
from xcltk.utils.base import assert_e, assert_n
from xcltk.utils.vcf import vcf_index, vcf_merge, vcf_split_chrom
from xcltk.utils.xlog import init_logging


APP = "python"
VERSION = "0.0.1"

COMMAND = "pipeline.py"

CELL_TAG = "CB"
N_CORES = 1
UMI_TAG = "UB"
MIN_COUNT = 20
MIN_MAF = 0.1


def usage(fp = sys.stdout):
    s =  "\n" 
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s [options]\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  --label STR        Task label.\n"
    s += "  --sam FILE         Comma separated indexed BAM/CRAM file(s).\n"
    s += "  --samList FILE     A file listing BAM/CRAM files, each per line.\n"
    s += "  --barcode FILE     A plain file listing all effective cell barcodes, for\n"
    s += "                     droplet-based data, e.g., 10x Genomics.\n"
    s += "  --sampleList FILE  A plain file listing sample IDs, one ID per BAM, for\n"
    s += "                     well-based data, e.g., SMART-seq.\n"
    s += "  --snpvcf FILE      A VCF file listing all candidate SNPs.\n"
    s += "  --outdir DIR       Output dir.\n"
    s += "  --gmap FILE        Path to genetic map provided by Eagle2\n"
    s += "                     (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz).\n"
    s += "  --eagle FILE       Path to Eagle2 binary file.\n"
    s += "  --paneldir DIR     Directory to phasing reference panel (BCF files).\n"
    s += "  --version          Print version and exit.\n"
    s += "  --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  --cellTAG STR      Cell barcode tag; Set to None if not available [%s]\n" % CELL_TAG
    s += "  --UMItag STR       UMI tag; Set to None if not available [%s]\n" % UMI_TAG
    s += "  --minCOUNT INT     Minimum aggregated UMI or read count [%d]\n" % MIN_COUNT
    s += "  --minMAF INT       Minimum minor allele frequency [%f]\n" % MIN_MAF
    s += "  --ncores INT       Number of threads [%d]\n" % N_CORES
    s += "\n"
    s += "Notes:\n"
    s += "1. One and only one of `--sam` and `--samlist` should be specified.\n"
    s += "2. For well-based data, the order of the BAM files (in `--sam` or `--samlist`)\n"
    s += "   and the sample IDs (in `--sampleList`) should match each other.\n"
    s += "3. For bulk data, the label (`--label`) will be used as the sample ID.\n"
    s += "\n"

    fp.write(s)


def pipeline_main(argv):
    if len(argv) <= 1:
        usage()
        sys.exit(0)

    init_logging(stream = sys.stdout)

    label = None
    sam_fn = sam_list_fn = None
    barcode_fn = sample_id_fn = None
    snp_vcf_fn = None
    out_dir = None
    gmap_fn = eagle_fn = panel_dir = None
    cell_tag, umi_tag = CELL_TAG, UMI_TAG
    min_count, min_maf = MIN_COUNT, MIN_MAF
    ncores = N_CORES

    opts, args = getopt.getopt(
        args = argv[1:],
        shortopts = "", 
        longopts = [
            "label=",
            "sam=", "samList=", 
            "barcode=", "sampleList=",
            "snpvcf=",
            "outdir=",
            "gmap=", "eagle=", "paneldir=",
            "version", "help",
            
            "cellTAG=", "UMItag=",
            "minCOUNT=", "minMAF=",
            "ncores="
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in ("--label"): label = val
        elif op in ("--sam"): sam_fn = val
        elif op in ("--samlist"): sam_list_fn = val
        elif op in ("--barcode"): barcode_fn = val
        elif op in ("--samplelist"): sample_id_fn = val
        elif op in ("--snpvcf"): snp_vcf_fn = val
        elif op in ("--outdir"): out_dir = val
        elif op in ("--gmap"): gmap_fn = val
        elif op in ("--eagle"): eagle_fn = val
        elif op in ("--paneldir"): panel_dir = val
        elif op in ("--version"): sys.stdout.write(VERSION + "\n"); sys.exit(0)
        elif op in ("--help"): usage(); sys.exit(0)

        elif op in ("--celltag"): cell_tag = val
        elif op in ("--umitag"): umi_tag = val
        elif op in ("--mincount"): min_count = int(val)
        elif op in ("--minmaf"): min_maf = float(val)
        elif op in ("--ncores"): ncores = int(val)     # keep it in `str` format.
        else:
            error("invalid option: '%s'." % op)
            return(-1)
        
    ret = pipeline_wrapper(
        label = label,
        sam_fn = sam_fn, sam_list_fn = sam_list_fn, 
        barcode_fn = barcode_fn, sample_id_fn = sample_id_fn,
        snp_vcf_fn = snp_vcf_fn,
        out_dir = out_dir,
        gmap_fn = gmap_fn, eagle_fn = eagle_fn, panel_dir = panel_dir,
        cell_tag = cell_tag, umi_tag = umi_tag,
        min_count = min_count, min_maf = min_maf,
        ncores = ncores
    )
    
    info("All Done!")

    return(ret)
        

def pipeline_wrapper(
    label,
    sam_fn = None, sam_list_fn = None, 
    barcode_fn = None, sample_id_fn = None,
    snp_vcf_fn = None,
    out_dir = None,
    gmap_fn = None, eagle_fn = None, panel_dir = None,
    cell_tag = "CB", umi_tag = "UB",
    min_count = 20, min_maf = 0.1,
    ncores = 1
):
    info("phasing pipeline starts ...")

    # check args
    info("check args ...")

    mode = None
    sample_id = None

    assert_n(label)
    assert (not sam_fn) ^ (not sam_list_fn)
    assert not (barcode_fn and sample_id_fn)
    if barcode_fn:
        assert_e(barcode_fn)
        mode = "droplet"
    elif sample_id_fn:
        assert_e(sample_id_fn)
        mode = "well"
    else:
        mode = "bulk"
        sample_id = label

    assert_e(snp_vcf_fn)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)
    script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)

    assert_e(gmap_fn)
    assert_e(eagle_fn)
    assert_e(panel_dir)
    for chrom in range(1, 23):
        assert_e(os.path.join(panel_dir, "chr%d.genotypes.bcf" % chrom))
        assert_e(os.path.join(panel_dir, "chr%d.genotypes.bcf.csi" % chrom))

    genome = "hg19" if "hg19" in gmap_fn else "hg38"

    # other args will be checked in pileup() and ref_phasing().
    info("run in '%s' mode (genome version '%s')" % (mode, genome))

    # pileup
    info("start pileup ...")

    pileup_dir = os.path.join(out_dir, "pileup")
    if not os.path.exists(pileup_dir):
        os.makedirs(pileup_dir, exist_ok = True)
    pileup_script_dir = os.path.join(script_dir, "pileup")
    os.makedirs(pileup_script_dir, exist_ok = True)
    pileup_script = os.path.join(pileup_script_dir, "run_pileup.sh")
    pileup_log_fn = os.path.join(pileup_script_dir, "pileup.log")

    pileup(
        sam_fn = sam_fn, sam_list_fn = sam_list_fn,
        barcode_fn = barcode_fn, sample_id_fn = sample_id_fn, 
        sample_id = sample_id,
        snp_vcf_fn = snp_vcf_fn,
        out_dir = pileup_dir,
        mode = mode,
        cell_tag = cell_tag, umi_tag = umi_tag,
        ncores = ncores,
        min_count = min_count, min_maf = min_maf,
        script_fn = pileup_script,
        log_fn = pileup_log_fn
    )

    pileup_vcf_fn = os.path.join(pileup_dir, "cellSNP.base.vcf.gz")
    assert_e(pileup_vcf_fn)

    info("pileup VCF is '%s'." % pileup_vcf_fn)

    # TODO: further filter SNPs given `minMAF` and `minCOUNT`, considering
    # only REF and ALT (AD & DP) alleles, but not OTH alleles.
    # cellsnp-lite (at least v1.2.3 and before) may keep some SNPs unexpectly,
    # e.g., SNP with `AD=9;DP=9;OTH=1` when `minMAF=0.1; minCOUNT=10`.

    # prepare VCF files for phasing
    info("prepare VCF files for phasing ...")

    phasing_dir = os.path.join(out_dir, "phasing")
    if not os.path.exists(phasing_dir):
        os.makedirs(phasing_dir, exist_ok = True)
    phasing_script_dir = os.path.join(script_dir, "phasing")
    os.makedirs(phasing_script_dir, exist_ok = True)
    phasing_script_prefix = os.path.join(phasing_script_dir, "run_phasing")
    phasing_log_prefix = os.path.join(phasing_script_dir, "phasing")

    # add genotypes
    info("add genotypes ...")

    genotype_vcf_fn = os.path.join(out_dir, "%s.genotype.vcf.gz" % label)
    vcf_add_genotype(
        in_fn = pileup_vcf_fn,
        out_fn = genotype_vcf_fn,
        sample = label,
        chr_prefix = True,     # add "chr" prefix
        sort = True,
        unique = True
    )

    # split VCF by chromosomes.
    info("split VCF by chromosomes ...")

    valid_chroms = []       # has "chr" prefix
    target_vcf_list = []
    res = vcf_split_chrom(
        fn = genotype_vcf_fn,
        out_dir = phasing_dir,
        label = label,
        chrom_list = ["chr" + str(i) for i in range(1, 23)],
        out_prefix_list = None,
        verbose = True
    )
    for chrom, n_variants, vcf_fn in res:
        if n_variants > 0:
            valid_chroms.append(chrom)
            target_vcf_list.append(vcf_fn)
            vcf_index(vcf_fn)

    info("%d chromosome VCFs are outputted with variants." % len(valid_chroms))
    #os.remove(genotype_vcf_fn)

    # reference phasing
    info("reference phasing ...")

    ref_vcf_list = [os.path.join(panel_dir, "%s.genotypes.bcf" % chrom) \
                    for chrom in valid_chroms]
    out_prefix_list = [os.path.join(phasing_dir, "%s_%s.phased" % \
                    (label, chrom)) for chrom in valid_chroms]
    ref_phasing(
        target_vcf_list = target_vcf_list,
        ref_vcf_list = ref_vcf_list,
        out_prefix_list = out_prefix_list,
        gmap_fn = gmap_fn,
        eagle_fn = eagle_fn,
        out_dir = phasing_dir,
        ncores = ncores,
        script_fn_prefix = phasing_script_prefix,
        log_fn_prefix = phasing_log_prefix,
        verbose = True
    )

    info("phased VCFs are in dir '%s'." % phasing_dir)

    # merge phased VCFs
    info("merge phased VCFs ...")

    phased_vcf_list = ["%s.vcf.gz" % prefix for prefix in out_prefix_list]
    for fn in phased_vcf_list:
        assert_e(fn)

    phased_vcf_fn = os.path.join(out_dir, "%s.phased.vcf.gz" % label)
    vcf_merge(phased_vcf_list, phased_vcf_fn, sort = True)

    info("merged VCF is '%s'." % phased_vcf_fn)


if __name__ == "__main__":
    pipeline_main(sys.argv)