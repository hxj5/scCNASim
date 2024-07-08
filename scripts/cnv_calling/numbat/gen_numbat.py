# gen_numbat.py - generate numbat scripts


import os
import stat


def assert_e(path):
    assert os.path.exists(path)


def gen_qsub_batch(
    sid_lst,
    mtx_dir_lst,
    cell_anno_fn_lst,
    ref_cell_type_lst,
    work_dir_lst,
    run_script,
    repo_dir,
    out_dir_lst = None,        # trick: ["$work_dir/numbat"] * n
    seq_platform = "10x",      # "10x" or "smartseq"
    genome = "hg38",           # "hg38" or "hg19"
    n_threads = 10,
    module = "combined",       # "combined" or "RDR"
    bam_fn_lst = None, barcode_fn_lst = None,
    gmap_fn = None, eagle_fn = None, snp_vcf = None, panel_dir = None,
    conda_env = "SCSC",
    qsub_q = "cgsd", 
    qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
    verbose = True
):
    n = len(sid_lst)
    assert len(mtx_dir_lst) == n
    assert len(cell_anno_fn_lst) == n
    assert len(ref_cell_type_lst) == n
    assert len(work_dir_lst) == n
    if out_dir_lst is None:
        out_dir_lst = [None] * n
    else:
        assert len(out_dir_lst) == n

    assert module in ("combined", "RDR")
    if module == "combined":
        assert len(bam_fn_lst) == n
        assert len(barcode_fn_lst) == n
    else:
        bam_fn_lst = [None] * n
        barcode_fn_lst = [None] * n

    if verbose:
        print("write to qsub scripts ...")

    qsub_script_lst = []
    for sid, mtx_dir, cell_anno_fn, ref_cell_type, work_dir, out_dir, \
        bam_fn, barcode_fn in \
        zip(sid_lst, mtx_dir_lst, cell_anno_fn_lst, ref_cell_type_lst, \
            work_dir_lst, out_dir_lst, bam_fn_lst, barcode_fn_lst
        ):
        qsub_script_fname = "%s.numbat.qsub.sh" % sid
        qsub_script_fpath = os.path.join(work_dir, qsub_script_fname)
        s = gen_qsub(
            sid = sid,
            mtx_dir = mtx_dir,
            cell_anno_fn = cell_anno_fn,
            ref_cell_type = ref_cell_type,
            repo_dir = repo_dir,
            out_dir = out_dir,
            seq_platform = seq_platform,
            genome = genome,
            n_threads = n_threads,
            module = module,
            bam_fn = bam_fn, barcode_fn = barcode_fn,
            gmap_fn = gmap_fn, eagle_fn = eagle_fn, 
            snp_vcf = snp_vcf, panel_dir = panel_dir,
            conda_env = conda_env,
            qsub_N = sid + "_numbat", 
            qsub_q = qsub_q, 
            qsub_l_nodes = qsub_l_nodes, qsub_l_ppn = qsub_l_ppn, 
            qsub_l_mem = qsub_l_mem, qsub_l_walltime = qsub_l_walltime,
            qsub_o = "%s_numbat.out" % sid, 
            qsub_e = "%s_numbat.err" % sid
        )
        with open(qsub_script_fpath, "w") as fp:
            fp.write(s)
        qsub_script_lst.append(qsub_script_fname)
        if verbose:
            print("    %s: %s" % (sid, qsub_script_fpath))
    if verbose:
        print("")

    # write to running script.
    if verbose:
        print("write to running script %s ..." % run_script)
    s  = "#!/bin/bash\n"
    s += "\n"
    for work_dir, qsub_script in zip(work_dir_lst, qsub_script_lst):
        s += "cd %s\n" % work_dir
        s += "qsub %s\n" % qsub_script
        s += "\n"
    s += "echo All Done!\n"
    with open(run_script, "w") as fp:
        fp.write(s)
    st = os.stat(run_script)
    os.chmod(run_script, st.st_mode  | stat.S_IXUSR)

    if verbose:
        print("All Done!")


def gen_qsub(
    sid,
    mtx_dir,
    cell_anno_fn,
    ref_cell_type,
    repo_dir,
    out_dir = None,
    seq_platform = "10x",      # "10x" or "smartseq"
    genome = "hg38",           # "hg38" or "hg19"
    n_threads = 10,
    module = "combined",       # "combined" or "RDR"
    bam_fn = None, barcode_fn = None,
    gmap_fn = None, eagle_fn = None, snp_vcf = None, panel_dir = None,
    conda_env = "SCSC",
    qsub_N = "numbat", qsub_q = "cgsd", 
    qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
    qsub_o = "numbat.out", qsub_e = "numbat.err"
):
    # check args
    assert_e(mtx_dir)
    assert_e(cell_anno_fn)
    assert_e(repo_dir)
    assert seq_platform in ("10x", "smartseq")
    assert genome in ("hg38", "hg19")
    assert module in ("combined", "RDR")
    if module == "combined":
        assert_e(bam_fn)
        assert_e(barcode_fn)
        assert_e(gmap_fn)
        assert_e(eagle_fn)
        assert_e(snp_vcf)
        assert_e(panel_dir)

    assert_e("/bin/bash")
    #assert_e("/home/.bashrc")
    if module == "combined":
        assert_e(os.path.join(repo_dir, "scripts/cnv_calling/numbat/pileup_and_phase.R"))
        assert_e(os.path.join(repo_dir, "scripts/cnv_calling/numbat/numbat.R"))
    else:
        assert_e(os.path.join(repo_dir, "scripts/cnv_calling/numbat/numbat.rdr.R"))
    assert_e("/usr/bin/time")


    # generate script
    s = ""

    s += '''#!/bin/bash\n'''
    s += '''#PBS -N %s\n''' % qsub_N
    s += '''#PBS -q %s\n''' % qsub_q
    s += '''#PBS -l nodes=%d:ppn=%d,mem=%dg,walltime=%d:00:00\n''' % (\
        qsub_l_nodes, qsub_l_ppn, qsub_l_mem, qsub_l_walltime)
    s += '''#PBS -o %s\n''' % qsub_o
    s += '''#PBS -e %s\n''' % qsub_e
    s += '''\n'''

    s += '''source ~/.bashrc\n'''
    s += '''conda activate %s\n''' % conda_env
    s += '''\n'''

    s += '''# run `set` after `source` & `conda activate` as the source file has an unbound variable\n'''
    s += '''set -eux\n'''
    s += '''\n'''

    s += '''repo_dir=%s\n''' % repo_dir
    s += '''\n'''

    s += '''work_dir=`cd $(dirname $0) && pwd`\n'''
    s += '''if [ -n "$PBS_O_WORKDIR" ]; then\n'''
    s += '''    work_dir=$PBS_O_WORKDIR\n'''
    s += '''fi\n'''
    s += '''\n'''

    if out_dir is None:
        s += '''out_dir=$work_dir/result\n'''
    else:
        s += '''out_dir=%s\n''' % out_dir      
    s += '''if [ ! -e "$out_dir" ]; then\n'''
    s += '''    mkdir -p $out_dir\n'''
    s += '''fi\n'''
    s += '''\n'''

    if module == "combined":
        s += '''cp  $repo_dir/scripts/cnv_calling/numbat/pileup_and_phase.R  $work_dir\n'''
        s += '''cp  $repo_dir/scripts/cnv_calling/numbat/numbat.R  $work_dir\n'''
    else:
        s += '''cp  $repo_dir/scripts/cnv_calling/numbat/numbat.rdr.R  $work_dir\n'''
    s += '''\n'''

    s += '''sid=%s\n''' % sid
    if module == "combined":
        s += '''platform=%s              # 10x, smartseq, or bulk\n''' % seq_platform
    s += '''genome=%s\n''' % genome
    s += '''ncores=%d\n''' % n_threads
    s += '''\n'''

    if module == "combined":
        s += '''bams=%s\n''' % bam_fn
        s += '''barcodes=%s\n''' % barcode_fn
    s += '''count_mtx_dir=%s\n''' % mtx_dir
    s += '''\n'''

    s += '''cell_anno_fn=%s\n''' % cell_anno_fn
    s += '''ref_cell_type="%s"\n''' % ref_cell_type
    s += '''\n'''

    if module == "combined":
        s += '''gmap=%s\n''' % gmap_fn
        s += '''eagle=%s\n''' % eagle_fn
        s += '''snpvcf=%s\n''' % snp_vcf
        s += '''paneldir=%s\n''' % panel_dir
        s += '''\n'''

    s += '''# preprocess args\n'''
    if module == "combined":
        s += '''prefix=${sid}.numbat\n'''
    else:
        s += '''prefix=${sid}.numbat_rdr\n'''
    s += '''\n'''

    if module == "combined":
        s += '''platform_option=\n'''
        s += '''if [ "$platform" = "10x" ]; then\n'''
        s += '''    platform_option=""\n'''
        s += '''elif [ "$platform" = "smartseq" ]; then\n'''
        s += '''    platform_option="--smartseq"\n'''
        s += '''elif [ "$platform" = "bulk" ]; then\n'''
        s += '''    platform_option = "--bulk"\n'''
        s += '''else\n'''
        s += '''    echo "invalid platform '$platform'"\n'''
        s += '''    exit 1\n'''
        s += '''fi\n'''
        s += '''\n'''

        s += '''#Run SNP pileup and phasing with 1000G\n'''
        s += '''if [ ! -e "$out_dir/allele" ]; then\n'''
        s += '''    mkdir -p $out_dir/allele\n'''
        s += '''fi\n'''
        s += '''\n'''

        s += '''#usage: pileup_and_phase.R [-h] --label LABEL --samples SAMPLES --bams BAMS\n'''
        s += '''#                          [--barcodes BARCODES] --gmap GMAP [--eagle EAGLE]\n'''
        s += '''#                          --snpvcf SNPVCF --paneldir PANELDIR --outdir OUTDIR\n'''
        s += '''#                          --ncores NCORES [--UMItag UMITAG]\n'''
        s += '''#                          [--cellTAG CELLTAG] [--smartseq] [--bulk]\n'''
        s += '''#Arguments:\n'''
        s += '''#  -h, --help           show this help message and exit\n'''
        s += '''#  --label LABEL        Individual label\n'''
        s += '''#  --samples SAMPLES    Sample names, comma delimited\n'''
        s += '''#  --bams BAMS          BAM files, one per sample, comma delimited\n'''
        s += '''#  --barcodes BARCODES  Cell barcode files, one per sample, comma delimited\n'''
        s += '''#  --gmap GMAP          Path to genetic map provided by Eagle2 (e.g.\n'''
        s += '''#                       Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)\n'''
        s += '''#  --eagle EAGLE        Path to Eagle2 binary file\n'''
        s += '''#  --snpvcf SNPVCF      SNP VCF for pileup\n'''
        s += '''#  --paneldir PANELDIR  Directory to phasing reference panel (BCF files)\n'''
        s += '''#  --outdir OUTDIR      Output directory\n'''
        s += '''#  --ncores NCORES      Number of cores\n'''
        s += '''\n'''

        s += '''/usr/bin/time -v Rscript $work_dir/pileup_and_phase.R   \\\n'''
        s += '''    --label  $sid    \\\n'''
        s += '''    --samples  $sid  \\\n'''
        s += '''    --bams  $bams   \\\n'''
        s += '''    --barcodes  $barcodes   \\\n'''
        s += '''    --gmap  $gmap   \\\n'''
        s += '''    --eagle  $eagle  \\\n'''
        s += '''    --snpvcf  $snpvcf  \\\n'''
        s += '''    --paneldir  $paneldir  \\\n'''
        s += '''    --outdir  $out_dir/allele  \\\n'''
        s += '''    --ncores  $ncores          \\\n'''
        s += '''    $platform_option\n'''
        s += '''\n'''
        s += '''\n'''

        s += '''#Rscript $work_dir/numbat.R  \\\n'''
        s += '''#  <allele file>       \\\n'''
        s += '''#  <count matrix dir>  \\\n'''
        s += '''#  <cell anno file>   \\\n'''
        s += '''#  <ref cell type>    \\\n'''
        s += '''#  <out dir>         \\\n'''
        s += '''#  <out prefix>      \\\n'''
        s += '''#  <genome version>   \\\n'''
        s += '''#  <sequencing platform>  \\\n'''
        s += '''#  <ncores>\n'''
        s += '''\n'''

        s += '''/usr/bin/time -v Rscript $work_dir/numbat.R  \\\n'''
        s += '''    $out_dir/allele/${sid}_allele_counts.tsv.gz    \\\n'''
        s += '''    $count_mtx_dir      \\\n'''
        s += '''    $cell_anno_fn    \\\n'''
        s += '''    "$ref_cell_type"       \\\n'''
        s += '''    $out_dir/cnv    \\\n'''
        s += '''    $prefix         \\\n'''
        s += '''    $genome         \\\n'''
        s += '''    $platform       \\\n'''
        s += '''    $ncores\n'''

    else:
        s += '''#Rscript $work_dir/numbat.rdr.R  \\\n'''
        s += '''#  <count matrix dir>  \\\n'''
        s += '''#  <cell anno file>   \\\n'''
        s += '''#  <ref cell type>    \\\n'''
        s += '''#  <out dir>         \\\n'''
        s += '''#  <out prefix>      \\\n'''
        s += '''#  <genome version>   \\\n'''
        s += '''#  <ncores>\n'''
        s += '''\n'''

        s += '''/usr/bin/time -v Rscript $work_dir/numbat.rdr.R  \\\n'''
        s += '''    $count_mtx_dir      \\\n'''
        s += '''    $cell_anno_fn    \\\n'''
        s += '''    "$ref_cell_type"       \\\n'''
        s += '''    $out_dir    \\\n'''
        s += '''    $prefix         \\\n'''
        s += '''    $genome         \\\n'''
        s += '''    $ncores\n'''

    s += '''\n'''
    s += '''set +ux\n'''
    s += '''conda deactivate\n'''
    s += '''echo All Done!\n'''

    return(s)


if __name__ == "__main__":
    # demo
    work_dir_lst = [
        "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_numbat/raw_1x",
        "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_numbat/simu_1x"
    ]
    sid_lst = ["raw_1x", "simu_1x"]
    n = len(work_dir_lst)
    gen_qsub_batch(
        sid_lst = sid_lst,
        mtx_dir_lst = [os.path.join(d, "data/matrix") for d in work_dir_lst],
        cell_anno_fn_lst = [os.path.join(d, "data/matrix/cell_anno.tsv") for d in work_dir_lst],
        ref_cell_type_lst = ["immune cells"] * n,
        work_dir_lst = work_dir_lst,
        run_script = "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_numbat/run.sh",
        repo_dir = "/home/xianjie/projects/scCNVSim",
        out_dir_lst = None,
        seq_platform = "10x",      # "10x" or "smartseq"
        genome = "hg38",
        n_threads = 10,
        module = "RDR",
        bam_fn_lst = None, barcode_fn_lst = None,
        gmap_fn = None, eagle_fn = None,
        snp_vcf = None, panel_dir = None,
        conda_env = "SCSC",
        qsub_q = "cgsd", 
        qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
        verbose = True
    )