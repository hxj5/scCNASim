# gen_infercnv.py - generate infercnv scripts


import os
import stat


def assert_e(path):
    assert os.path.exists(path)


def gen_qsub_batch(
    sid_lst,
    mtx_dir_lst,
    cell_anno_fn_lst,
    ref_cell_type_lst,
    gene_anno_fn_lst,
    work_dir_lst,
    run_script,
    repo_dir,
    out_dir_lst = None,        # trick: ["$work_dir/infercnv"] * n
    seq_platform = "10x",      # "10x" or "smartseq"
    n_threads = 10,
    conda_env = "XCLBM",
    qsub_q = "cgsd", 
    qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
    verbose = True
):
    n = len(sid_lst)
    assert len(mtx_dir_lst) == n
    assert len(cell_anno_fn_lst) == n
    assert len(ref_cell_type_lst) == n
    assert len(gene_anno_fn_lst) == n
    assert len(work_dir_lst) == n
    if out_dir_lst is None:
        out_dir_lst = [None] * n
    else:
        assert len(out_dir_lst) == n

    if verbose:
        print("write to qsub scripts ...")

    qsub_script_lst = []
    for sid, mtx_dir, cell_anno_fn, ref_cell_type, gene_anno_fn, work_dir, \
        out_dir in \
        zip(sid_lst, mtx_dir_lst, cell_anno_fn_lst, ref_cell_type_lst, \
            gene_anno_fn_lst, work_dir_lst, out_dir_lst
        ):
        qsub_script_fname = "%s.infercnv.qsub.sh" % sid
        qsub_script_fpath = os.path.join(work_dir, qsub_script_fname)
        s = gen_qsub(
            sid = sid,
            mtx_dir = mtx_dir,
            cell_anno_fn = cell_anno_fn,
            ref_cell_type = ref_cell_type,
            gene_anno_fn = gene_anno_fn,
            repo_dir = repo_dir,
            out_dir = out_dir,
            seq_platform = seq_platform,
            n_threads = n_threads,
            conda_env = conda_env,
            qsub_N = sid + "_infercnv", 
            qsub_q = qsub_q, 
            qsub_l_nodes = qsub_l_nodes, qsub_l_ppn = qsub_l_ppn, 
            qsub_l_mem = qsub_l_mem, qsub_l_walltime = qsub_l_walltime,
            qsub_o = "%s_infercnv.out" % sid, 
            qsub_e = "%s_infercnv.err" % sid
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
    gene_anno_fn,
    repo_dir,
    out_dir = None,
    seq_platform = "10x",      # "10x" or "smartseq"
    n_threads = 10,
    conda_env = "XCLBM",
    qsub_N = "infercnv", qsub_q = "cgsd", 
    qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
    qsub_o = "infercnv.out", qsub_e = "infercnv.err"
):
    # check args
    assert_e(mtx_dir)
    assert_e(cell_anno_fn)
    assert_e(gene_anno_fn)
    assert_e(repo_dir)
    assert seq_platform in ("10x", "smartseq")

    assert_e("/bin/bash")
    #assert_e("/home/.bashrc")
    assert_e(os.path.join(repo_dir, "scripts/cnv_calling/infercnv/infercnv.R"))
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

    s += '''cp  $repo_dir/scripts/cnv_calling/infercnv/infercnv.R  $work_dir\n'''
    s += '''\n'''

    s += '''#Rscript $work_dir/infercnv.R \\\n'''
    s += '''#  <sample id>       \\\n'''
    s += '''#  <matrix dir>      \\\n'''
    s += '''#  <cell anno file>     \\\n'''
    s += '''#  <ref cell type>      \\\n'''
    s += '''#  <gene anno file>     \\\n'''
    s += '''#  <out dir>            \\\n'''
    s += '''#  <sequencing platform>: "10x" or "smartseq"   \\\n'''
    s += '''#  <number of threads>\n'''
    s += '''\n'''

    s += '''/usr/bin/time -v Rscript $work_dir/infercnv.R \\\n'''
    s += '''    "%s"  \\\n''' % sid
    s += '''    %s  \\\n''' % mtx_dir
    s += '''    %s  \\\n''' % cell_anno_fn
    s += '''    "%s"     \\\n''' % ref_cell_type
    s += '''    %s  \\\n''' % gene_anno_fn
    s += '''    $out_dir        \\\n'''
    s += '''    "%s"           \\\n''' % seq_platform
    s += '''    %d\n''' % n_threads
    s += '''\n'''
    s += '''set +ux\n'''
    s += '''conda deactivate\n'''
    s += '''echo All Done!\n'''

    return(s)


if __name__ == "__main__":
    # demo
    work_dir_lst = [
        "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_infercnv/raw_1x",
        "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_infercnv/simu_1x"
    ]
    sid_lst = ["raw_1x", "simu_1x"]
    n = len(work_dir_lst)
    gen_qsub_batch(
        sid_lst = sid_lst,
        mtx_dir_lst = [os.path.join(d, "data/matrix") for d in work_dir_lst],
        cell_anno_fn_lst = [os.path.join(d, "data/matrix/cell_anno.tsv") for d in work_dir_lst],
        ref_cell_type_lst = ["immune cells"] * n,
        gene_anno_fn_lst = ["/groups/cgsd/xianjie/data/refapp/infercnv/hg38_gene_note_noheader_unique.txt"] * n,
        work_dir_lst = work_dir_lst,
        run_script = "/groups/cgsd/xianjie/debug/test-sccnvsim/test_gen_infercnv/run.sh",
        repo_dir = "/home/xianjie/projects/scCNVSim",
        out_dir_lst = None,
        seq_platform = "10x",      # "10x" or "smartseq"
        n_threads = 10,
        conda_env = "XCLBM",
        qsub_q = "cgsd", 
        qsub_l_nodes = 1, qsub_l_ppn = 5, qsub_l_mem = 100, qsub_l_walltime = 100,
        verbose = True
    )
