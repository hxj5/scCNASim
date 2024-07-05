# main.py - cmdline interface.


import time

from logging import info, error


def main():
    pass


def main_core(conf):
    ret = -1

    # preprocessing


    # allele-specific feature counting
    conf.afc.out_dir = conf.g.out_dir

    


def main_run(conf):
    ret = -1

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        ret = main_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return(ret)
