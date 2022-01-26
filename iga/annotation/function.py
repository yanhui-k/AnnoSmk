"""
functional annotation: iprscan
"""

# 0 input pep
# 1 threads
from iga.apps.base import bsub, emain, waitjob
import os.path as op

# 0 pep file
# 1 threads
iprscan_sh = """
interproscan.sh -i {0} -b {0}.out -cpu {1} -f tsv -iprlookup --goterms -pa -dp"""


def iprscan(pep_file=None, threads=30):
    """
    interproscan wrapper
    :param pep_file:
    :return:
    """
    job_name = "ipr." + op.basename(pep_file)
    cmd = iprscan_sh.format(pep_file, threads)
    job = bsub(cmd, cpus=threads, name=job_name, queue='Q64C1T_X4')
    #waitjob(job)


if __name__ == "__main__":
    emain()
