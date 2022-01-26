"""
blast related wrappers
"""
from collections import defaultdict

from iga.apps.base import emain, bsub, waitjob, sh

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 ref
# 1 qry
# 2 output
blastp_sh = r"""
makeblastdb -in {0} -dbtype prot
blastp  -num_threads 5 -db {0} -query {1} -out {2}.bln -evalue 1e-5 -outfmt 7
"""


def blastp(ref=None, qry=None, threads=5):
    r"""
    :param query:
    :param ref:
    :return:
    """
    cmd = blastp_sh.format(ref, qry, qry)
    job = bsub(cmd, cpus=threads)
    waitjob(job)
    return 0


# 0 ref
# 1 qry
# 2 output
blastn_sh = r"""
makeblastdb -in {0} -dbtype nucl
blastn  -num_threads 5 -db {0} -query {1} -out {2}.bln -evalue 1e-5 -outfmt 7
"""


def blastn(ref=None, qry=None, threads=5):
    r"""
    :param qry:
    :param ref:
    :return:
    """
    cmd = blastn_sh.format(ref, qry, qry)
    job = bsub(cmd, name='blast', cpus=5)
    #waitjob(job)
    return 0


def blast2bed(bln=None):
    """
    %s bln > bln.bed
    :param bln:
    :return:
    """
    cmd = r"""
awk '{print $1"\t"$7"\t"$8"\t"$2"\t"$12"\t."}' {0} > {0}.bed
""".format(bln)
    sh(cmd)

def filter_reciprocal_best(bln=None):
    qry_best = {}
    qry_line = {}
    ref_best = {}
    highest = defaultdict(int)
    with open(bln) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            mylist = line.rstrip().split()
            qry = mylist[0]
            ref = mylist[1]
            try:
                bitscore = float(mylist[-1])
            except ValueError:
                logger.error(line)
                continue
            if (bitscore > highest[qry]):
                highest[qry] = bitscore
                qry_best[qry] = ref
                qry_line[qry] = line.rstrip()
            if (bitscore > highest[ref]):
                highest[ref] = bitscore
                ref_best[ref] = qry
    for k in qry_best:
        this_best_ref = qry_best[k]
        if (ref_best[this_best_ref] == k):
            # Reciprocal best
            print(qry_line[k])


if __name__ == "__main__":
    emain()
