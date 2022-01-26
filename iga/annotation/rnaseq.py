"""
Isoseq relevant utils
"""
import logging
import os

from iga.apps.base import sh, emain, conda_act, abspath_list, get_prefix, bsub, waitjob


# 0: Reference (Hisat2-build)
# 1: reads1
# 2: reads2
# 3: threads
# Output: reads1.bam

def prepare():
    print("""
    
    python -m iga.annotation.rnaseq reads_align_assembly "elumb.lncRNA.EuG11_1.clean.fq.gz elumb.lncRNA.EuG11_2.clean.fq.gz" $REF --threads 20
    """)

# 0 ref
# 1 reads1
# 2 read2
# 3 thread
align_sh = r"""

# if [ ! -e {0}.1.ht2 ]
# then
#     echo "hisat2-index not exists building
#     hisat2-build {0} {0}
# fi

hisat2 --rna-strandness RF --mp 3,1 -p {3} -x {0} -1 {1} -2 {2} | samtools sort -@ {3} -o {1}.bam 
# && stringtie -p $threads -o $newID.with_novel.gtf  $newID.bam
"""

# 0: bam
# 1: CPU
# 2: prefix
trinity_sh = r"""
Trinity --SS_lib_type RF \
        --genome_guided_bam {0} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU {1} \
          --output {2}_trinity 
"""


def reads_align_assembly(ref=None, reads=None, threads=30, output=''):
    """
    Align RNA-Seq reads to reference and assemble with tirinity

    Make sure hisat2-index exists
    :param reads: "root_1.fq.gz root_2.fq.gz"
    :param ref: AA.genome
    :param threads: 30 default
    :return:
    """
    if ' ' in reads or '\t' in reads:
        (read1, read2) = reads.strip().split()
    else:
        logging.warning("Single End reads found!")
    if output == '':
        prefix = get_prefix(read1)
    else:
        prefix = output
    align_cmd = align_sh.format(ref, read1, read2, threads)
    trinity_cmd = trinity_sh.format(read1 + ".bam", threads, read1)
    if os.path.exists(read1 + ".bam"):
        cmd = trinity_cmd
    else:
        cmd = align_cmd + trinity_cmd
    jobid = bsub(cmd, cpus=threads, name="Trinity_Guided_read1")
    logging.info("Results is {}_trinity/Trinity-GG.fasta".format(prefix))
    waitjob(jobid)


if __name__ == "__main__":
    emain()
