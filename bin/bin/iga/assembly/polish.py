"""
Polish work scripts
"""

from iga.apps.base import conda_act, get_prefix, bsub, emain, abspath_list

import os.path as op

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 contig.fa [abs path]
# 1 sgs.fq.gzs [abs path]
# 2 prefix
# 3 threads recommand 30
nextpolish_sh = r"""
export PATH=/ds3200_1/users_root/yitingshuang/lh/anaconda2/bin:$PATH

export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextPolish/bin:$PATH
export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/NextPolish:$PATH

WORKDIR="workdir_nextpolish_"{2}
mkdir -p $WORKDIR
cd $WORKDIR
PREFIX={2}
contig={0}

touch sgs.fofn ; rm sgs.fofn
touch $PREFIX.0.fa ; rm $PREFIX.0.fa
     
ls {1} > sgs.fofn
ln -s $contig $PREFIX.0.fa

#for i in `seq 0 1`
for i in 0
do
    j=`expr $i + 1`
    GENOME=$PREFIX.$i.fa
    OUTPUT=$PREFIX.$j.fa
    WORKDIR=./03_polish.sgs_round$i
    echo "[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = {3}
multithread_jobs = 1
genome = $GENOME
genome_size = auto
workdir = $WORKDIR
polish_options = -p {{multithread_jobs}}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
">polish_sgs.$i.cfg
#Run
    nextPolish polish_sgs.$i.cfg
    cat $WORKDIR/03.kmer_count/05.polish.ref.sh.work/polish_genome*/*.fasta > $OUTPUT
done
"""


def nextpolish(contig=None, fastq=None, threads=30, queue='Q104C512G_X4'):
    r"""
    :param contig: to be polished contigs
    :param fastq: fastqs, if multiple ,input with quotes like "a.fq b.fq"
    :return:
    """
    fastq_list = fastq.split()
    abspath_list(fastq_list)
    contig = op.abspath(contig)
    logger.warning(contig)
    logger.warning(fastq_list)
    fastq = " ".join(fastq_list)
    prefix = get_prefix(contig)
    cmd = nextpolish_sh.format(contig, fastq, prefix, threads)
    bsub(cmd, name="nextpolish" + prefix, cpus=threads, queue=queue)


# 0 contig
# 1 bam
# 2 threads
gcpp_sh = """
set -euxo pipefail
pbmm2 align -j {2} --sort  {0} {1} {0}.bam
samtools index {0}.bam
nproc={2}
gcpp --algorithm=arrow -x 5 -X 120 -q 0 -j $nproc \
        -r {0} {0}.bam \
        -o {0}.polish.fasta,{0}.polish.fastq,{0}.polish.vcf
rm {0}.bam
"""


def gcpp(contig=None, bam=None, threads=100):
    r"""
    Polish Pacbio Assembly with Pacbio Reads
    %s contig bam
    :param contig: Raw assembly
    :param bam: subreads bam
    :param threads: default=100"
    :return:
    """
    # assembly, subreads
    prefix = get_prefix(contig)
    cmd = conda_act.format('falcon') + gcpp_sh.format(contig, bam, threads)
    #    cmd = conda_act
    bsub(cmd, cpus=threads, name="gcpp" + prefix)


# 0 contig (abs path)
# 1 subreads.fasta (abs path)
# 2 threads
# 3 prefix
purge_dups_sh=r"""
pri_asm={0}

mkdir -p workdir_purge_dups_{3}

minimap2 -t {2} -xmap-pb $pri_asm {1} | gzip -c - > align.paf.gz

pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs 2>calcults.log
#        Notice If you have a large genome, please set minimap2 -I option to ensure the genome can be indexed once, otherwise read depth can be wrong.
#        Step 1. Split an assembly and do a self-self alignment. Commands are following:
split_fa $pri_asm > $pri_asm.split
minimap2 -t {2} -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

#Step 2. Purge haplotigs and overlaps with the following command.
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

#Step 3. Get purged primary and haplotig sequences from draft assembly.
get_seqs dups.bed $pri_asm
"""

def purge_dups(contig=None, fasta=None, threads=40, queue='Q104C512G_X4'):
    prefix = get_prefix(contig)
    contig = op.abspath(contig)
    fasta = op.abspath(fasta)
    cmd = purge_dups_sh.format(contig, fasta, threads, prefix)
    bsub(cmd, cpus=threads, name="purgedup"+prefix, queue=queue)

# def main():
#     prog_name = "Polish Pacbio Assembly with Pacbio Reads"
#     usage = "Polish Pacbio Assembly with Pacbio Reads"
#
#     parser = argparse.ArgumentParser(
#         prog=prog_name,
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description=textwrap.dedent(usage),
#         epilog="")
#     parser.add_argument("CONTIG", help="Raw assembly")
#     parser.add_argument("BAM", help="subreads bam")
#     parser.add_argument("-t", "--threads", default=100, type=int, help="flanking distance default (1000)")
#     args = parser.parse_args()
#
#     ctg_file = args.CONTIG
#     bam_file = args.BAM
#     gcpp(ctg_file, bam_file)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()
