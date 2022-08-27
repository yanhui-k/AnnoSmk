"""
Isoseq relevant utils
"""
import os

from iga.apps.base import sh, emain, conda_act, abspath_list, get_prefix, bsub

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 subreads.bam
# 1 workdir
# 2 output.bam
# def sam2gff(sam, gff=""):
#
# isoseq_official_sh=r"""
# ccs [movie].subreads.bam [movie].ccs.bam --min-rq 0.9
# lima --isoseq --dump-clips --no-pbi --peek-guess -j 24 ccs.bam primers.fasta demux.bam
# isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta flnc.bam
# bamtools convert -format fastq -in flnc.bam > flnc.fastq
# """

# 0 subreads.bam
# 1 workdir
isoseq_sh = r"""export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/isoseq3/BGI-Full-Length-RNA-Analysis-Pipeline/bin:$PATH
export PERL5LIB=""
WORKDIR={1}
mkdir ${{WORKDIR}}
ccs {0} ${{WORKDIR}}/ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.9
cd ${{WORKDIR}}
samtools view ccs.bam | awk '{{print ">"$1"\n"$10}}' > ccs.fa
echo ">primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT">primer.fa
makeblastdb -in primer.fa -dbtype nucl
blastn -num_threads 32 -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7
classify_by_primer -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 16 -min_isolen 200 -outdir ./
flnc2sam ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam
samtools view -bS isoseq_flnc.sam > isoseq_flnc.bam
isoseq3 cluster isoseq_flnc.bam unpolished.bam --verbose --use-qvs
#isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 10
# pbindex subreads.bam
# for i in `seq 0 9`
# do
#     isoseq3 polish unpolished.${{i}}.bam ${{ROOT}}/input/*.subreads.bam polished.${{i}}.bam --verbose &
# done
# wait
# samtools merge -@ 20 polished_total.bam polished.*.bam
# isoseq3 summarize polished_total.bam summary.csv
# Result is polished_total.bam.fastq 
"""


def isoseq_bgi(subreads=None, workdir=''):
    r"""
    isoseq subreads.fasta

    Wrapper for `isoseq`
    """
    if (type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)


# 0 workdir
# 1 subreads.bam
# 2 primer.fa
# 3 threads
isoseq_pb_sh = r"""mkdir -p {0}
cd {0}
ln -s ../{1}
if [ ! -e {2} ]
then
    ln -s ../{2}
else
    echo "{2} already exist, will not try re-link"
fi

ROOT=$PWD
INPUTBAM={1}
PRIMER={2}
THREADS={3}
ccs ${{INPUTBAM}} ${{INPUTBAM}}.ccs.bam --min-rq 0.9 -j $THREADS
        #--min-rq                  FLOAT  Minimum predicted accuracy in [0, 1]. [0.99]
        #-j,--num-threads          INT    Number of threads to use, 0 means autodetection. [0]
    #lima is used to remove primer sequence
    #but can it be used to identify reads containing primer sequence as full length reads?
lima ${{INPUTBAM}}.ccs.bam ${{PRIMER}} ${{INPUTBAM}}.fl.bam --isoseq  --peek-guess
        #--isoseq                              Activate specialized IsoSeq mode.
        #--peek-guess                          Try to infer the used barcodes subset, by peeking at the first 50,000 ZMWs,
    #                                     whitelisting barcode pairs with more than 10 counts and mean score <A1><DD> 45.
isoseq3 refine ${{INPUTBAM%.bam}}.fl.*.bam ${{PRIMER}} ${{INPUTBAM%.bam}}.flnc.bam --require-polya -j $THREADS
        #refine     Remove polyA and concatemers from FL reads and generate FLNC transcripts (FL to FLNC)
        #--min-polya-length  INT   Minimum poly(A) tail length. [20]
        #--require-polya           Require FL reads to have a poly(A) tail and remove it.
        #-j,--num-threads    INT   Number of threads to use, 0 means autodetection. [0]

isoseq3 cluster ${{INPUTBAM%.bam}}.flnc.bam ${{INPUTBAM%.bam}}.clustered.bam --verbose --use-qvs -j $THREADS

"""


def isoseq_pb(subreads=None, primer=None, workdir='', threads=50):
    r"""
    convert Isoseq(pacbio standard) subreads.bam to flnc.fastq
    :param subreads: Multiple bam inputs not supported currently, a TODO
    :param primer: Primer fasta
    :param workdir: if not given, default is base name of first subreads
    :param threads:  threads
    :return:
    """
    # subreads = subreads.replace(' ', ',')
    # subreads = subreads.split(',')

    prefix = get_prefix(subreads)

    # abspath_list(subreads)

    # subreads = ','.join(subreads)

    if (workdir == ''):
        workdir = "workdir_isoseq_" + prefix
    # logging.debug(workdir)
    cmd = conda_act.format('isoseq3') + isoseq_pb_sh.format(workdir, subreads, primer, threads)
    bsub(cmd, name="isoseq3" + prefix, cpus=threads)


# 0 input

isoseq3_cluster_sh = """
isoseq3 cluster {0} {1} --verbose --use-qvs -j {2}
"""

# def isoseq3_cluster(BAM=None, threads=20):
#
#     sam_merge = "samtools merge {}"
#     cmd = conda_act.format('isoseq3') + iso


if __name__ == "__main__":
    emain()
