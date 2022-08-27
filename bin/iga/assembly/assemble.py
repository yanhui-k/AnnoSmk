"""
Various assembler wrapper
"""
import os

import os.path as op
import time

from iga.apps.base import conda_act, Config, mkdir, get_prefix, sh, bsub, emain, abspath_list, waitjob, mv

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


# 0 reads.fasta
# 1 prefix, like altr
# 2 threads
hifiasm_sh='''
READS={0}
PREFIX={1}
THREADS={2}

mkdir -p workdir_hifiasm_{1}

cd workdir_hifiasm_{1}

hifiasm -o ${{PREFIX}} -t${{THREADS}}  ${{READS}} 2> ${{PREFIX}}.log
awk '/^S/{{print ">"$2;print $3}}' ${{PREFIX}}.p_ctg.gfa > ${{PREFIX}}.p_ctg.fa  # get primary contigs in FASTA
'''


def hifiasm(ccs_reads=None, threads=64, prefix='', submit='T', queue='Q104C512G_X4'):
    r"""
    flye runs on single machine
    :param ccs_reads:
    :param threads: default is 64, dependes how many available of host machine
    :param prefix: (species name)
    :param submit: T stands for submit this job to lsf, other value indicate output shell script but do not submit
    :param queue: Default is Q104C512G_X4, could also be Q64C1T_X4
    :return:
    """
    ccs_reads = op.abspath(ccs_reads)
    if prefix == '':
        prefix = get_prefix(ccs_reads)
    cmd_sh = hifiasm_sh.format(ccs_reads, prefix, threads)
    bsub(cmd_sh, queue=queue, name='hifiasm.' + prefix, submit=submit, cpus=threads)


# 0: prefix
# 1: corrected reads.fasta
# 2: genome size, 123m or 1g
# 3: threads
flye_sh = '''
INPUT={1}
WORKDIR=workdir_flye_{0}
flye --pacbio-corr ${{INPUT}} --out-dir ${{WORKDIR}} --genome-size {2} --threads {3}
'''


def flye(corrected_reads=None, genome_size=None, threads=64, prefix='', submit='T', queue='Q104C512G_X4'):
    r"""
    flye runs on single machine
    :param corrected_reads: corrected pacbio reads
    :param genome_size: 100m stands for 100 Mb, 1gb is also supported
    :param threads: default is 64, dependes how many available of host machine
    :param prefix: (species name)
    :param submit: T stands for submit this job to lsf, other value indicate output shell script but do not submit
    :param queue: Default is Q104C512G_X4, could also be Q64C1T_X4
    :return:
    """
    corrected_reads = op.abspath(corrected_reads)

    logger.debug(corrected_reads)

    if prefix == '':
        prefix = get_prefix(corrected_reads)

    cmd_sh = conda_act.format('flye')
    cmd_sh += flye_sh.format(prefix, corrected_reads, genome_size, threads)

    bsub(cmd_sh, queue=queue, name='flye.' + prefix, submit=submit, cpus=threads)


# 0: prefix
# 1: corrected reads.fasta
# 2: genome size, 123m or 1g
# 3: threads
wtdbg_sh = """
PREFIX={0}
LONGREADS={1}
WORKDIR=workdir_wtdbg_${{PREFIX}}
# 12m or 1g
GENOMESIZE={2}
THREADS={3}

#rs sq ont ccs
LIB=ccs

echo -n "Start at "
date
#Pipeline
wtdbg2.pl -t ${{THREADS}} -x ${{LIB}} -g ${{GENOMESIZE}} -o ${{PREFIX}} ${{LONGREADS}} >assemble.log 2>assemble.err

echo -n "End at "
date
"""


def wtdbg(corrected_reads=None, genome_size=None, threads=64, prefix='', submit='T', queue='Q104C512G_X4'):
    r"""
    :param corrected_reads: corrected pacbio reads
    :param genome_size: 100m stands for 100 Mb, 1gb is also supported
    :param threads: default is 64, dependes how many available of host machine
    :param prefix: (species name)
    :param submit: T stands for submit this job to lsf, other value indicate output shell script but do not submit
    :param queue: Default is Q104C512G_X4, could also be Q64C1T_X4
    :return:
    """
    corrected_reads = op.abspath(corrected_reads)

    logger.debug(corrected_reads)

    if prefix == '':
        prefix = get_prefix(corrected_reads)

    workdir = 'workdir_wtdbg_{}'.format(prefix)
    if op.exists(workdir):
        mv(workdir, workdir + str(time.time()).replace('.', ''))
    if not mkdir(workdir):
        logger.error("Workdir existing, exiting...")
        exit(1)

    os.chdir(workdir)

    cmd_sh = wtdbg_sh.format(prefix, corrected_reads, genome_size, threads)

    bsub(cmd_sh, queue=queue, name='wtdbg.' + prefix, submit=submit, cpus=threads)


# threads.config
canu_threads_config = """
maxMemory=500g
maxThreads=64

merylThreads=64
merylMemory=200

cormhapThreads=8
cormhapConcurrency=1
obtmhapThreads=8
utgmhapThreads=8

cnsThreads=8
corThreads=8

corovlThreads=8
obtovlThreads=8
utgovlThreads=8
"""

# 0: prefix, usually species name, will be used as directory name
# 1: full length to fasta.gz
# 2: genome size
# 3: type, could be pacbio, nanopore, and pacbio-hifi
canu_sh = r'''
# CANU=/ds3200_1/users_root/yitingshuang/lh/bin/canu/canu-2.0/Linux-amd64/bin/canu
CANU=/ds3200_1/users_root/yitingshuang/lh/bin/canu-2.1.1/bin/canu
PREFIX={0}
WORKDIR="workdir_canu_"${{PREFIX}}_{2}
INPUT={1}
#pacbio or nanopore
TYPE={3}

mkdir -p ${{WORKDIR}}

#polyploid_param
# corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \

${{CANU}} \
 -d ${{WORKDIR}} -p ${{PREFIX}}  \
 -s threads.config \
 correctedErrorRate=0.035 \
 utgOvlErrorRate=0.065 \
 trimReadsCoverage=2 \
 trimReadsOverlap=500 \
genomeSize={2} \
useGrid=true \
minReadLength=1000 \
minOverlapLength=600 \
-${{TYPE}} \
${{INPUT}} \
  >>${{WORKDIR}}/assemble.log 2>>${{WORKDIR}}/assemble.err
'''


def canu(subreads=None, genome_size=None, prefix='', type='pacbio', etc='', submit='T', queue='Q104C512G_X4'):
    r"""
    :param subreads: pacbio DNA subreads (FASTA or FASTQ, can be uncompressed, gzip, bzip2 or xz compressed)
    :param genome_size: 100m stands for 100 Mb
    :param prefix: (species name)
    :param etc: Deprecated for now
    :param submit: T stands for submit this job to lsf, other value indicate output shell script but do not submit
    :param type: 'pacbio|nanopore|pacbio-hifi'
    :param queue: Default is Q104C512G_X4, could also be Q64C1T_X4
    :return:
    """
    logger.debug(subreads)

    if '.bam' in subreads:
        subreads = bam2fastq(subreads)

    subreads = op.abspath(subreads)

    logger.debug(subreads)

    if prefix == '':
        prefix = get_prefix(subreads)

    with open('threads.config', 'w') as fh:
        fh.write(canu_threads_config)

    # workdir = 'workdir_canu_{}'.format(prefix)
    # if op.exists(workdir):
    #     mv(workdir, workdir + str(time.time()).replace('.', ''))
    # if not mkdir(workdir):
    #     logger.error("Workdir existing, exiting...")
    #     exit(1)

    cmd_sh = canu_sh.format(prefix, subreads, genome_size, type)

    # change lib_file
    lib_file = '/ds3200_1/users_root/yitingshuang/lh/bin/canu/canu-2.0/Linux-amd64/lib/site_perl/canu/Grid_LSF.pm'
    original = r'setGlobalIfUndef("gridEngineSubmitCommand",.*'
    target = r'setGlobalIfUndef("gridEngineSubmitCommand", "bsub -q {}");'.format(queue)
    sh("""sed -i 's/{0}/{1}/' {2}""".format(original, target, lib_file))

    bsub(cmd_sh, queue=queue, name='canu.' + prefix, submit=submit)


def bam2fastq(subreads=None, submit='T'):
    """
    Input zeins.bam zeins.bam.pbi
    Output zeins.fasta.gz
    :param subreads:
    :param submit: [T/F] whether to use bsub to submit
    :return:
    """
    if type(subreads) == list:
        abspath_list(subreads)
        subreads = " ".join(subreads)
    else:
        subreads = op.abspath(subreads)

    if not op.exists(subreads + '.pbi'):
        logger.error("{}.pbi is needed! Exiting..".format(subreads))
        return 1

    prefix = get_prefix(subreads)
    cmd = 'bam2fasta -o {0} {1}'.format(prefix, subreads)
    if submit == 'T':
        job = bsub(cmd, name='bam2fasta')
        waitjob(job)
    else:
        sh(cmd)
    return '{0}.fasta.gz'.format(prefix)


# 0 cfg_file
falcon_sh = r"""
fc_run {} >>falcon_run.out 2>>falcon_run.err
"""


def falcon(subreads=None, genome_size=None, prefix='', etc='', submit='T'):
    r"""
    :param subreads: pacbio DNA subreads (bam accepted)
    :param genome_size: (genomesize in bp)
    :param prefix: (species name)
    :param etc: (other fields need to be updated in falcon cfg)
    :param submit: T stands for submit this job to lsf, other value indicate output shell script but do not submit
    :return:
    """
    logger.debug(subreads)

    if '.bam' in subreads:
        subreads = bam2fastq(subreads)

    # Change Mb to bp
    if 'M' or 'm' in genome_size:
        genome_size = genome_size.replace('M', '')
        genome_size = genome_size.replace('m', '')
        genome_size = genome_size + '0'*6

    # Change Gb to bp
    if 'G' or 'g' in genome_size:
        genome_size = genome_size.replace('g', '')
        genome_size = genome_size.replace('G', '')
        genome_size = genome_size + '0'*9

    logger.debug(subreads)

    if type(subreads) == list:
        abspath_list(subreads)
        subreads = " ".join(subreads)
    else:
        subreads = op.abspath(subreads)

    if prefix == '':
        prefix = get_prefix(subreads)

    workdir = 'workdir_falcon_{}'.format(prefix)
    if op.exists(workdir):
        mv(workdir, workdir + str(time.time()).replace('.', ''))
    if not mkdir(workdir):
        logger.error("Workdir existing, exiting...")
        exit(1)
    os.chdir(workdir)

    fofn_file = 'bam.fofn'
    cfg_file = "falcon.cfg"
    fofn_file = op.abspath(fofn_file)
    sh("ls {0} > {1}".format(subreads, fofn_file))
    cfg = Config('falcon')
    cfg.update('[General]input_fofn={0}'.format(fofn_file))
    cfg.update('[General]genome_size={0}'.format(genome_size))
    if etc != '':
        cfg.update(etc)
    # logging.debug(cfg.get_text())
    cfg.write_to_file(cfg_file)

    cmd = conda_act.format('falcon') + falcon_sh.format(cfg_file)
    bsub(cmd, name='falcon' + prefix, submit=submit)


# 0 reference contig (abs path)
# 1 subreads.fasta (abs path)
purge_dups_sh = r"""
pri_asm={0}
subreads={1}


minimap2 -t 40 -xmap-pb $pri_asm $subreads | gzip -c - > align.paf.gz

pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs 2>calcults.log
#        Notice If you have a large genome, please set minimap2 -I option to ensure the genome can be 
# indexed once, otherwise read depth can be wrong.

#        Step 1. Split an assembly and do a self-self alignment. Commands are following:
split_fa $pri_asm > $pri_asm.split
minimap2 -t 40 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

#Step 2. Purge haplotigs and overlaps with the following command.
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

#Step 3. Get purged primary and haplotig sequences from draft assembly.
get_seqs dups.bed $pri_asm
"""


def purge_dups(contig=None, subreads=None, prefix=''):
    r"""
    purge dups wrapper
    :param contig:  contig fasta
    :param subreads: subreads.fasta[.gz]
    :param prefix: usually name of contig.fasta
    :return:
    """
    contig = op.abspath(contig)
    subreads = op.abspath(subreads)
    if prefix == '':
        prefix = get_prefix(subreads)

    workdir = 'workdir_purge_dups_{}'.format(prefix)
    mkdir(workdir)
    os.chdir(workdir)

    cmd = purge_dups_sh.format(contig, subreads)
    bsub(cmd, name="purge_dups_" + prefix, cpus=4)
    return 0


#0 fastq1
#1 fastq2
fastp_sh = """
fastp -i {0} -I {1} -o {0}.clean.fq.gz -O {1}.clean.fq.gz
"""


# 0 fastq[s]
# 1 threads
# 2 prefix
platanus_sh = """
platanus assemble -f {0} -t {1} -o {2}_assembly.fa -m {3}
"""


def platanus(fastq=None, clean='F', threads=20, mem=400):
    r"""
    assemble with platanus
    :param fastq:
    :param clean: [T/F] if T, use fastp to clean
    :param mem: (Gb)
    :return:
    """
    if type(fastq) == list:
        fastq = " ".join(fastq)
    prefix = get_prefix(fastq.split()[0])
    cmd = platanus_sh.format(fastq, threads, prefix, mem)
    bsub(cmd, name='platanus', cpus=threads)


if __name__ == '__main__':
    emain()
