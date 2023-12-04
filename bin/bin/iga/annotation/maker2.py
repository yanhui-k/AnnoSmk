"""
maker relevant utils
"""
import os
import os.path as op
import re
import sys
import time
from collections import defaultdict

from parse import parse

from iga.annotation.gff import GFF
from iga.apps.base import sh, conda_act, Config, abspath_list, split_fasta, mkdir, \
    mv, waitjob, bsub, emain, get_prefix

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

# 0 repeatmask.gff
prep_repeat_sh = r"""
GFF={0}
PREFIX=${{GFF%.gff}}
grep -v -e "Satellite" -e ")n" -e "-rich" ${{GFF}} \
    > ${{PREFIX}}.complex.gff3
# reformat to work with MAKER
cat ${{PREFIX}}.complex.gff3 | \
    perl -ane '$id; if(!/^\#/){{@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; 
    $_ = join("\t", @F)."\n"}}; print $_'  > ${{PREFIX}}.complex.reformat.gff3
/ds3200_1/users_root/yitingshuang/lh/bin/bundle/MAKER/rename_repeat_masker_gff_for_maker.pl \
${{PREFIX}}.complex.reformat.gff3 > ${{PREFIX}}.complex.reformat.MAKER.gff3
"""


def prep_repeat_gff(GFF=None):
    """
    :param GFF: repeat mask gff
    :return: Output
    """
    cmd = ''
    if re.search(r'.out$', GFF):
        cmd = "rmOutToGFF3.pl {0} > {0}.gff\n".format(GFF)
    cmd += prep_repeat_sh.format(GFF + ".gff")

    bsub(cmd, name="format_repeat_mask")
    return 0


# 0 genome.fasta
# 1 protein.fasta
# 2 prefix
prep_genblast_sh = """
REF={0}
QRY={1}
PREFIX={2}

ln -s ${{REF}}
# ln -s ${{REF}}.nin
# ln -s ${{REF}}.nhr
# ln -s ${{REF}}.nsq
ln -s ${{QRY}}
REF=`basename ${{REF}}`
QRY=`basename ${{QRY}}`

#genblast -p genblastg -q $QRY -t $REF -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
# -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o $PREFIX.genblast
# change to -h 1 if the above command failed

python -m iga.annotation.genblast2 run $QRY $REF $PREFIX --rank {3}

python -m iga.annotation.genblast2 filter_genblast $PREFIX.genblast*.gff > $PREFIX.slim.genblast.gff

python -m iga.annotation.genblast2 filter_early_stop $PREFIX.genblast*.pro > $PREFIX.genblast.noearly_stop.id

selectGFF.pl $PREFIX.genblast.noearly_stop.id $PREFIX.slim.genblast.gff > $PREFIX.filter.genblast.gff

sed 's/transcript/protein_match/; s/coding_exon/match_part/' $PREFIX.filter.genblast.gff > $PREFIX.final.gff
"""


def prep_genblast(genome=None, protein=None, rank=3, chunk=100, output=''):
    """
    Run genblast from protein to genome, output maker compatible gffs as well as normal gffs
    Relative path
    :param genome: genome fasta
    :param protein: protein fasta
    :param chunk: how many files the genome is splitted\
    :param output: the output gff file
    :return:
    """
    # test ok in /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/input/genblast/falcon_genblast/test
    abs_ref = op.abspath(genome)
    rel_pt = op.relpath(protein)
    workdir = genome + ".genblast." + rel_pt
    if op.exists(workdir):
        mv(workdir, workdir + str(time.time()).replace('.', ''))
    fasta_list = split_fasta(protein, workdir, chunk)
    job_list = []
    # if not (op.exists(abs_ref + ".nin") and op.exists(abs_ref + ".nsq") and op.exists(abs_ref + ".nhr")):
    #     sh('formatdb -p F -i {}'.format(abs_ref))
    workdir = op.abspath(workdir)
    for protein_i in fasta_list:
        os.chdir(workdir)
        mkdir(protein_i + '.run')
        # mv(protein_i, protein_i + '.run/' + protein_i)
        os.chdir(protein_i + '.run')
        final_prefix = protein_i
        cmd = prep_genblast_sh.format(abs_ref, "../" + protein_i, final_prefix, rank)
        job_list.append(bsub(cmd, name='genblast'))
    waitjob(job_list)
    if output == '':
        output = ".".join([abs_ref, rel_pt, '.gff'])
    sh('cat {}/*.run/*.final.gff > {}'.format(workdir, output))
    logging.debug("The resulting gff is {}".format(output))
    return 0


# collect_genblast_gff_sh = """
# PREFIX=${{{}%.gff}}
#
# python -m iga.annotation.genblast filter_genblast $PREFIX.gff > $PREFIX.slim.genblast.gff
#
# python -m iga.annotation.genblast filter_early_stop $PREFIX.pro > $PREFIX.genblast.noearly_stop.id
#
# selectGFF.pl $PREFIX.genblast.noearly_stop.id $PREFIX.slim.genblast.gff > $PREFIX.filter.genblast.gff
#
# sed 's/transcript/protein_match/; s/coding_exon/match_part/' $PREFIX.filter.genblast.gff > $PREFIX.final.gff
# """
#
#
# def collect_genblast_gff(gff=None):
#     sh(process_genblast_gff_sh.format(gff))


def prep_fastx2gff(fastx=None, genome=None, output='', workdir=''):
    r"""
    Prepare est evidence to maker compatible format
    :param fastx: fasta or fastq
    :param genome:
    :param output:
    :param workdir:
    :return:
    """
    fastq2gff(fastx, genome, output, workdir)


# 0 ref genome
# 1 qry fasta
fastq2gff_sh = """
minimap2 -t20 -C5 -ax splice {0} {1} |samtools view -F 256 -b >{1}.bam
bedtools bamtobed -split -i {1}.bam > {1}.raw.bed 
awk '$3-$2 > 1' {1}.raw.bed  > {1}.bed 
gt bed_to_gff3 {1}.bed | sort -k9,9 -k1,1 -k7,7 -k4,4n  > {1}.rawgff
"""


def fastq2gff(fastq=None, genome=None, output='', workdir=''):
    r"""
    align fastqs generated by Trinity or Isoseq to
    references and transform to maker acceptable gffs
    eg: %s fastq genome
    """
    # fastq = fastq)
    # genome = str(genome)
    logging.debug(fastq)
    logging.debug(genome)
    #    if (workdir == ''):
    #        workdir = "workdir_fastq2gff_" + fastq
    # workdir_sh.format(workdir) +
    cmd = conda_act.format('EDTA') + \
          fastq2gff_sh.format(genome, fastq)
    job_id = bsub(cmd, name='fastq2gff')
    waitjob(job_id)
    rawgff = fastq + '.rawgff'
    gff = format_gt_gff_to_maker_gff(rawgff)
    if output != '':
        mv(gff, output)
        logging.debug(output)
    else:
        logging.debug("Output file is {}".format(gff))
    return gff


def format_gt_gff_to_maker_gff(gff=None, max_intron_size=20000):
    """
    Format gff generated by gt to maker acceptable format
    Input Example:
    000034F|arrow_np1212    .       BED_feature     1899770 1900081 60      -       .       Name=TRINITY_GG_10000_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2360009 2360227 60      -       .       Name=TRINITY_GG_10001_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2359766 2359901 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360103 2360189 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360284 2360402 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    Output Example:
    000034F|arrow_np1212    match   Trinity_Minimap 1899770 1900081 0       -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 1899770 1900081 60      -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0-exon-1;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2360009 2360227 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 2360009 2360227 60      -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0-exon-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2359766 2361129 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0
    """
    qry1_file = gff
    output_file = gff + ".gff"
    output_buff = ""
    prefix = re.sub(r'\..*', '', qry1_file)
    intron_cutoff = int(max_intron_size)
    source = 'Trinity_Minimap'
    last_chr = ''
    last_end = ''
    last_strand = ''
    last_pos = []
    last_lines = []
    last_name = ''
    count = 0
    feat_count = defaultdict(int)
    with open(qry1_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            mylist = line.rstrip().split()
            # print(mylist[-1])
            mylist[-1] = mylist[-1].replace('/', '_')
            # Name feats
            this_chr = mylist[0]
            this_type = 'match_part'
            this_start = mylist[3]
            this_end = mylist[4]
            this_score = mylist[5]
            this_strand = mylist[6]
            this_phase = mylist[7]
            feat_parse = parse("Name={name}", mylist[-1])
            try:
                feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(
                    feat_count[feat_parse['name']])  # + this_chr.replace('|', '') + this_start
            except TypeError:
                logger.warning("TypeError on feat {}".format(mylist[-1]))
            this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
            line = "\t".join(
                [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])

            if (last_chr == '' or
                    last_name == feat_name and
                    last_chr == this_chr and
                    last_strand == this_strand and
                    abs(last_end - int(this_start)) < intron_cutoff):
                count += 1

                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_chr = this_chr
                last_end = int(this_end)
                last_pos.append(int(this_start))
                last_pos.append(int(this_end))
                last_strand = this_strand
                last_lines.append(line)
                last_name = feat_name

            else:
                if last_name == feat_name and last_chr == this_chr and last_strand == this_strand:
                    logger.warning("LargeIntron {} on {}".format(str(abs(last_end - int(this_start))), feat_name))
                # Prepare print
                match_chr = last_chr
                match_type = "match"
                match_start = str(min(last_pos))
                match_end = str(max(last_pos))
                match_score = '0'
                match_strand = last_strand
                match_phase = '.'
                match_feat = "ID={};Name={}".format(last_name, last_name)
                match_line = "\t".join(
                    [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                     match_feat])
                last_lines.insert(0, match_line)
                # print("\n".join(last_lines))
                output_buff += "\n".join(last_lines) + "\n"
                # Initialize
                count = 0
                last_chr = this_chr
                last_end = int(this_end)
                last_pos = [int(this_start), int(this_end)]
                last_strand = this_strand
                # Add count of same name so no duplicate name would occur
                feat_name_slim = last_name.split('-+-')[1]
                # exit()
                feat_count[feat_name_slim] += 1
                # print(feat_name_slim)
                # print(feat_count[feat_name_slim])
                try:
                    feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(feat_count[feat_parse['name']])
                except TypeError:
                    logger.warning("TypeError on feat {}".format(mylist[-1]))
                # print(feat_parse['name'])
                # print(feat_count[feat_parse['name']])
                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_lines = [line]
                last_name = feat_name
        else:
            match_chr = last_chr
            match_type = "match"
            match_start = str(min(last_pos))
            match_end = str(max(last_pos))
            match_score = '0'
            match_strand = last_strand
            match_phase = '.'
            match_feat = "ID={};Name={}".format(last_name, last_name)
            match_line = "\t".join(
                [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                 match_feat])
            last_lines.insert(0, match_line)
            # print("\n".join(last_lines))
            output_buff += "\n".join(last_lines) + "\n"
    with open(output_file, 'w') as fh:
        fh.write(output_buff)
    return output_file


def cat_est(est_files=None):
    """
    Sometimes different tissue will have duplicated fasta header name. This will rename those ests and
    cat_est workdir_isoseq_elumb/*subreads.clustered.hq.fasta.gz > total_flnc.fasta
    :param est_files:
    :return:
    """
    if ' ' in est_files:
        est_files = est_files.split()

    if 'isoseq' in est_files[0]:
        # workdir_isoseq_elumb/elumb.flcdna.pb.EuY3.leaf4.subreads.clustered.hq.fasta.gz
        match_from = r'.*.pb.(.*).subreads.*'
    elif 'Trinity' in est_files[0]:
        # elumb.lncRNA.EuG11_1.clean.fq.gz_trinity/Trinity-GG.fasta
        match_from = "(.*).clean.*"
    if est_files[0].endswith('.gz'):
        cat = 'zcat'
    else:
        cat = 'cat'
    for t in est_files:
        result = re.search(match_from, t)
        if result is not None:
            prefix = result[1]
        else:
            logging.error('Cant find prefix in {0} with pattern {1}'.format(t, match_from))
        cmd = '{} {} | sed "s/>/>{}/;s/\//_/; s/\./_/; s/\s.*//;"'.format(cat, t, prefix)
        result = sh(cmd, warning='F')
        print(result)


# 0 reference genome.fasta
# 1 protein seqeunce.fasta
# 2 est.bam
# 3 workdir
# 4 species prefix

breaker_sh = r"""
wd={3}

if [ -d $wd ]; then
    rm -r $wd
fi

REF={0}
PEP={1}
ESTBAM={2}

( time braker.pl --genome=${{REF}} --prot_seq=${{PEP}} --prg=gth --bam=${{ESTBAM}} --gth2traingenes \
--softmasking --workingdir=$wd --species {4}_braker ) &> $wd.log
"""


def braker(genome=None, protein=None, estbam=None, workdir='', send_to_augustus='F', species=''):
    """
    :param genome: genome.fasta
    :param protein: protein.fasta
    :param estbam: est.bam
    :param workdir: default is workdir_braker_genome
    :param send_to_augustus: whether to copy the trained hmm to augustus lib dir
    :return: workdir/species/Sp_5/
    """
    if workdir == '':
        workdir = 'workdir_braker_{}'.format(genome)
    if species == '':
        prefix = get_prefix(genome)
    else:
        prefix = species
    cmd = conda_act.format('braker2')
    cmd += breaker_sh.format(genome, protein, estbam, workdir, prefix)
    bsub(cmd, name='braker_train')
    return 0


def fix_comma_in_parent():
    r"""
    change Parent=CORNET00004845-1,CORNET00004845-2 to
    ...Parent=CORNET00004845-1
    ...Parent=CORNET00004845-2
    Read from stdin and print to stdout
    :return:
    """
    import fileinput
    fileinput.sys.argv = sys.argv[2:]
    for l in fileinput.input():
        found_error = re.search(r'(.*;)Parent=(.*,.*)', l)
        if found_error:
            leading = found_error[1]
            lagging = found_error[2].rstrip(';').split(',')
            for lag in lagging:
                print(leading + "Parent=" + lag + ';')
        else:
            print(l, end='')


# 0 reference_genome.fa
# 1 transcripts.fa
# 2 maker.gff
pasa_refine_sh = r"""

#singularity shell /ds3200_1/users_root/yitingshuang/lh/bin/pasapipeline_latest.sif

export PASAHOME=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/bin/PASApipeline.v2.4.1
export PATH=$PASAHOME/bin/:~/lh/anaconda3/bin/:$PATH
export PERL5LIB=$PASAHOME/PerlLib:$PASAHOME/SAMPLE_HOOKS:/ds3200_1/users_root/yitingshuang/lh/anaconda3/lib/site_perl/5.26.2/

touch {0}.sqlite

echo "## templated variables to be replaced exist as <__var_name__>

# database settings
DATABASE=$PWD/{0}.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
" >   pasa.alignAssembly.sqlite.txt

#-+- Align transcript to genome and create the original sqlite db (~3h in singularity gmap mode)
# Transcripts' name should not contain '/' character
sed 's/\///' {1} > {1}.rename
$PASAHOME/Launch_PASA_pipeline.pl \
    -c pasa.alignAssembly.sqlite.txt -C -R -g {0} \
    -t {1}.rename  \
     --ALIGNERS gmap --CPU 20

#      --ALIGNERS gmap,blat --CPU 20

#-+- Rename maker.gff(necessary?)

#-+- Sort maker.gff in gene, mRNA, exon, CDS order (other type is not needed by pasa)
awk '$3=="gene"' {2} > {2}.gff3
awk '$3=="mRNA"' {2} >> {2}.gff3
awk '$3=="CDS"' {2} >> {2}.gff3
# Fix the Parent=mRNA1,mRNA2 issue. changing them into two lines. Maker is so weird :(
awk '$3=="exon"' {2} |python -m iga.annotation.maker fix_comma_in_parent >> {2}.gff3

#-+- Load GFF < 1h 
# TODO: Also need to remove tRNA, or the program will give an error
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
    -c pasa.alignAssembly.sqlite.txt \
    -g {0} \
    -P {2}.gff3

echo "# database settings
DATABASE=$PWD/{0}.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = script_name + : + parameter
# assign a value as done above.


#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>
cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>
" > pasa.annotCompare.config

#-+- Refine genome gff, (Very slow, 6h in singularity mode)
$PASAHOME/Launch_PASA_pipeline.pl \
    -c pasa.annotCompare.config -A \
    -g {0} \
    -t {1}.rename

"""


def pasa_refine(genome=None, transcript=None, gff=None, use_grid='F'):
    r"""
    pasa_refine ref.fa flnc_rna.fasta genome.maker.gff  --use_grid T
    :param genome: the assembled genome (fasta)
    :param transcript: the assembled transcripts (fasta)
    :param gff: the existing annotation (gff3)
    :param use_grid: whether to use bsub (T or F)
    :return:
    """
    cmd = pasa_refine_sh.format(genome, transcript, gff)
    if use_grid == 'T':
        jobid = bsub(cmd, direct_submit='F', cpus=5, name='pasa_refine')
        waitjob(jobid)
    else:
        sh(cmd)


# 0 input.gff
# 1 prefix, like A188 or CORNE
maker_rename_sh = r"""
maker_map_ids --abrv_gene '' --prefix {1} --justify 8 --suffix '-t' --iterate 1 {0}  > {0}.map.txt
#Caution map_gff_ids will rewrite the file instead of generating a new one
first_only.pl  {0}.map.txt >  {0}.map_uniq.txt
# cp {0} {0}.format.gff
# map_gff_ids {0}.map_uniq.txt {0}.format.gff
"""


def maker_rename_gff(gff=None, prefix='MAKER', comment_print='F'):
    output = gff.replace('.gff3', '').replace('.gff', '') + '.format.gff'
    cmd = maker_rename_sh.format(gff, prefix)
    sh(cmd)
    gff_rename_tab_file = gff + ".map_uniq.txt"
    gff_rename_dict = {}
    with open(gff_rename_tab_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            gff_rename_dict[mylist[0]] = mylist[1]
    with open(gff) as fh, open(output, 'w') as fh_out:
        for line in fh:
            if line.startswith('#'):
                if comment_print == 'T':
                    fh_out.write(line)
                continue
            if line.strip() == '':
                continue
            result = re.search(r'ID=(.*?);', line)
            if result is not None and result[1] is not None and result[1] in gff_rename_dict:
                keyword = result[1]
                line = re.sub(keyword, gff_rename_dict[keyword], line)
            result = re.search(r'Parent=(.*?);', line.replace('\n', ';'))
            if result is not None and result[1] is not None and result[1] in gff_rename_dict:
                keyword = result[1]
                line = re.sub(keyword, gff_rename_dict[keyword], line)
            fh_out.write(line)
    return 0


def which(executable, path=None):
    """Find if 'executable' can be run. Looks for it in 'path'
    (string that lists directories separated by 'os.pathsep';
    defaults to os.environ['PATH']). Checks for all executable
    extensions. Returns full path or None if no command is found.
    """
    # https://gist.github.com/techtonik/4368898
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    for p in paths:
        f = os.path.join(p, executable)
        if os.path.isfile(f):
            return f
    else:
        return None


# environment for maker
maker_env_sh = r"""
ZOE_HMM_DIR=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/

"""

# 0 workdir
# 1 prev_round
# 2 current_round
# 3 cfg_file
# 4 cfg

maker_run_sh = r"""
cd {}
maker *ctl >> maker.out 2>> maker.err
"""

def split_fasta1(fasta = None, dir = None, size = 4000000):
    #file_list = []
    z = 1
    with open(fasta, 'r') as f:
        data = f.read().split('>')
        if os.path.exists(dir) == "True":
            #mv(dir, dir + str(time.time()).replace('.', ''))
            #mkdir(dir)
            pass
        else:
            mkdir(dir)
        z = 1
        file = f"{z}.fa"
        new_file_dir = os.path.join(dir,file)
        open(new_file_dir, "w")
        for i, j in enumerate(data[1:], start=1):
            (header, content) = j.split('\n', 1)
            #new_file_name = f'{header}.fasta'
            #new_file_dir = os.path.join(dir,new_file_name)
            #file_list.append(new_file_name)
            if os.path.getsize(new_file_dir) < size:
                with open(new_file_dir, 'a') as f:
                    f.write('>' + j)
            else:
                z += 1
                file = f"{z}.fa"
                new_file_dir = os.path.join(dir, file)
                with open(new_file_dir,"w") as f:
                    f.write(">" + j)


# split_fasta1("test/altra.genome.fa.masked.fa", dir="F:/R1")
# 分了1000多个文件

def prepare_cfg(estgff=None, pepgff=None,
              rmgff=None, round=1, species='', use_grid='T', cpus=2,
              augustus_species='', snap_hmm='', queue='Q104C512G_X4', update='', workdir = ''):
    snap_hmm_dir = '/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/'
    # default returned a string with file names, changing it into list type

    [estgff, pepgff, rmgff] = abspath_list([estgff, pepgff, rmgff])
    cfg_exe = Config('maker_exe')
    # makeblastdb=/home/yanhui/anaconda3/envs/repeat/bin/makeblastdb
    # blastn=/home/yanhui/anaconda3/envs/repeat/bin/blastn
    # blastx=/home/yanhui/anaconda3/envs/repeat/bin/blastx
    # tblastx=/home/yanhui/anaconda3/envs/repeat/bin/tblastx
    # RepeatMasker=/home/yanhui/anaconda3/envs/repeat/bin/RepeatMasker
    # exonerate=/home/yanhui/anaconda3/envs/repeat/bin/exonerate
    # snap=/home/yanhui/anaconda3/envs/repeat/bin/snap
    # augustus=/home/yanhui/anaconda3/envs/repeat/bin/augustus
    # evm=/home/yanhui/anaconda3/envs/repeat/bin/evidence_modeler.pl
    # tRNAscan-SE=/home/yanhui/anaconda3/envs/repeat/bin/tRNAscan-SE
    maker_exes = ["makeblastdb", "blastn", "blastx", "tblastx", "RepeatMasker", "exonerate", "snap", "augustus", "evm", "tRNAscan-SE"]
    abbr_to_exe = {"makeblastdb":"makeblastdb",
        "blastn":"blastn",
        "blastx":"blastx",
        "tblastx":"tblastx",
        "RepeatMasker":"RepeatMasker",
        "exonerate":"exonerate",
        "snap":"snap",
        "augustus":"augustus",
        "evm":"evidence_modeler.pl",
        "tRNAscan-SE":"tRNAscan-SE"}
    maker_exe_paths = ['']
    for exe in maker_exes:
        exe_path = which(abbr_to_exe[exe])
        if exe_path is None:
            logging.error("Can't find {} in system PATH".format(exe))
            exit(1)
        maker_exe_paths.append(exe_path)
        cfg_exe.update('{}={}'.format(exe, exe_path))

    cfg_bopts = Config('maker_bopts')
    cfg = Config('maker')
    cfg.update('est_gff={};protein_gff={};rm_gff={}'.format(estgff, pepgff, rmgff))
    cfg.write_to_file(op.join(workdir, "maker_opts.ctl"))
    cfg_exe.write_to_file(op.join(workdir, 'maker_exe.ctl'))
    cfg_bopts.write_to_file(op.join(workdir, 'maker_bopts.ctl'))
    cmd = maker_run_sh.format(workdir)


def maker_run_single_fasta(genome=None, opts=None, exe=None, bopts=None):
    """maker --genome *ctl """
    pass


def maker_run(genome=None, estgff=None, pepgff=None,
              rmgff=None, round=1, species='', use_grid='T', cpus=2,
              augustus_species='', snap_hmm='', queue='Q104C512G_X4', update=''):
    """
    Give genome and evidence, run maker gene prediction in parallel
    :param genome:
    :param estgff:
    :param pepgff:
    :param rmgff:
    :param round:
    :param species:
    :param use_grid: whether to use LSF to submit jobs
    :param augustus_species: species for augustus, if
        from round1, then it's directory name, like coriaria_contig.fa_R1
    :param snap_hmm: hmm file for snap, if from round1,
        then it's directory name, like coriaria_contig.fa_R1
    :return:
    """
    snap_hmm_dir = '/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/'
    workdir = ''
    if species == '':
        workdir = genome + '_R' + str(round)
    if os.path.exists(workdir):
        rnd = str(time.time())
        mv(workdir, workdir + rnd)
    # logger.warning(workdir)
    # exit(1)
    # Split genome and placing into working directory like:
    # coriaria_round1:
    #   chunk.1/1.fa
    #   chunk.2/2.fa
    fa_list = split_fasta(genome, workdir, 100)
    # default returned a string with file names, changing it into list type
    # logger.debug(fa_list)
    # change estgff file name in to absolute path
    # logger.warning(os.getcwd())
    [estgff, pepgff, rmgff] = abspath_list([estgff, pepgff, rmgff])
    #    logger.warning([estgff, pepgff, rmgff])
    # Preparing cfg files
    cfg_exe = Config('maker_exe')
    # makeblastdb=/home/yanhui/anaconda3/envs/repeat/bin/makeblastdb
    # blastn=/home/yanhui/anaconda3/envs/repeat/bin/blastn
    # blastx=/home/yanhui/anaconda3/envs/repeat/bin/blastx
    # tblastx=/home/yanhui/anaconda3/envs/repeat/bin/tblastx
    # RepeatMasker=/home/yanhui/anaconda3/envs/repeat/bin/RepeatMasker
    # exonerate=/home/yanhui/anaconda3/envs/repeat/bin/exonerate
    # snap=/home/yanhui/anaconda3/envs/repeat/bin/snap
    # augustus=/home/yanhui/anaconda3/envs/repeat/bin/augustus
    # evm=/home/yanhui/anaconda3/envs/repeat/bin/evidence_modeler.pl
    # tRNAscan-SE=/home/yanhui/anaconda3/envs/repeat/bin/tRNAscan-SE
    maker_exes = ["makeblastdb", "blastn", "blastx", "tblastx", "RepeatMasker", "exonerate", "snap", "augustus", "evm", "tRNAscan-SE"]
    abbr_to_exe = {"makeblastdb":"makeblastdb",
        "blastn":"blastn",
        "blastx":"blastx",
        "tblastx":"tblastx",
        "RepeatMasker":"RepeatMasker",
        "exonerate":"exonerate",
        "snap":"snap",
        "augustus":"augustus",
        "evm":"evidence_modeler.pl",
        "tRNAscan-SE":"tRNAscan-SE"}
    maker_exe_paths = ['']
    for exe in maker_exes:
        exe_path = which(abbr_to_exe[exe])
        if exe_path is None:
            logging.error("Can't find {} in system PATH".format(exe))
            exit(1)
        maker_exe_paths.append(exe_path)
        cfg_exe.update('{}={}'.format(exe, exe_path))
    cfg_bopts = Config('maker_bopts')
    cfg = Config('maker')
    cfg.update('est_gff={};protein_gff={};rm_gff={}'.format(estgff, pepgff, rmgff))
    if round == 1:
        # Only first round will be ran in direct predict mode
        cfg.update('est2genome=1;protein2genome=1')
    else:
        cfg.update('est2genome=0;protein2genome=0')
        if augustus_species != '':
            cfg.update('augustus_species={}'.format(augustus_species))
        if snap_hmm != '':
            if '.hmm' not in snap_hmm:
                # In case .hmm extension was not added in input
                snap_hmm = snap_hmm + '.hmm'
            cfg.update('snaphmm={}'.format(op.join(snap_hmm_dir, snap_hmm)))
    if update != '':
        cfg.update(update)
    #
    # get abs path of all fasta files
    os.chdir(workdir)
    fa_list = abspath_list(fa_list)
    # abspath_list(fa_list)
    # job list for storing submitted job IDs
    job_list = []
    for i in fa_list:
        fa_name = op.basename(i)
        workdir_sep = i + '.run/'
        mkdir(workdir_sep)
        mv(i, workdir_sep)
        workdir_sep = op.abspath(workdir_sep)
        # fasta = op.join(workdir, fa_name)
        cfg.update('genome={}'.format(fa_name))
        cfg.write_to_file(op.join(workdir_sep, "maker_opts.ctl"))
        cfg_exe.write_to_file(op.join(workdir_sep, 'maker_exe.ctl'))
        cfg_bopts.write_to_file(op.join(workdir_sep, 'maker_bopts.ctl'))
        cmd = maker_run_sh.format(workdir_sep)
        # sh(cmd)
        if use_grid == 'T':
            job_name = op.basename(workdir_sep)
            job_id = bsub(cmd, queue=queue, cpus=cpus, name="maker_{}".format(job_name))
            job_list.append(job_id)
            time.sleep(30)
        else:
            job_list.append(cmd)
    if use_grid == 'T':
        logger.warning("Submitted jobs:")
        logger.warning(job_list)
        waitjob(job_list)
        logger.warning("Submmited job finished, check log files to make sure they really finished")
    else:
        sh(job_list, parallel='T', cpus=cpus)


prepare_pipe_sh = r"""#!/bin/bash

set -eo

# User defined variables
PREFIX=CLPLA
REF=clpla.contig.fa
PEP=arath_med_sprot.pep
ISOSEQDIR=workdir_isoseq_clpla
# User defined variables end

python -m iga.annotation.repeat repeatmasker --species Viridiplantae --denovo F ${{REF}} &
python -m iga.annotation.repeat repeatmasker --species '' --denovo T ${{REF}} &

hisat2-build ${{REF}} ${{REF}}

python -m iga.annotation.rnaseq  reads_align_assembly ${{REF}} "cechi.ssRNAbgi.CC_Y5_1.clean.fq.gz cechi.ssRNAbgi.CC_Y5_2.clean.fq.gz" &
python -m iga.annotation.rnaseq  reads_align_assembly ${{REF}} "cechi.ssRNAbgi.CC_Y7_1.clean.fq.gz cechi.ssRNAbgi.CC_Y7_2.clean.fq.gz" &
python -m iga.annotation.rnaseq  reads_align_assembly ${{REF}} "cechi.ssRNAbgi.CC_Y9_1.clean.fq.gz cechi.ssRNAbgi.CC_Y9_2.clean.fq.gz" &

wait
python -m iga.annotation.repeat post_repeatmasker workdir_repeatmask_${{REF}} ${{REF}}

ln -s $PWD/workdir_repeatmask_${{REF}}/Full_mask/full_mask.complex.reformat.gff3 repeat.gff

REPEAT_GFF=repeat.gff

MASKEDREF=${{REF}}.masked.fa

if [ ! -e  ${{MASKEDREF}} ]
then
    echo "Error, masked genome not exists, check log"
    exit
fi

python -m iga.annotation.maker prep_genblast ${{MASKEDREF}} ${{PEP}}

ln -s ${{MASKEDREF}}.${{PEP}}..gff genblast.gff

python -m iga.annotation.maker cat_est ${{ISOSEQDIR}}/*hq.fasta.gz > flnc.fasta
python -m iga.annotation.maker cat_est *trinity/Trinity-GG.fasta > rnaseq.fasta

python -m iga.annotation.maker fastq2gff flnc.fasta ${{REF}} # output flnc.fasta.rawgff.gff
python -m iga.annotation.maker fastq2gff rnaseq.fasta ${{REF}}  # output rnaseq.fasta.rawgff.gff

cat flnc.fasta.rawgff.gff rnaseq.fasta.rawgff.gff > total_est.gff

FLNCFASTA=flnc.fasta
FLNCGFF=flnc.fasta.rawgff.gff
ESTGFF=total_est.gff
REPEAT_GFF=repeat.gff
PEPGFF=genblast.gff

python -m iga.annotation.maker maker_pipe --ref_genome ${{MASKEDREF}} \
    --flnc_GFF ${{FLNCGFF}} --est_gff ${{ESTGFF}} --flcdna_fasta ${{FLNCFASTA}} \
    --pep_gff ${{PEPGFF}} --repeat_gff ${{REPEAT_GFF}} --prefix $PREFIX > maker.sh

"""


def prepare_pipe():
    print(prepare_pipe_sh.format())


# 0 reference.fa
# 1 est.gff
# 2 flcdna.gff
# 3 pep.gff
# 4 repeat.gff
maker_pipe_sh = """#!/bin/bash
set -eo

REF={0}
FLNCESTGFF={1}
ESTGFF={2}
CDNAFASTA={3}
PEPGFF={4}
REPEATGFF={5}
PREFIX={6}


#-+-First Round
ROUND=1
#python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${{REF}} ${{FLNCESTGFF}} ${{PEPGFF}} ${{REPEATGFF}} --cpus 2
python -m iga.annotation.maker maker_check_resub ${{REF}}_R1
python -m iga.annotation.maker maker_collect     ${{REF}}_R1
python -m iga.annotation.maker maker_train       ${{REF}}_R1 --cdna_fasta ${{CDNAFASTA}}  --snap 'T' --augustus F
cd ${{REF}}_R${{ROUND}}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..


#-+-Second Round
#python -m iga.annotation.maker deploy_augustus
PREV=1
ROUND=2
python -m iga.annotation.maker maker_run         ${{REF}} ${{ESTGFF}} ${{PEPGFF}} ${{REPEATGFF}} --round ${{ROUND}} \
--augustus_species ${{REF}}_R${{PREV}}_direct --snap_hmm ${{REF}}_R${{PREV}} --update "alt_splice=1"
python -m iga.annotation.maker maker_check_resub ${{REF}}_R${{ROUND}}
python -m iga.annotation.maker maker_collect     ${{REF}}_R${{ROUND}}
python -m iga.annotation.maker maker_train       ${{REF}}_R${{ROUND}} --cdna_fasta ${{CDNAFASTA}} --augustus F
cd ${{REF}}_R${{ROUND}}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Third Round #augustus is directly trained
ROUND=3
PREV=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${{REF}} ${{ESTGFF}} ${{PEPGFF}} ${{REPEATGFF}} --round ${{ROUND}} \
--augustus_species ${{REF}}_R${{PREV}}_direct --snap_hmm ${{REF}}_R${{PREV}} --update "trna=1;alt_splice=1"
#--queue Q64C1T_X4
python -m iga.annotation.maker maker_check_resub ${{REF}}_R${{ROUND}}
python -m iga.annotation.maker maker_collect     ${{REF}}_R${{ROUND}}
cd ${{REF}}_R${{ROUND}}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Final Pasa Refine
cd ${{REF}}_R${{ROUND}}
python -m iga.annotation.maker pasa_refine ref.fa  ../$CDNAFASTA genome.maker.gff
chmod -w ref.fa.sqlite.gene_structures_post_PASA_updates.*.gff3
cp ref.fa.sqlite.gene_structures_post_PASA_updates.*.gff3 ref.fa.pasa.gff3
python -m iga.annotation.maker maker_rename_gff --prefix $PREFIX ref.fa.pasa.gff3
grep trna genome.maker.gff > trna.gff
cat ref.fa.pasa.format.gff trna.gff > ${{REF}}.gene_structure.gff3
gff_genome_to_genes.pl ${{REF}}.gene_structure.gff3 ref.fa > ${{REF}}.gene_structure.cds
cds2aa.pl ${{REF}}.gene_structure.cds > ${{REF}}.gene_structure.pep
python -m iga.assembly.assess busco --mode prot ${{REF}}.gene_structure.pep

#-+-Functional annotation
python -m iga.annotation.maker func_anno ${{REF}}.gene_structure.pep


"""


def maker_pipe(ref_genome='', flnc_GFF='', est_gff='', flcdna_fasta='', pep_gff='', repeat_gff='', prefix=''):
    """
    :param ref_genome: fasta format of reference genome
    :param est_gff: GFF format of aligned transcripts
    :param flcdna_fasta: fasta format of full lenth cDNA used to train augustus
    :param pep_gff: GFF format of aligned homologous proteins
    :param repeat_gff: GFF format of repeat elements
    :return: print maker pipeline commands
    """
    cmd = maker_pipe_sh.format(ref_genome, flnc_GFF, est_gff, flcdna_fasta, pep_gff, repeat_gff, prefix)
    print(cmd)


maker_resub_sh = r"""
cd {}
mkdir -p rm 
mv *.maker.output rm/
rm -rf rm &
maker *ctl > maker.out 2> maker.err
"""


def maker_resub(dir_list=None, queue="Q104C512G_X4", cpus=4):
    r"""
    Resubmit failed jobs by directory name
    :param dir_list:
    :param queue:
    :return:
    """
    if type(dir_list) == str:
        dir_list = [dir_list]
    logger.warning(dir_list)
    # logger.debug(queue)
    # exit(1)
    job_list = []
    for i in dir_list:
        cmd = maker_resub_sh.format(i)
        job_id = bsub(cmd, queue=queue, cpus=cpus)
        job_list.append(job_id)
        time.sleep(30)
    logger.warning("Submitted jobs:")
    logger.warning(job_list)
    waitjob(job_list)
    logger.warning("Submmited job finished, check log files to make sure they really finished")


def maker_check(workdir=None):
    """
    check whether all partitions finished
    :param workdir:
    :return:
    """
    subdir = os.listdir(workdir)
    error_list = []
    unfinished_list = []

    for sd in subdir:
        if op.isdir(sd):
            maker_log = op.join(workdir, sd, 'maker.err')
            maker_log_buff = ''
            fail_mark = 1
            try:
                with open(maker_log) as fh:
                    maker_log_buff = fh.read()
                if 'Maker is now finished!!!' in maker_log_buff:
                    if not 'ERROR' in maker_log_buff and not 'Fail' in maker_log_buff:
                        pass
                    else:
                        error_list.append(sd)
                else:
                    unfinished_list.append(sd)
            except FileNotFoundError:
                logger.error("Can't find maker error log for {}".format(sd))
                unfinished_list.append(sd)
    if len(unfinished_list) + len(error_list) == 0:
        logger.warning("Cheers! All finished without errors!")
    else:
        if len(unfinished_list) > 0:
            logger.warning("Unfinished chunks are:")
            [logger.warning(l) for l in unfinished_list]
        if len(error_list) > 0:
            logger.warning("Chunks with errors are:")
            [logger.warning(l) for l in error_list]
        # exit(1)

    return error_list + unfinished_list


def maker_check_resub(workdir=None, queue="Q104C512G_X4"):
    i = 1
    current_dir = op.abspath(os.curdir)
    # absworkdir = op.abspath(workdir)
    while i < 3:
        i += 1
        # at most resub two times
        os.chdir(current_dir)
        failed_list = maker_check(workdir)
        os.chdir(workdir)
        if len(failed_list) > 1:
            maker_resub(failed_list, queue=queue, cpus=i)
        else:
            return 0
    logger.error("Failed too many times, stop submitting")
    exit(1)
    return 1


def deploy_augustus():
    r"""
    deploy augustus config dir to /tmp/lh/ for maker use
    :return:
    """
    node_list = {}
    node_list['Q64C1T_X4'] = ['node02', 'node03', 'node04', 'node05']
    node_list['Q104C512G_X4'] = ['node10', 'node11', 'node12', 'node13']
    augustus_config_dir = '/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config'
    local_dir = '/tmp/lh_config'
    cmd = 'touch {0} && rm -fr {0} && cp -fr {1} {0}'.format(local_dir, augustus_config_dir)
    for q in node_list:
        for node in node_list[q]:
            bsub(cmd, q + ' -m {}'.format(node))
    time.sleep(120)
    return 0


# 0 working directory
collect_maker_sh = r"""
set -euxo pipefail
cd {}
#9.fa.run/9.maker.output/9_master_datastore_index.log
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       STARTED
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       FINISHED
touch total_master_datastore_index.log
rm    total_master_datastore_index.log
touch total_master_datastore_index.log
for i in `ls -d *.fa.run/`
do
    echo $i
    j=${{i%.fa.run/}}
    cat ${{i}}/${{j}}.maker.output/${{j}}_master_datastore_index.log |sed "s/\t/\t${{j}}.fa.run\/${{j}}.maker.output\//" \
    >>total_master_datastore_index.log
done


#b. Merge fasta
fasta_merge -d total_master_datastore_index.log
gff3_merge -o genome.all.gff -d total_master_datastore_index.log
gff3_merge -n -o genome.all.noseq.gff -d total_master_datastore_index.log

echo "Merge completed succefully:"

awk '$2=="maker"' genome.all.noseq.gff > genome.maker.gff
# echo "##FASTA" >> genome.maker.gff
# cat *.run/*.fa >> genome.maker.gff
cat *.run/*.fa > ref.fa

date"""


def maker_collect(workdir=None, use_grid='F'):
    """
    Collect maker result from a paralleled run in workdir
    :param workdir:
    :return:
    """
    cmd = collect_maker_sh.format(workdir)
    if use_grid == 'T':
        job = bsub(cmd, direct_submit='F', name='maker_collect')
        waitjob(job)
    else:
        res = sh(cmd)
        logger.warning(res)
    return 0


busco_export_sh = r"""
export BUSCO_CONFIG_FILE=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/busco/myconfig.ini
"""

# 0 workdir
# 1 PREFIX of this model
train_snap_sh = r"""export HMMDIR=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/
mkdir -p {0}/train_snap
cd {0}/train_snap
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi
maker2zff -x 0.25 -l 50  genome.all.gff
fathom -gene-stats genome.ann genome.dna >gene-stats.log 2>&1
fathom -validate genome.ann genome.dna >validate.log 2>&1

fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna

mkdir -p params

cd params

forge ../export.ann ../export.dna >../forge.log 2>&1

cd ..
hmm-assembler.pl snap_trained params > snap_trained.hmm

#if [ -d ${{HMMDIR}}/{1}.hmm ]
#then
#    RND=$(date +%s%N)
#    mv ${{HMMDIR}}/{1}.hmm ${{HMMDIR}}/{1}.hmm.$RND
#fi

cp snap_trained.hmm ${{HMMDIR}}/{1}.hmm

echo "Train SNAP completed succefully:"
echo "${{HMMDIR}}/{1}.hmm"
date
"""

# added optimise step as BUSCO will fail this step naturally, maybe need 5 hours to finish
# Error: IO.c: loadable library and perl binaries are mismatched (got handshake key 0xdb80080, needed 0xde00080)
# 0 workdir
# 1 prefix
train_augustus_sh = r"""
mkdir -p {0}/train_augustus
cd {0}/train_augustus
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi

if [ ! -e ref.fa ]
then
    ln -s ../ref.fa
fi

awk -v OFS="\t" '{{ if ($3 == "mRNA") print $1, $4, $5 }}' genome.all.gff | \
  awk -v OFS="\t" '{{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }}' | \
  bedtools getfasta -fi ref.fa -bed - -fo total.all.maker.transcripts1000.fasta

#Do not need it in current HPC environment
#export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

LINEAGE=embryophyta_odb10
THREADS=104
INPUT=total.all.maker.transcripts1000.fasta
OUTPUT={1}
NEWMODEL={1}
#TODO
AUGUSTUS_SPECIES=arabidopsis
AUGUSTUS_CONFIG_PATH_ORIGINAL=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config
if [ -d $OUTPUT ]
then
    rm -rf $OUTPUT
fi
busco -i $INPUT  -o $OUTPUT  -l $LINEAGE \
  -m genome -c $THREADS --long --augustus_species $AUGUSTUS_SPECIES \
  --augustus_parameters='--progress=true' >busco.out 2>busco.err

cd $OUTPUT/run_${{LINEAGE}}/augustus_output/retraining_parameters/BUSCO_${{OUTPUT}}/

rename "BUSCO_" "" *

sed -i 's/BUSCO_//g' {1}_parameters.cfg

if [ -d $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL ]
then
    RND=$(date +%s%N)
    mv $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL.$RND
fi
mkdir -p $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL
cp ./${{OUTPUT}}*  $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/{1}/

#/ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/busco/scripts/optimize_augustus.pl --cpus=8 \
#--species={1} ../../training_set.db
#
echo "Train Augustus completed succefully"
echo "Augustus species: {1}"
date

"""


def filter_gff_by_aed(gff=None, gff_out='', aed='0.2'):
    r"""
    filter maker_gff by aed value
    output as gff.filter as default
    :param gff:
    :param aed:
    :return:
    """
    buff_list = []
    buff = ''
    if gff_out == '':
        gff_out = gff + '.filter'
    with open(gff) as fh:
        for line in fh:
            buff += line
            mylist = line.split()
            if mylist[2] == 'gene':
                buff_list.append(buff)
                buff = ''
    result = ''
    for bf in buff_list:
        pattern_result = re.search(r'_AED=(.*?);', bf)[1]
        if float(pattern_result) <= float(aed):
            result += bf
    with open(gff_out) as fh:
        fh.write(result)
    return 0


# training augustus without BUSCO
# run after snap is finished
# single/8 thread, default in local run.
# 0 workdir
# 1 prefix
# 2 absolute path to full length fasta
train_augustus_direct_sh = r"""
if [ -d {0}/train_augustus_direct ]
then
    rm -rf {0}/train_augustus_direct
fi
mkdir -p {0}/train_augustus_direct
cd {0}/train_augustus_direct
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi

if [ ! -e ref.fa ]
then
    ln -s ../ref.fa
fi

NUMFOUND=500
NUMSPLIT=250
CDNA_FASTA={2}
AUGUSTUS_SPECIES_NAME={1}_direct
WORKING_DIR=$PWD
ROOT=$PWD

if [ -d $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME ]
then
    RND=$(date +%s%N)
    mv $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME.$RND
fi

ln -s ../train_snap/uni.ann
ln -s ../train_snap/uni.dna 

/ds3200_1/users_root/yitingshuang/lh/bin/GC_specific_MAKER/fathom_to_genbank.pl --annotation_file uni.ann \
 --dna_file uni.dna  --genbank_file augustus.gb \
 --number ${{NUMFOUND}}
perl -e  'while (my $line = <>){{ if ($line =~ /^LOCUS\s+(\S+)/) {{ print "$1\n"; }} }}'  ${{WORKING_DIR}}/augustus.gb \
 >  ${{WORKING_DIR}}/genbank_gene_list.txt
/ds3200_1/users_root/yitingshuang/lh/bin/GC_specific_MAKER/get_subset_of_fastas.pl  \
 -l  ${{WORKING_DIR}}/genbank_gene_list.txt   \
 -f ${{WORKING_DIR}}/uni.dna  -o  ${{WORKING_DIR}}/genbank_gene_seqs.fasta

perl ~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/randomSplit.pl ${{WORKING_DIR}}/augustus.gb ${{NUMSPLIT}}

~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME \
--genome=${{WORKING_DIR}}/genbank_gene_seqs.fasta --trainingset=${{WORKING_DIR}}/augustus.gb --cdna=$CDNA_FASTA  \
--noutr

# Failed due to mysql issue
# --pasa --useGMAPforPASA

cd ./autoAug/autoAugPred_abinitio/shells

x=1
while [ -e ./aug${{x}} ]
do
    echo "A.  $x"
    ./aug${{x}} &
    let x=x+1
done

wait

cd $WORKING_DIR

~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME \
--genome=${{WORKING_DIR}}/genbank_gene_seqs.fasta --useexisting --hints=${{WORKING_DIR}}/autoAug/hints/hints.E.gff \
 -v -v -v  --index=1

cd ${{WORKING_DIR}}/autoAug/autoAugPred_hints/shells/

let x=1
while [ -e ./aug${{x}} ]
do
    echo "B.  $x"
    ./aug${{x}} &
    let x=x+1
done

wait

echo "Successfully finished"
"""


def maker_train(workdir=None, prefix='', augustus='T', snap='T', use_grid='F', augustus_direct='T',
                cdna_fasta=''):
    """
    :param workdir:
    :return:
    """
    cmd = ''
    workdir = op.abspath(workdir)
    set_workdir = 'cd {};'.format(workdir)
    if prefix == '':
        prefix = op.basename(workdir)
    if snap == 'T':
        cmd += set_workdir + "\n" + train_snap_sh.format(workdir, prefix)
    if augustus == 'T':
        # BUSCO 4.1.2 failed to retrain augustus, even specified augustus config dir to local
        # The error was no exon_probs.pbl file produced.
        # I don't know why. For now, I used busco v4.0.1(with bug manual fixed).
        # cmd += set_workdir + "\n" + busco_export_sh + train_augustus_sh.format(workdir, prefix)
        cmd += set_workdir + "\n" + conda_act.format('busco') + busco_export_sh + \
               train_augustus_sh.format(workdir, prefix)
    if augustus_direct == 'T':
        if cdna_fasta != '':
            cdna_fasta = op.abspath(cdna_fasta)
            logger.warning(workdir)
            logger.warning(cdna_fasta)
            # TODO in future: use pasa to train augustus, may need to add pasa_export_sh before this line
            cmd += train_augustus_direct_sh.format(workdir, prefix, cdna_fasta)
        else:
            logger.error("Provide cdna.fasta before train augustus_direct")
            exit(1)
    if use_grid == 'T':
        joblist = bsub(cmd, direct_submit='F', name='maker_train', cpus=2)
        waitjob(joblist)
    else:
        sh(cmd)
    return 0


def check_repeat_overlap(maker_gff=None, repeat_gff=None):
    sh('bedtools covergae -a {} -b {}'.format(maker_gff, repeat_gff))


def check_aed_dist(maker_gff=None):
    result = sh('AED_cdf_generator.pl -b 0.025 {}'.format(maker_gff), )
    print(result)


def str_to_class(str1):
    return getattr(sys.modules[__name__], str1)


# liftover
# Require RaGOO
liftover_sh = r"""
lift_over.py 
"""


def liftover_by_agp(gff=None, agp=None):
    """
    liftover script from RaGOO
    lift over gff files based on agp file, this is used when annotation genes on contig level and need to transfer to
    chromosome level
    :param gff:
    :param agp:
    :return:
    eg:
    python -m iga.annotation.maker liftover_by_agp genome.maker.gff CORNE_v1.0.agp > genome.lifted.gff
    gff_genome_to_genes.pl genome.lifted.gff.cds CORNE_v1.0.chr.fa > genome.lifted.gff.cds.fa
    cds2aa.pl genome.lifted.gff.cds.fa > genome.lifted.gff.cds.fa.pep
    """
    # from RaGOO
    # chr name of contig
    reverse_strand = {'-': '+', '+': '-'}
    chrd = {}
    # start on chr
    startd = {}
    # end on chr
    endd = {}
    # contig strand on chr
    strandd = {}
    # contig length
    lengthd = {}
    with open(agp) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                mylist = line.split()
                contig_name = mylist[5]
                if mylist[4] == 'W':
                    chrd[contig_name] = mylist[0]
                    startd[contig_name] = int(mylist[1])
                    endd[contig_name] = int(mylist[2])
                    strandd[contig_name] = mylist[-1]
                    lengthd[contig_name] = int(mylist[-2])

    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                print(line, end='')
            else:
                if line.strip() == '':
                    continue
                mylist = line.rstrip().split()
                try:
                    this_contig = mylist[0]
                except IndexError:
                    logger.error(line)
                    exit(1)
                this_start = int(mylist[3])
                this_end = int(mylist[4])
                this_strand = mylist[6]
                # transforming
                new_chr = chrd[this_contig]
                if strandd[this_contig] == '-':
                    new_strand = reverse_strand[this_strand]
                    # new start coordinate = seqLength - endCoord
                    this_start = lengthd[this_contig] - (this_start - 1)
                    this_end = lengthd[this_contig] - (this_end - 1)
                else:
                    new_strand = this_strand
                # Add offset, because loci in AGP is 1-based, so minus 1 is the real offset
                new_start = this_start + startd[this_contig] - 1
                new_end = this_end + startd[this_contig] - 1
                # in case we need swap loci like Chr0 . gene 800 1
                if new_end < new_start:
                    (new_start, new_end) = (new_end, new_start)
                mylist[0] = new_chr
                mylist[3] = str(new_start)
                mylist[4] = str(new_end)
                mylist[6] = new_strand
                new_line = "\t".join(mylist)
                print(new_line)


# 0 pep.fasta

anno_func_sh = r"""
# blast swissprot: 
# blast tair: 
# Eggnog: 
# Interproscan: Get GO/KEGG/InterPro domains
"""


def func_anno(pep=None, interproscan='T', eggnog='F', tair='T', medtr5='T', swissprot='T'):
    """
    :param pep: fasta format of peptide to be annotated with function
    :param interproscan:
    :param eggnog: currently not stable
    :param tair:
    :param medtr5: medicagov5 pep database
    :param swissprot:
    :return:
    """
    blast_cpus = '12'
    db_database = {
        'swissprot': '/ds3200_1/users_root/yitingshuang/lh/database/function/swissprot_viridiplantae+AND+reviewed_yes_20200518.pep',
        'swissprot_desc': '/ds3200_1/users_root/yitingshuang/lh/database/function/swissprot_viridiplantae+AND+reviewed_yes_20200518.pep.description',
        'tair10': '/ds3200_1/users_root/yitingshuang/lh/database/function/tair10.pep',
        'medtr5': '/ds3200_1/users_root/yitingshuang/lh/database/function/medtrv5.pep'}
    joblist = []
    if interproscan == 'T':
        ipr_out = pep + '.ipr.tsv'
        if not os.path.exists("{0}.ipr.tsv".format(pep)):
            cmd = "interproscan.sh -f tsv -dp  -pa  -goterms -i {0} -b {0}.ipr".format(pep)
            job = bsub(cmd, name='iprscan', cpus=4)
            joblist.append(job)
    # Currently eggnog bin has a bug, so skip this step
    # if eggnog == 'T':
    #     cmd = 'emapper.py -m diamond --target_taxa "Viridiplantae" -o {0}.eggnog -i {0} --cpu 6'.format(pep)
    #     job = bsub(cmd, name='eggnog', cpus=6)
    #     joblist.append(job)
    if tair == 'T':
        tair_out = pep + ".tair10.bln"
        if not os.path.exists("{0}.tair10.bln".format(pep)):
            db = db_database['tair10']
            cmd = "blastp -subject {0} -query {1} -out {1}.tair10.bln " \
                  "-evalue 1e-5 -outfmt 6 -num_threads {2}".format(db, pep, blast_cpus)
            job = bsub(cmd, name='tair10', cpus=6)
            joblist.append(job)
    if swissprot == 'T':
        swissprot_out = pep + ".swissprot.bln"
        if not os.path.exists("{0}.swissprot.bln".format(pep)):
            db = db_database['swissprot']
            cmd = "blastp -subject {0} -query {1} -out {1}.swissprot.bln -evalue 1e-5 -outfmt 6 -num_threads " + blast_cpus
            cmd = cmd.format(db, pep)
            job = bsub(cmd, name='swissprot', cpus=6)
            joblist.append(job)
    if medtr5 == 'T':
        medtr5_out = pep + ".medtr5.bln"
        if not os.path.exists("{0}.medtr5.bln".format(pep)):
            db = db_database['medtr5']
            cmd = "blastp -subject {0} -query {1} -out {1}.medtr5.bln -evalue 1e-5 -outfmt 6 -num_threads " + blast_cpus
            cmd = cmd.format(db, pep)
            job = bsub(cmd, name='medtr5', cpus=6)
            joblist.append(job)
    waitjob(joblist)
    # func_dict = {}
    result = parse_func_result(ipr_file=ipr_out, medtr_bln=medtr5_out, tair_bln=tair_out, swissprot_bln=swissprot_out)
    print("Gene\tGO\tKEGG\tTAIR10\tMedtr5\tIPR\tPfam\tNote")
    for gene_name in result:
        template = "{}\t" * 8
        template = template.rstrip()
        g = result[gene_name]
        print(template.format(gene_name, g.go, g.kegg, g.tair, g.medtr, g.ipr, g.pfam, g.swissprot))
    return 0


class IprClass:
    def __init__(self):
        self.go = ''
        self.ipr = ''
        self.pfam = ''
        self.kegg = ''
        self.tair = ''
        self.medtr = ''
        self.swissprot = ''

    def add_ipr(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.ipr:
                if self.ipr != '':
                    self.ipr += ","
                self.ipr += t

    def add_go(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.go:
                if self.go != '':
                    self.go += ","
                self.go += t

    def add_kegg(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.kegg and t.startswith('KEGG'):
                if self.kegg != '':
                    self.kegg += ","
                self.kegg += t

    def add_pfam(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.pfam:
                if self.pfam != '':
                    self.pfam += ","
                self.pfam += t

    def add_tair(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.tair:
                if self.tair != '':
                    self.tair += ","
                self.tair += t

    def add_medtr(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.medtr:
                if self.medtr != '':
                    self.medtr += ","
                self.medtr += t

    def add_swissprot(self, token):
        tokenl = token.split('|')
        for t in tokenl:
            if t != '-' and t not in self.swissprot:
                if self.swissprot != '':
                    self.swissprot += ","
                self.swissprot += t


def best_hit_from_blast(bln=None, print_out='F'):
    """
    Get best (First hits from blast result (format 6)
    :param bln:
    :return:
    """
    hit_ortho = {}

    if not os.path.exists(bln):
        logging.error("{} do not exist".format(bln))
        return hit_ortho

    with open(bln) as fh:
        for line in fh:
            mylist = line.strip().split()
            if mylist[0] not in hit_ortho:
                hit_ortho[mylist[0]] = mylist[1]

    if print_out == 'T':
        for k in hit_ortho:
            print("{}\t{}".format(k, hit_ortho[k]))
    return hit_ortho


### parsing the IPRScan output
### making gene-model wise annotation list
def parse_func_result(ipr_file=None, medtr_bln='', tair_bln='', swissprot_bln=''):
    # 1-based column
    # 4,5
    # Pfam    PF00657
    #
    # 12
    # IPR001087
    #
    # 14
    # GO
    #
    # 15
    # KEGG
    ipr_out = defaultdict(IprClass)
    with open(ipr_file) as fh:
        for line in fh:
            line = line.strip()
            token = line.split('\t')
            gene_id = token[0]
            token_ipr = token[11]
            token_pfam = token[4]
            if not token_pfam.startswith('PF'):
                token_pfam = '-'
            # GO:12312|GO:12321 -
            try:
                token_go = token[13]
            except IndexError:
                token_go = '-'
            # KEGG: 21|KEGG: 32 -
            try:
                token_kegg = token[14]
            except IndexError:
                token_kegg = '-'
            if token_ipr != '-':
                ipr_out[gene_id].add_ipr(token_ipr)
                ipr_out[gene_id].add_go(token_go)
                ipr_out[gene_id].add_kegg(token_kegg)
                ipr_out[gene_id].add_pfam(token_pfam)
    if tair_bln != "":
        tair_ortho = best_hit_from_blast(tair_bln)
        for g in tair_ortho:
            ipr_out[g].add_tair(tair_ortho[g])
    if medtr_bln != "":
        medtr_ortho = best_hit_from_blast(medtr_bln)
        for g in medtr_ortho:
            ipr_out[g].add_medtr(medtr_ortho[g])
    if swissprot_bln != "":
        swiss_ortho = best_hit_from_blast(swissprot_bln)
        for g in swiss_ortho:
            ipr_out[g].add_swissprot(swiss_ortho[g])
    return ipr_out


def add_func(gff=None, table=None, tag='GO', pos='2'):
    r"""
    Embed function information into gff files:
    :param gff:
    :param table: geneID\tGO:123\tKO123\tAcetalysase
    :param tag: The type of functional annotation, like GO or KO or note, support multiple value like "GO,KO,Note"
    :param pos: The # column (starting from 1) of tag items, support multiple pos like "2,3,4"
    :return: the content also print to screen
    """
    # Gene is default in first column
    gene_pos = 0
    # Convert str to list
    if ',' in tag:
        tag_list = tag.split(',')
        pos_list = pos.split(',')
        if len(tag_list) != len(pos_list):
            logger.error("Error: Number of tags {} do not equal number of positions {}".format(
                len(tag_list), len(pos_list)))
    else:
        tag_list = [tag]
        pos_list = [pos]
    # Input from command line in 1-based
    # But python is 0-based. so this block is to fix 0 and 1 issue
    for iter in range(0, len(pos_list)):
        pos_list[iter] = int(pos_list[iter]) - 1
    # Pre format complete
    gff_db = GFF(gff)
    logger.warning("Reading GFF complete")
    with open(table) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            mylist = line.rstrip('\n').split('\t')
            gene_id = mylist[gene_pos]
            if gene_id not in gff_db.GFF_dict:
                continue
            for iter in range(0, len(pos_list)):
                this_tag = tag_list[iter]
                this_pos = pos_list[iter]
                try:
                    real_val = mylist[this_pos]
                    if ';' in real_val:
                        real_val = real_val.replace(';', '')
                    gff_db.GFF_dict[gene_id].update_tag(this_tag, real_val)
                except IndexError:
                    logger.error("Error on line {}, list {} and pos {}".format(line, mylist, this_pos))
                    continue
    logger.warning("Modifying GFF complete")
    gff_db.print_out()


# function

# def main():
#     """
#     the main function
#     """
#     emain()
#     # actions = (
#     #     ('isoseq', 'extract isoseq flnc reads from subreads.bam')
#     #     ('fastq2gff', 'map fastq to reference genome and get gff files'),
#     # )
#
#     # p = ActionDispatcher(actions)
#     # p.dispatch(globals())
#
#     # print(__file__)
#     # print(__doc__)
#     # exit()
#     # prog_name = "busco_wrapper"
#     # usage = "run busco on selected GENOME"
#     #
#     # parser = argparse.ArgumentParser(
#     #     prog=prog_name,
#     #     formatter_class=argparse.RawDescriptionHelpFormatter,
#     #     description=textwrap.dedent(usage),
#     #     epilog="")
#     # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
#     # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
#     # args = parser.parse_args()
#     #
#     # busco(args.GENOME)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()
