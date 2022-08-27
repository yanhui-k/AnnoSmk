"""
code used in A188 project
"""
import copy
import os
import re
from collections import defaultdict
from statistics import mean

import pandas as pd
from parse import parse

from iga.annotation.gff import Loci, Bed, GFF
from iga.apps.base import emain, qsub, get_prefix, sh, split_fasta

import logging
import coloredlogs
import os.path as op
from Bio import SeqIO

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


# # 0 ref fasta
# # 1 qry fasta
# nucmer_sh = r"""
# # Whole genome alignment. Any other alignment can also be used.
# WORKDIR={0}.{1}.nucmer
# mkdir -p $WORKDIR
# cd $WORKDIR
# ln -s ../{0}
# ln -s ../{1}
# export PATH=/lustre/home/liuhui/bin/mummer4/bin:$PATH
# # nucmer --maxmatch -c 100 -b 500 -l 50 {0} {1}
# nucmer --batch 1 -c 100 -b 500 -l 50 {0} {1}
# # Remove small and lower quality alignments
# delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
# # Convert alignment information to a .TSV format as required by SyRI
# show-coords -THrd out.filtered.delta > out.filtered.coords
# """


# syri_sh = r"""
# SYRI=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/syri
# PLOTSR=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/plotsr
# python3 $SYRI -c out.filtered.coords -d out.filtered.delta -r {0} -q {1}
# python3 $PLOTSR syri.out {0} {1} -H 8 -W 5
# """


# def nucmer(ref=None, qry=None, threads=3):
#     cmd = nucmer_sh.format(ref, qry)
#     prefix = get_prefix(ref)
#     prefix += get_prefix(qry)
#     if len(ref.split('.')) > 2:
#         chr_id = ref.split('.')[-1]
#         prefix += chr_id
#     qsub(cmd, cpus=threads, name='nucmer.'+prefix, sub=False)


# def merge_nucmer_result():
#     pass


# def syri(ref=None, qry=None, threads=4):
#     cmd = nucmer_sh.format(ref, qry) + '\nconda activate syri\n' + syri_sh.format(ref, qry)
#     prefix = get_prefix(ref)
#     prefix += get_prefix(qry)
#     if len(ref.split('.')) > 2:
#         chr_id = ref.split('.')[-1]
#         prefix += chr_id
#     qsub(cmd, cpus=threads, name='syri.'+prefix)


class LociPE:
    """
    Two loci objects
    """

    def __init__(self, left_chr, left_start, left_end, left_strand,
                 right_chr, right_start, right_end, right_strand, name):
        self.left = Loci(left_chr, left_start, left_end, name, '.', left_strand)
        self.right = Loci(right_chr, right_start, right_end, name, '.', right_strand)

    def get_line(self):
        """
        return string of this object
        :return:
        """
        result = "\t".join([self.left.chr, str(self.left.start), str(self.left.end),
                            self.right.chr, str(self.right.start), str(self.right.end),
                            self.right.name, '.', self.left.strand, self.right.strand]) + "\n"
        return result


class BedPE:
    """
    parsing lastz or Syri result into BedPE object
    """

    def __init__(self, input_file='', type='bedpe'):
        # Nested storage system
        # level1: [dict] chromosome
        # level2: [list] Loci
        self.bedpe_db = defaultdict(list)
        self.bedpe_db_right_chr = defaultdict(list)
        self.type = type
        if input_file != '':
            self.read(input_file)

    def read(self, input_file):
        """
        Read a syri file and store
         A188_2  1307    3091    -       -       B73_2   12548   14331   SYNAL1  SYN1    SYNAL   -
        :return:
        """
        with open(input_file) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                mylist = line.rstrip().split()
                if self.type == 'syri':
                    (left_chr, left_start, left_end, undef, undef,
                     right_chr, right_start, right_end, name, name1, align_type, undef) = mylist
                    left_strand = '+'
                    right_strand = '+'
                    # Default input for syri is 1-based for both start and end, so have a change to bedpe type
                    right_start = int(right_start) - 1
                    left_start = int(left_start) - 1
                elif self.type == 'bedpe':
                    # [0-based start, 0-based end)
                    (left_chr, left_start, left_end, right_chr, right_start, right_end, name,
                     undef, left_strand, right_strand) = mylist[0:10]
                elif self.type == 'lastz':
                    # lastz is also compatible with bed format: [0-based start, 0-based end)
                    # name1  zstart1 end1 name2   strand2 zstart2+  end2+ identity idPct coverage covPct  cigarx-
                    (left_chr, left_start, left_end, right_chr, right_strand, right_start, right_end,
                     identity, idcPct, coverage, covPct, cigarx) = mylist
                    left_strand = '.'
                    name = '.'

                this_lp = LociPE(left_chr, left_start, left_end, left_strand,
                                 right_chr, right_start, right_end, right_strand, name)
                if self.type == 'lastz':
                    # if len(self.bedpe_db[left_chr]) >= 2 and \
                    #         self.bedpe_db[left_chr][-1].left.start > this_lp.left.start:
                    #     break
                    # TODO just a tempory fix for nonincrement alignment
                    if len(self.bedpe_db[left_chr]) >= 2 and \
                            (self.bedpe_db[left_chr][-1].left.end > this_lp.left.start or
                             self.bedpe_db[left_chr][-1].right.end > this_lp.right.start) and \
                            (self.bedpe_db[left_chr][-2].left.end < this_lp.left.start and
                             self.bedpe_db[left_chr][-2].right.end < this_lp.right.start):
                        self.bedpe_db[left_chr].pop()
                    elif len(self.bedpe_db[left_chr]) < 1 or \
                            self.bedpe_db[left_chr][-1].left.end < this_lp.left.start and \
                            self.bedpe_db[left_chr][-1].right.end < this_lp.right.end:
                        self.bedpe_db[left_chr].append(this_lp)
                else:
                    self.bedpe_db[left_chr].append(this_lp)
                    self.bedpe_db_right_chr[right_chr].append(this_lp)
                    # logging.debug(left_chr)
                    # logging.debug(right_chr)
                    # exit()

    def stat(self, short):
        """
        stat bedpe size
        :param short:, whether to use short format, T is use
        :return:
        """
        size_list_left = []
        size_list_right = []
        size_left_ins = []
        size_right_ins = []
        size_unknown_left = []
        size_unknown_right = []
        threshold = 50
        for chr_id in self.bedpe_db:
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                left_size = lp.left.get_size(bed_format=True)
                right_size = lp.right.get_size(bed_format=True)
                if left_size == 1:
                    # Skip loci with size of 1
                    left_size = 0
                if right_size == 1:
                    right_size = 0
                size_list_left.append(left_size)
                size_list_right.append(right_size)
                if left_size < threshold < right_size:
                    size_right_ins.append(right_size)
                elif left_size > threshold > right_size:
                    size_left_ins.append(left_size)
                elif left_size > threshold and right_size > threshold:
                    size_unknown_left.append(left_size)
                    size_unknown_right.append(right_size)

        pd.set_option('display.float_format', lambda x: '%.0f' % x)

        header = ['Left', 'Right', "LeftIns", "RightIns", "Left Mosaic", "Right Mosaic"]
        tables = [size_list_left, size_list_right, size_left_ins, size_right_ins, size_unknown_left, size_unknown_right]

        # df = pd.DataFrame(np.array(tables), columns=header)
        # print(df.describe())
        if short == 'T':
            print("{}\t{}".format(sum(tables[0]), sum(tables[1])))
        else:
            print("\t", end='')
            for h in header:
                print("{:<15}".format(h), end='')
            print('')

            print("max", end="\t")
            for i in range(0, len(header)):
                try:
                    mt = max(tables[i])
                except ValueError:
                    mt = 0
                print("{:<15}".format(mt), end='\t')
            print('')

            print("mean", end="\t")
            for i in range(0, len(header)):
                try:
                    mm = int(mean(tables[i]))
                except ValueError:
                    mm = 0
                print("{:<15}".format(mm), end='\t')
            print('')

            print("sum", end="\t")
            for i in range(0, len(header)):
                try:
                    sm = sum(tables[i])
                except ValueError:
                    sm = 0
                print("{:<15}".format(sm), end='\t')
            print('')
        return 0

        # for i, v in enumerate(header):
        #     print(header[i])
        #     if len(tables[i]) < 1:
        #         continue
        #     df = pd.DataFrame(tables[i])
        #     print(df.describe())

    def get_mosaic(self, outtable=''):
        """
        Get mosaic regions of this bedpe file
        CHANGE0128: Use increment alignment only
        :return:
        """
        # logger.debug("chrid {}".format(self.bedpe_db.keys()))
        complement_db = BedPE()
        count = 0
        for chr_id in self.bedpe_db:
            # logger.debug("chrid {}".format(chr_id))
            chr_lp = self.bedpe_db[chr_id]
            for i, lp in enumerate(chr_lp):
                if i > 0:
                    left_start = chr_lp[i - 1].left.end + 1
                    left_end = chr_lp[i].left.start - 1
                    right_start = chr_lp[i - 1].right.end + 1
                    right_end = chr_lp[i].right.start - 1
                    if left_end <= left_start:
                        left_end = left_start
                    if right_end <= right_start:
                        right_end = right_start
                    # Skipping overlaping bedpes, rare but exists
                    if chr_lp[i - 1].left.end >= chr_lp[i].left.start and \
                            chr_lp[i - 1].right.end >= chr_lp[i].left.start:
                        continue
                    # From now, 1-based is transformed into 0-based start and 1-based end.
                    # Not really, 0610.2021
                    name = chr_lp[i - 1].right.name
                    if name == '.':
                        name = "SYN" + str(count)
                    chr_numeric_id = re.sub('.*_', '', chr_lp[i - 1].left.chr)
                    name = chr_numeric_id + "_" + "NOT" + name
                    new_lp = LociPE(chr_lp[i - 1].left.chr, left_start, left_end,
                                    chr_lp[i - 1].left.strand,
                                    chr_lp[i - 1].right.chr, right_start, right_end,
                                    chr_lp[i - 1].right.strand, name)
                    complement_db.bedpe_db[chr_id].append(new_lp)
                    count += 1
        complement_db.write_to_table(outtable)

    def write_to_table(self, table='', left_only=False, right_only=False):
        """
        write bedpe object into a table
        :param table: output table name
        :parameter left_only: print left bed only
        :parameter right_only: print right bed only
        :return:
        """
        result = ''
        for k in sorted(self.bedpe_db.keys()):
            # logger.debug(k)
            for i in self.bedpe_db[k]:
                if left_only:
                    result += i.left.get_line()
                elif right_only:
                    result += i.right.get_line()
                else:
                    result += i.get_line()
        if table != '':
            with open(table, 'w') as fh:
                fh.write(result)
        else:
            print(result, end='')

    def exists(self, bedpe_loci, wobble=100, type='ll', output='l'):
        """
        Test whether bedpe_loci also exists in in self.bedpe_db, with flexible boundary defined as wobble,
        and comparison type defined with type (left to left or right to right)
        :param bedpe_loci: the bedpe class object
        :param wobble: int, the boundaries can be flexible with plus or minus wobble
        :param type: the alignment type, could be ll, lr, rl, rr.
        :return:
        """
        abbrev = {"l": "left", "r": "right"}
        result = ''
        # bedpe_loci = LociPE()
        # bed_loci = Loci()
        if type[0] == 'l':
            # qry = bedpe_loci.abbrev[type[0]]
            qry = bedpe_loci.left
        else:
            qry = bedpe_loci.right

        if type[1] == 'l':
            pe_dict = self.bedpe_db
        else:
            pe_dict = self.bedpe_db_right_chr
        # logging.debug(qry.get_line())
        # logging.debug(qry.chr)
        # logging.debug(pe_dict.keys())
        # exit()
        for ref_pe in pe_dict[qry.chr]:
            # logging.debug(ref_pe.get_line())
            # exit()
            if type[1] == 'l':
                ref = ref_pe.left
            else:
                ref = ref_pe.right
            if abs(ref.start - qry.start) <= wobble and \
                    abs(ref.end - qry.end) <= wobble:
                if 'l' in output:
                    result += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    result += ref_pe.get_line() + "\t"
                result = result.rstrip() + '\n'
        return result

    def hotspot(self, bedpe_loci, wobble=100, type='ll', output='l'):
        """
        Test whether bedpe_loci also exists in in self.bedpe_db, with flexible boundary defined as wobble,
        and comparison type defined with type (left to left or right to right)
        :param bedpe_loci: the bedpe class object
        :param wobble: int, the boundaries can be flexible with plus or minus wobble
        :param type: the alignment type, could be ll, lr, rl, rr.
        :return:
        """
        abbrev = {"l": "left", "r": "right"}
        start = ''
        end = ''
        # bedpe_loci = LociPE()
        # bed_loci = Loci()
        if type[0] == 'l':
            # qry = bedpe_loci.abbrev[type[0]]
            qry = bedpe_loci.left
        else:
            qry = bedpe_loci.right

        if type[1] == 'l':
            pe_dict = self.bedpe_db
        else:
            pe_dict = self.bedpe_db_right_chr
        # logging.debug(qry.get_line())
        # logging.debug(qry.chr)
        # logging.debug(pe_dict.keys())
        # exit()
        for ref_pe in pe_dict[qry.chr]:
            # logging.debug(ref_pe.get_line())
            # exit()
            if type[1] == 'l':
                ref = ref_pe.left
            else:
                ref = ref_pe.right
            if abs(ref.start - qry.start) <= wobble:
                # logging.debug('start found')
                if 'l' in output:
                    start += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    start += ref_pe.get_line() + "\t"
                start = start.rstrip() + "\n"
                # logging.debug(start)
            if abs(ref.end - qry.end) <= wobble:
                # logging.debug('end found')
                if 'l' in output:
                    end += bedpe_loci.get_line().rstrip() + "\t"
                if 'r' in output:
                    end += ref_pe.get_line() + "\t"
                end = end.rstrip() + "\n"
                # logging.debug(end)
        # logging.debug([start, end])
        return [start, end]


def bedpe_stat(bedpe_file=None, short='F'):
    """
    stat bedpe file
    :param bedpe_file:
    :return:
    """
    bedpe = BedPE(bedpe_file)
    bedpe.stat(short)


def bed_stat(bed_file=None, short='F'):
    """
    stat bed file
    :param bed_file:
    :return:
    """
    bed = Bed(bed_file)
    bed.stat(short)


# SyRI related utils
syri_sh = r"""
REF={0}
QRY={1}
WORKDIR={0}.{1}.syri
mkdir -p $WORKDIR
cd $WORKDIR
ln -s ../{0}
ln -s ../{1}
minimap2 -ax asm5 --eqx $REF $QRY | samtools view -bS  > $REF.$QRY.bam
SYRI=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/syri
PLOTSR=/lustre/home/liuhui/project/buzzo/syri/bin/syri-1.3/syri/bin/plotsr

python3 $SYRI -c $REF.$QRY.bam -r $REF -q $QRY -k -F B --lf $REF.$QRY.log --prefix $REF.$QRY.
# python3 $SYRI -c out.filtered.coords -d out.filtered.delta -r {0} -q {1}

for i in .snps.txt .notAligned.txt .ctxOut.txt .sv.txt .synOut.txt .dupOut.txt .invDupOut.txt .invTLOut.txt .TLOut.txt .invOut.txt
do 
    rm $REF.$QRY$i
done
python3 $PLOTSR $REF.$QRY.syri.out {0} {1} -H 8 -W 5
"""


def syri(ref=None, qry=None, threads=6, submit='T', bychr='F', node='rock0[12]'):
    """
    Syri Wrapper
    :param ref:
    :param qry:
    :param threads:
    :param submit:
    :param bychr: [F|T], whether to split into seperate chromosomes and submit.
    :return:
    """
    cmd = 'conda activate syri' + syri_sh.format(ref, qry)
    prefix = get_prefix(ref)
    prefix += get_prefix(qry)
    if len(ref.split('.')) > 2:
        chr_id = ref.split('.')[-1]
        prefix += chr_id
    if bychr == 'F':
        if submit == 'T':
            qsub(cmd, cpus=threads, name='syri.' + prefix)
        else:
            sh(cmd)
    else:
        ref_list = split_fasta(ref, '', bypart='T')
        qry_list = split_fasta(qry, '', bypart='T')
        for chr_ref, chr_qry in zip(ref_list, qry_list):
            prefix = get_prefix(chr_ref)
            cmd = 'conda activate syri' + syri_sh.format(chr_ref, chr_qry)
            if node != '':
                qsub(cmd, cpus=threads, name='syri.' + prefix, node=node)
            else:
                qsub(cmd, cpus=threads, name='syri.' + prefix, node=node)


def syri_batch(genome_list=None, bychr='F', node='rock0[12]', threads=8):
    """
    batch syri wrapper
    :param genome_list: multiple genome fasta, like a.genome b.genome c.genome
    :return:
    """
    if type(genome_list) is list:
        pass
    elif ' ' in genome_list:
        genome_list = genome_list.split()
    else:
        logging.error("Genome input error: {0}".format(str(genome_list)))
        exit(1)
    for i in range(0, len(genome_list)):
        for j in range(i+1, len(genome_list)):
            qry = genome_list[i]
            ref = genome_list[j]
            if qry == ref:
                continue
            else:
                syri(qry, ref, bychr=bychr, node=node, threads=threads)


def synal_to_mosaic(synal_file=None, syriout='F', syn_tag='SYNAL', output=''):
    """
    synal_to_mosaic syri.synal.txt
    or synal_to_mosaic syri.out
    convert syri.out result to bedpe like result by complementing syntenic regions
    Input eg:
    A188_2  1307    3091    -       -       B73_2   12548   14331   SYNAL1  SYN1    SYNAL   -
    A188_2  12474   16885   -       -       B73_2   31212   35640   SYNAL2  SYN2    SYNAL   -
    A188_2  16893   17371   -       -       B73_2   37259   37738   SYNAL3  SYN2    SYNAL   -
    A188_2  18572   20443   -       -       B73_2   39049   41277   SYNAL4  SYN2    SYNAL   -
    A188_2  25932   26692   -       -       B73_2   41807   42567   SYNAL5  SYN2    SYNAL   -
    Output eg:
    A188_2  3092    12473    -       -       B73_2   14332   31211   NONSYNAL1  SYN1    SYNAL   -
    A188_2  16886   16892   -       -       B73_2   35641   37258    NONSYNAL2  SYN2    SYNAL   -
    :param synal_file: Input alignment file
    :param syriout: [F/T] whether this is syri.out file
    :return:
    """
    # logging.debug('abc')
    if re.search(r'syri.out$', synal_file):
        syriout = 'T'
    if syriout == 'T':
        sh("""awk '$11=="{1}"' {0} > {0}.{1}""".format(synal_file, syn_tag))
        synal_file += '.' + syn_tag
    bedpe = BedPE(synal_file, type='syri')
    if output != '':
        bedpe.get_mosaic(outtable=output)
    else:
        output = synal_file + '.mosaic'
        bedpe.get_mosaic(output)


def sv_within_contig():
    r"""
    Find SV within contigs.
    
    :return:
    """
    pass


# /lustre/home/liuhui/bin/bundle/fagap2bed.py fasta > gap.bed
# bedtools intersect -b A188.gap.bed -a A188_B73.total.SYN.mosaic >  A188_B73.total.SYN.mosaic.with_gap
# bedtools intersect -b A188.gap.bed -a A188_B73.total.SYN >  A188_B73.total.SYN.with_gap
# unselectItem.pl 8 8 A188_B73.total.SYN.with_gap A188_B73.total.SYN > A188_B73.total.SYN.no_gap
# unselectItem.pl 6 6 A188_B73.total.SYN.mosaic.with_gap A188_B73.total.SYN.mosaic > A188_B73.total.SYN.mosaic.no_gap
# grep NOT combined.bedpe.filter  |awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'


def format_syri_offset(offset1=None, offset2=None, syri_file=None, pos1='2,3', pos2='7,8', print_out='T'):
    r"""
    %s 100 10 > syri.out. coordinates will be added (left +99, right +9) automatically
    offset is where the segment is cut. like fa_pos.pl 100 110 > qry.fa. then it's 100 for offset1
    Why to use:
        Sometimes aligners will fail for any reason, in this way we cut the failed fasta and align
        them seperatly. Usually this will fix the previous problem. But the new coordinate is on
        the segment, so we use this utility to map coordinates to original coordinate systems.
    :param offset1: The original start of segment on original left fasta. s 100
    :param offset2: The original start of segment on original left fasta
    :param syri_file: (1-based alignment file), works for syri.out now, should support more files in future
    :param pos1: Comma delimited column number. for which to change with offset1, like 2,3 are column 2 and 3 (1-based)
    :param pos2: the same as above, except those columns will be changed with offset2
    :param print: whether to print to screen [T] or return [F]
    :return: str
    """
    pos1_list = pos1.split(',')
    pos2_list = pos2.split(',')
    # 100-1000 → 1:901, offset = 100, 1+100-1 = 100, 901 + 100 - 1 = 1000
    offset1 = int(offset1) - 1
    offset2 = int(offset2) - 1
    for i in range(0, len(pos1_list)):
        pos1_list[i] = int(pos1_list[i]) - 1
    for i in range(0, len(pos2_list)):
        pos2_list[i] = int(pos2_list[i]) - 1
    # id_increment = 10000000
    buffer = ""
    with open(syri_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            for c1 in pos1_list:
                # logging.debug(c1)
                try:
                    mylist[c1] = str(int(mylist[c1]) + offset1)
                except ValueError:
                    pass
            for c2 in pos2_list:
                # logging.debug(c2)
                try:
                    mylist[c2] = str(int(mylist[c2]) + offset2)
                except ValueError:
                    pass
            # if 'SYN' in mylist[8]:
            for col in range(8, 10):
                # Rename SYN1 to SYN10001 to avoid same name conflict with previous SYN tags
                search = re.search(r'(\D+)(\d+)', mylist[col])
                if search:
                    tag = search[1]
                    id = search[2]
                    # id = str(id_increment + int(id))
                    id = 'p' + id
                    mylist[col] = tag + id
            line = "\t".join(mylist)
            if print_out == 'T':
                print(line)
            else:
                buffer += line + "\n"
    return buffer


def synal_to_paf(synal_file=None):
    """
    :param synal_file:
    :return:
    """
    # Description
    # 1 	string 	Query sequence name
    # 2 	int 	Query sequence length
    # 3 	int 	Query start (0-based; BED-like; closed)
    # 4 	int 	Query end (0-based; BED-like; open)
    # 5 	char 	Relative strand: "+" or "-"
    # 6 	string 	Target sequence name
    # 7 	int 	Target sequence length
    # 8 	int 	Target start on original strand (0-based)
    # 9 	int 	Target end on original strand (0-based)
    # 10 	int 	Number of residue matches
    # 11 	int 	Alignment block length
    # 12 	int 	Mapping quality (0-255; 255 for missing)
    bedpe = BedPE(synal_file, type='syri')
    for chr in bedpe.bedpe_db:
        for lp in bedpe.bedpe_db[chr]:
            paf_list = [lp.left.chr, 1000000, lp.left.start, lp.left.end, '+',
                        lp.right.chr, 1000000, lp.right.start, lp.right.end,
                        min(lp.left.end - lp.left.start, lp.right.end - lp.right.start),
                        max(lp.left.end - lp.left.start, lp.right.end - lp.right.start),
                        0]
            for i in range(0, len(paf_list)):
                paf_list[i] = str(paf_list[i])
            print("\t".join(paf_list))
    return 0


def fix_syri_end(syri_out=None, qry_fa=None, ref_fa=None, postprocess='F', threads=1):
    """
    SyRI fails to find syntenic region at the end of chromosome1. resulting 40 Mb FP unsyntenic regions
    This script is used to fix this problem. The input argument is directory where syri was executed
    postprocess: [T]/F, First round run with F, Second round run with T. T is used to generated curated result
    :return:
    """
    # A188.genome.chr1.B73.genome.chr1.syri
    with open(syri_out) as f:
        for line in f:
            mylist = line.split()
            if mylist[10] == "SYN":
                last_line = line
    mylist = last_line.rstrip().split()
    qry_offset = int(mylist[2]) + 1
    qry_chr = mylist[0]
    ref_offset = int(mylist[7]) + 1
    ref_chr = mylist[5]

    qry_tail_fa = "{}_tail.fa".format(qry_fa)
    ref_tail_fa = "{}_tail.fa".format(ref_fa)

    # Do not wipe previous runs
    putative_dir = "{}.{}.syri".format(qry_tail_fa, ref_tail_fa)
    if not op.exists(putative_dir) or postprocess != "T":
        qry_fadt = SeqIO.to_dict(SeqIO.parse(qry_fa, "fasta"))
        ref_fadt = SeqIO.to_dict(SeqIO.parse(ref_fa, "fasta"))
        with open(qry_tail_fa, 'w') as out_qry:
            out_qry.write(qry_fadt[qry_chr][qry_offset - 1:].format('fasta'))
        with open(ref_tail_fa, 'w') as out_ref:
            out_ref.write(ref_fadt[ref_chr][ref_offset - 1:].format('fasta'))
        syri(qry_tail_fa, ref_tail_fa, submit='T', threads=threads, node='')
    if postprocess == "T":
        formated_SYN = format_syri_offset(qry_offset, ref_offset,
                                          "{}/{}.{}.syri.out".format(putative_dir, qry_tail_fa, ref_tail_fa),
                                          print_out='F')
        curated_SYN_file = syri_out + '.curated'
        with open(curated_SYN_file, 'w') as fo, \
                open(syri_out, 'r') as fi:
            original_syri = fi.read()
            fo.write(original_syri)
            fo.write(formated_SYN)
    # synal_to_mosaic(curated_SYN_file, syri_out='T', output=curated_SYN_file + ".mosaic")


###
def split_paf(paf_file=None, bed_file=None, bin_size=1000000, offset='T'):
    """
    split minimap paf file into 1M-window seperated, for convenient plot with gggenome
    :param paf_file: arranged in A188 B73, B73 Mo17 order
    :param bed_file: the normal gene bed file. also transforming bed_file into split gff
    :return:
    """
    # Store split result in a list
    window_list = ['']
    # Store delimeter in a dict, like dict['a']=[1000000, 2000000]
    boundary_dict = {}
    bin_size = int(bin_size)
    longest_chr = 100000000000
    max_end = bin_size * int(longest_chr / bin_size)
    last_left_chr = ''
    last_right_chr = ''
    with open(paf_file) as fh:
        line = fh.readline()
        boundary_dict[line.split()[0]] = []
        boundary_dict[line.split()[5]] = []
        last_left_chr = line.split()[0]
        last_right_chr = line.split()[5]
        for i in range(bin_size, max_end, bin_size):
            boundary_dict[line.split()[0]].append(i)
    # Window number
    window_id = 0
    # Used to switch the two side
    side_list = ['left', 'right']
    known_side = 'left'
    unknown_side = 'right'
    last_unknown_end = 0
    with open(paf_file) as fh:
        for line in fh:
            this_line = defaultdict(dict)
            line_list = line.split()
            this_line['left']['chr'] = line_list[0]
            this_line['left']['start'] = line_list[2]
            this_line['left']['end'] = line_list[3]
            this_line['right']['chr'] = line_list[5]
            this_line['right']['start'] = line_list[7]
            this_line['right']['end'] = line_list[8]
            # (this_line['left']['chr'], undef1, this_line['left']['start'], this_line['left']['end'], undef2,
            #  this_line['right']['chr'], undef3, this_line['right']['start'], this_line['right']['end'], undef,
            #  undef, undef) = line.split()
            chr_id = this_line[known_side]['chr']
            # logging.debug(line)
            # logging.debug(chr_id)
            # logging.debug(this_line[known_side]['end'])
            # logging.debug(boundary_dict[chr_id][window_id])
            flag = False
            if last_left_chr == line_list[0] and last_right_chr == line_list[5]:
                flag = True
            if flag:
                # logging.debug(line)
                # logging.debug(chr_id)
                # logging.debug(window_id)
                # TODO: in current version start - offset can be negative
                if int(this_line[known_side]['end']) <= boundary_dict[chr_id][window_id]:
                    # window_list[window_id] += line
                    last_unknown_end = this_line[unknown_side]['end']
                else:
                    boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))
                    window_id += 1
                    if len(window_list) <= window_id:
                        window_list.append('')
            else:
                # First judge whether this is a continuing block
                # Then judge which side is known side
                # But first finish the unfinished work
                boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))
                if this_line[known_side]['chr'] not in boundary_dict:
                    # from A188 Mo17 to B73 Mo17
                    (known_side, unknown_side) = (unknown_side, known_side)
                if this_line[known_side]['chr'] not in boundary_dict:
                    logging.debug("Error: No known window offset in both columns")
                    exit(1)
                boundary_dict[this_line[unknown_side]['chr']] = []
                window_id = 0
                # window_list[window_id] += line
            if offset == 'T':
                if window_id > 0:
                    loff_set = boundary_dict[line_list[0]][window_id - 1]
                    roff_set = boundary_dict[line_list[5]][window_id - 1]
                else:
                    loff_set = 0
                    roff_set = 0
                line_list[3] = str(int(line_list[3]) - int(loff_set))
                line_list[2] = str(int(line_list[2]) - int(loff_set))
                line_list[8] = str(int(line_list[8]) - int(roff_set))
                line_list[7] = str(int(line_list[7]) - int(roff_set))
                line = "\t".join(line_list).rstrip() + "\n"
            window_list[window_id] += line
            last_left_chr = line_list[0]
            last_right_chr = line_list[5]
        # After loop end, appending last window offset
        boundary_dict[this_line[unknown_side]['chr']].append(int(last_unknown_end))

    for wd in range(0, len(window_list)):
        out_paf = '{}.{}.paf'.format(paf_file, wd)
        out_bed = '{}.{}.bed'.format(paf_file, wd)
        out_gff = '{}.{}.gff'.format(paf_file, wd)
        out_fake_fa = '{}.{}.fa'.format(paf_file, wd)
        out_fai = '{}.{}.fa.fai'.format(paf_file, wd)
        with open(out_paf, 'w') as fh:
            fh.write(window_list[wd])
        with open(out_bed, 'w') as fh:
            buffer = ''
            for k in boundary_dict:
                if wd == 0:
                    start = 0
                else:
                    start = boundary_dict[k][wd - 1]
                # logging.debug(k)
                # logging.debug(wd)
                # logging.debug(boundary_dict[k])
                try:
                    buffer += "{}\t{}\t{}\n".format(k, start, boundary_dict[k][wd])
                except IndexError:
                    logging.debug([k, start, wd])
            fh.write(buffer)
        sh('bedtools intersect -a {} -b {} -wb |cut -f4,5,6,7,8,9 |sort -k1,1 -k2,2n |uniq >{} '.format(
            out_bed, bed_file, out_bed + 'ist'))

        if offset == 'T':
            intersect_bed = Bed(out_bed + 'ist')
            with open(out_fai, 'w') as fh:
                for k in boundary_dict:
                    if wd == 0:
                        start = 0
                    else:
                        start = boundary_dict[k][wd - 1]
                    intersect_bed.change_offset(k, start)
                    try:
                        chr_size = boundary_dict[k][wd] - start
                        fh.write('{0}\t{1}\t{2}\t{3}\n'.format(k, chr_size, start, boundary_dict[k][wd]))
                    except IndexError:
                        logging.debug([k, start, wd])
                intersect_bed.write(out_bed + 'ist')
            sh('touch {}'.format(out_fake_fa))
        bed_to_gff(out_bed + 'ist', out_gff)
        # debug
        # break
    return 0


def lastz_to_mosaic(lastz_file=None):
    """
    %s lastz.txt(sort -k1,1) > lastz.mosaic.txt
    :param lastz_file: Input alignment file
    :return: STDOUT
    """
    bedpe = BedPE(lastz_file, type='lastz')
    bedpe.get_mosaic()


def mosaic_ratio(fai=None, stat=None):
    """
    Calculate mosaic ratio based on mosaic region size and chromosome size
    :param fai: chromosome size
    :param stat: mosaic ratio (format: left_mosaic_size\tright_mosaic_size)
    :return:
    """
    chr_size = {}
    with open(fai) as fh:
        for line in fh:
            mylist = line.split()
            chr_size[mylist[0]] = int(mylist[1])
    with open(stat) as fh:
        line = fh.readline()
        (left_mosaic_size, right_mosaic_size) = line.strip().split()
    # PH207.genome.chr10.W22.genome.chr10.syri.out.mosaic.stat
    mylist = parse("{}.genome.chr{}.{}.genome.{}.syri.out.mosaic.stat", stat)
    tag = "{}-{}-{}".format(mylist[0], mylist[2], mylist[1])
    left_chr = "{}_{}".format(mylist[0], mylist[1])
    right_chr = "{}_{}".format(mylist[2], mylist[1])
    ratio = (int(left_mosaic_size) + int(right_mosaic_size)) / (chr_size[left_chr] + chr_size[right_chr])
    print("{}\t{}\t{}\t{}\t{}\t{:.1%}".format(tag, left_mosaic_size, right_mosaic_size,
                                              chr_size[left_chr], chr_size[right_chr], ratio))


def mosaic_ratio2(mosaic=None, synal=None):
    r"""
    Calculate mosaic ratio based on mosaic region size and chromosome size
    :param mosaic_file:
    :param synal_file:
    :return:
        xx-yy: mosaic%
        yy-xx: mosaic%
    """
    mosaic_left = 0
    mosaic_right = 0
    synal_left = 0
    synal_right = 0
    with open(mosaic) as fh:
        for line in fh:
            #A188_1  3453829 3454224 B73_1   1757    1996
            mylist = line.strip().split()
            mosaic_left += int(mylist[2]) - int(mylist[1]) + 1
            mosaic_right += int(mylist[5]) - int(mylist[4]) + 1
    with open(synal) as fh:
        for line in fh:
            #A188_1  5430894 5450765 -       -       B73_1   198985  218839  SYNAL9  SYN4    SYNAL   -
            mylist = line.strip().split()
            synal_left += int(mylist[2]) - int(mylist[1]) + 1
            synal_right += int(mylist[7]) - int(mylist[6]) + 1
    #prefix = get_prefix(mosaic)
    prefix = mosaic.replace('.genome', '')
    prefix = prefix.replace('.syri.out.SYNAL.mosaic', '')
    prefix = prefix.split('.')
    print("{0}-{2}\t{1}\n{2}-{0}\t{3}".format(prefix[0], mosaic_left/(mosaic_left+synal_left),
                                            prefix[1], mosaic_right/(mosaic_right+synal_right)))


def chromosome_level_ratio(stat=None):
    """
    Run this function after
    `$for i in *stat; do python -m iga.project.a188 mosaic_ratio total.fai $i >> mosaic_ratio_stat.txt; done`
    Input eg:
        A188-B73-10     66760136        70522917        146705344       150982314       46.1%
        A188-Mo17-10    70611181        72223243        146705344       149041351       48.3%
        A188-PH207-10   76993877        74673644        146705344       145727586       51.9%
        A188-SK-10      73901986        75129013        146705344       148363080       50.5%
    :param stat: mosaic_ratio_stat.txt
    :return:
    """
    cmp_dict1 = defaultdict(int)
    cmp_dict2 = defaultdict(int)
    cmp_dict3 = defaultdict(int)
    cmp_dict4 = defaultdict(int)
    with open(stat) as fh:
        for line in fh:
            (tag, l_mosaic_size, r_mosaic_size, l_chr_size, r_chr_size, ratio) = line.split()
            (l_id, r_id, chr_id) = tag.split('-')
            chr_tag = l_id + '-' + r_id
            cmp_dict1[chr_tag] += int(l_mosaic_size)
            cmp_dict2[chr_tag] += int(r_mosaic_size)
            cmp_dict3[chr_tag] += int(l_chr_size)
            cmp_dict4[chr_tag] += int(r_chr_size)
    for k in cmp_dict1:
        ratio = (cmp_dict1[k] + cmp_dict2[k]) / (cmp_dict3[k] + cmp_dict4[k])
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2%}".format(k, cmp_dict1[k], cmp_dict2[k],
                                                        cmp_dict3[k], cmp_dict4[k], ratio))


def split_by_tag(table=None, column=None):
    """
    split table by column
    :param table: 0-based
    :param colum:
    :return:
    """
    table_dict = defaultdict(str)
    column = int(column)
    with open(table) as fh:
        for line in fh:
            mylist = line.strip().split()
            table_dict[mylist[column]] += line
    for k in table_dict.keys():
        new_name = table + '.' + k
        with open(new_name, 'w') as fh:
            fh.write(table_dict[k])


def bedpe_intersect(bed1=None, bed2=None):
    """
    intersection of two bed pe file, like A188.B73.mosaic  and B73.Mo17.mosaic, default in left to right order
    :param bed1: A188.B73.mosaic
    :param bed2: B73.Mo17.mosaic
    :return:
    """
    bed1_buf = BedPE(bed1)
    bed2_buf = BedPE(bed2)

    bed1_out = bed1 + '.cut'
    bed2_out = bed2 + '.cut'

    bed1_buf.write_to_table(table=bed1_out, right_only=True)
    bed2_buf.write_to_table(table=bed2_out, left_only=True)

    sh("bedtools intersect -a {} -b {} > {}".format(bed1_out, bed2_out, bed1_out + bed2_out))

    intesect = Bed(bed1_out + bed2_out)

    intersect_sum = intesect.sum_size()
    print(intersect_sum)
    return intersect_sum


def bed_size(bed=None):
    """
    return bed size
    :param bed:
    :return:
    """
    bed = Bed(bed)
    print(bed.sum_size())
    return bed.sum_size()


def intersect_bedpe(bed1=None, bed2=None, type='ll', output='lr', wobble=100):
    """
    Get mosaic intersections
    :param bed1:
    :param bed2:
    :param type: could be ll lr rl rr. defining which side to be compared
    :param output: could be 'lr, l, r',defining pring both beds, left bed or right bed only
    :param wobble: the number of nucleotide that can be wobbled for the boundaries of mosaic region
    :return:
    """
    bed1_obj = BedPE(bed1)
    bed2_obj = BedPE(bed2)

    # self.bedpe_db[left_chr].append(this_lp)
    result = ''

    for chr in bed1_obj.bedpe_db:
        for bedpe_loci in bed1_obj.bedpe_db[chr]:
            search_result = bed2_obj.exists(bedpe_loci, wobble=int(wobble), type=type, output=output)
            if search_result != '':
                result += search_result
    print(result, end='')
    return result


def breakpoint_hotspot(bed1=None, bed2=None, type='ll', output='lr', wobble=100):
    """
    Get breakpoint hotspot between two bedPE files
    :param bed1:
    :param bed2:
    :param type: could be ll lr rl rr. defining which side to be compared
    :param output: could be 'lr, l, r',defining pring both beds, left bed or right bed only
    :param wobble: the number of nucleotide that can be wobbled for the boundaries of mosaic region
    :return:
    """
    bed1_obj = BedPE(bed1)
    bed2_obj = BedPE(bed2)

    # self.bedpe_db[left_chr].append(this_lp)
    start_all = ''
    end_all = ''

    for chr in bed1_obj.bedpe_db:
        for bedpe_loci in bed1_obj.bedpe_db[chr]:
            [start, end] = bed2_obj.hotspot(bedpe_loci, wobble=int(wobble), type=type, output=output)
            start_all += start
            end_all += end
    with open(bed1 + bed2 + '.start', 'w') as fh:
        fh.write(start_all)
    with open(bed1 + bed2 + '.end', 'w') as fh:
        fh.write(end_all)
    return 0


def bed_to_gff(bed=None, output=''):
    """
    Transform bed to gff format
    :param bed:
    :return:
    """
    bed_obj = Bed(bed)
    gff_obj = GFF()
    for k in bed_obj.bed_list:
        # logging.debug(k.name)
        if k.strand == '.':
            k.strand = '+'
        gff_obj.append(chr=k.chr, start=k.start, end=k.end, strand=k.strand, name=k.name)
    if output == '':
        gff_obj.print_out()
    else:
        with open(output, 'w') as fh:
            fh.write(gff_obj.to_str())


def annotate_block(anchor_file=None, qbed='', sbed=''):
    """
    alias: Ortho to bedpe
    Add ID, qstart, qend, sstart, send information for each syntenic block
    :param anchor_file: .anchor file generated by jcvi.compara.ortholog
    :param qbed: query bed file
    :param sbed: subject bed file
    :return:
    """
    qbed_obj = Bed(qbed)
    sbed_obj = Bed(sbed)
    # The ID for syntenic block
    syn_id = -1
    # The dictionary storing content for each syntenic block
    syn_content = defaultdict(str)
    syn_start = {}
    syn_end = {}
    with open(anchor_file) as fh:
        for line in fh:
            if line.startswith('#'):
                # starting of a syntenic block
                syn_id += 1
            else:
                (qgene, sgene, score) = line.rstrip().split()
                qloci = qbed_obj.select_name(qgene, format='loci')
                sloci = sbed_obj.select_name(sgene, format='loci')
                syn_content[syn_id] += qloci.get_line().strip() + "\t" + sloci.get_line()
                if syn_id not in syn_start:
                    syn_start[syn_id] = [qloci, sloci]
                if True:
                    syn_end[syn_id] = [qloci, sloci]
    for i in range(0, syn_id + 1):
        header = '#'
        qry_chr = syn_start[i][0].chr
        qry_start = syn_start[i][0].start
        qry_end = syn_end[i][0].end
        sub_chr = syn_start[i][1].chr
        sub_start = syn_start[i][1].start
        sub_end = syn_end[i][1].end
        strand = (sub_end - sub_start) * (qry_end - qry_start)
        if strand < 0:
            strand = '-'
        else:
            strand = '+'
        block_name = "Block{0}".format(i)
        tmp = "{}\t" * 8
        tmp = tmp.rstrip()
        header += tmp.format(qry_chr, qry_start, qry_end, sub_chr, sub_start, sub_end,
                             block_name, strand)
        content = syn_content[i]
        print(header + "\n" + content, end='')


def join_adjacent_bed(bed=None):
    """
    %s join_adjacent_bed manual_reviewed_subgenome.bed > manual_reviewed_subgenome.joined.bed
    same name and same chr are joined
    Input:
    1	6275775	6999144	ALN1
    1	7257371	15121456	ALN1
    Output:
    1 6275775 15121456 ALN1
    :param bed:
    :return:
    """
    bed_read = Bed(bed)
    loci = copy.copy(bed_read.bed_list[0])
    for i in range(1, len(bed_read.bed_list)):
        loci_i = bed_read.bed_list[i]
        if loci.chr == loci_i.chr and loci.name == loci_i.name:
            loci.start = min(loci.start, loci_i.start)
            loci.end = max(loci.end, loci_i.end)
        else:
            print(loci.get_line(), end='')
            loci = copy.copy(loci_i)
    print(loci.get_line(), end='')


def sums(input_list):
    result = 0
    for l in input_list:
        result += int(l)
    return result


def join_contiguous_bed(bed=None):
    """
    %s join_contiguous_bed A188.1bp.B73.Mo17.bed > A188.1bp.B73.Mo17.bed.joined.bed
    same chr and same haplotype are joined, with syntenic status tagged at last position
    Input (chr start end syn_B73 syn_Mo17 syn_W22...):
    chr1	0   1	0   0   0   0   0
    chr1	1   2   0   0   0   0   0
    chr1   2   3   1   0   0   0   0
    chr1   3   4   0   1   0   1   0
    chr1   4   5   0   0   1   1   0
    chr1   5   6   1   1   1   1   1
    chr1   6   7   1   1   1   1   1
    Output:
    chr1   0    2   SMR
    chr1    2   5   OMR
    chr1    5   7   CSR
    :param bed:
    :return:
    """
    prev_chrid = ''
    prev_start = 0
    prev_end = 0
    tag_dict = {0: "MOSAIC", 1: "OMOSAIC", 2: "OMOSAIC", 3: "OMOSAIC", 4: "OMOSAIC", 5: "CSYN"}
    num_comparison = 5
    prev_tag = ''
    prev_line = ''
    with open(bed) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            (chrid, start, end) = mylist[0:3]
            taglist = mylist[3:]
            for i in range(0, len(taglist)):
                taglist[i] = str(min(1, int(taglist[i])))
            # this_sum = min(num_comparison, sums(taglist))
            this_sum = sums(taglist)
            try:
                this_tag = "\t".join(taglist + [tag_dict[this_sum]])
                # this_tag = tag_dict[this_sum]
            except KeyError:
                logging.debug(this_sum)
                logging.debug(tag_dict)
                logging.debug(line)
                exit(1)
            if prev_chrid != '' and prev_chrid == chrid and prev_tag == this_tag:
                prev_end = end
                prev_line = "\t".join([prev_chrid, prev_start, prev_end, prev_tag]) + "\n"
            else:
                print(prev_line, end='')
                prev_tag = this_tag
                prev_start = start
                prev_chrid = chrid
                prev_end = end
                prev_line = "\t".join([prev_chrid, prev_start, prev_end, prev_tag]) + "\n"
        print(prev_line, end='')


def calc_percent_bed_intersect(bedwo=None):
    """
    %s A188_B73.total.SYN.mosaic.B73.bed.leaf_atac.wo > A188_B73.total.SYN.mosaic.B73.bed.leaf_atac.wo.perc
    Get relative position of intersection on Feature A
    Input:
        B73_1   2097669 2104245 NOTSYN26        B73_1   2104191 2105241 genic   54
        B73_1   2566205 2629140 NOTSYN34        B73_1   2628589 2628839 proximal        250
        B73_1   2636499 2746464 NOTSYN35        B73_1   2635661 2636711 genic   212
        B73_1   2892522 2894589 NOTSYN39        B73_1   2894533 2894758 proximal        56
    Output:
         B73_1   2097669 2104245 NOTSYN26        B73_1   2104191 2105241 genic   54      99.18%  100.00%
         B73_1   2566205 2629140 NOTSYN34        B73_1   2628589 2628839 proximal        250     99.12%  99.52%
         B73_1   2636499 2746464 NOTSYN35        B73_1   2635661 2636711 genic   212     0.00%   0.19%
         B73_1   2892522 2894589 NOTSYN39        B73_1   2894533 2894758 proximal        56      97.29%  100.00%
         B73_1   2973072 2974460 NOTSYN41        B73_1   2974176 2974451 genic   275     79.54%  99.35%
    :param bedwo:
    :return:
    """
    min_mosaic_size = 1000
    with open(bedwo) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            region_size = int(mylist[2]) - int(mylist[1])
            if region_size < min_mosaic_size:
                continue
            offset = int(mylist[1])
            rel_start = int(mylist[5]) - offset
            rel_end = int(mylist[6]) - offset
            start_percent = max(0, rel_start / region_size)
            end_percent = min(1, rel_end / region_size)
            print(line.rstrip() + "\t{:.2%}\t{:.2%}".format(start_percent, end_percent))


def percent_to_range(percent_file=None):
    r"""
    pass [a,b] to a a+1..b like format, for hist use
    Input:
        97.18% 100.00%
    Output:
        97
        98
        99
        100
    :param percent_file:
    :return:
    """
    import re
    with open(percent_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            for i in range(0, len(mylist)):
                mylist[i] = re.sub(r'\..*', '', mylist[i])
                mylist[i] = int(mylist[i])
            for a in range(mylist[0], mylist[1] + 1):
                print(a)


def percent_to_range2(percent_file=None):
    r"""
    pass [a,b] to a a+1..b like format, for hist use
    :param percent_file:
    :return:
    """
    # Input:
    #     B73_1   2431547 2462809 1_NOTSYNAL117   B73_1   2435718 2435990 Mo17rnd-5_family-12     0       -       272     13.34% 14.21%
    # Output:
    #     B73_1   2431547 2462809 1_NOTSYNAL117   B73_1   2435718 2435990 Mo17rnd-5_family-12     0       -       272     13%
    #     B73_1   2431547 2462809 1_NOTSYNAL117   B73_1   2435718 2435990 Mo17rnd-5_family-12     0       -       272     14%
    #     98
    #     99
    #     100
    import re
    with open(percent_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            start_perc = round(float(mylist.pop(-2).replace('%', '')))
            end_perc = round(float(mylist.pop(-1).replace('%', '')))
            # logging.warning(start_perc)
            # logging.warning(end_perc)
            # exit(1)
            new_line = "\t".join(mylist)
            for i in range(start_perc, end_perc + 1):
                print("{}\t{}".format(new_line, i))
                # mylist[i] = re.sub(r'\..*', '', mylist[i])
                # mylist[i] = int(mylist[i])
            # for a in range(mylist[0], mylist[1] + 1):
            #     print(a)


def breakpoint_screen(depth=None, highcutoff=100, lowcutoff=5):
    """
    %s A.depth.gz > A.breakpoint.txt
    Read depth file, and identify waterfall site that had a significant coverage drop
    :param highcutoff: larger than that is normal
    :param lowcutoff: lower than that is waterfall regions
    :param depth:
    :return:
    """
    import gzip
    higher_cutoff = highcutoff
    low_cutoff = lowcutoff
    with gzip.open(depth) as fh:
        prev_depth = 0
        prev_loci = 0
        prev_chr = ''
        buffer = None
        start_flag = False
        for line in fh:
            (chr_id, loci, depth) = line.decode().rstrip().split()
            if int(loci) % 10000000 == 0:
                logging.debug("Running on chr {} loci {}".format(chr_id, loci))
            depth = int(depth)
            if chr_id != prev_chr:
                if buffer is not None:
                    print(buffer.get_line(), end='')
            elif prev_depth > higher_cutoff and depth < low_cutoff:
                start_flag = 1
                buffer = Loci(chr_id, loci, loci, '.', '.', '.')
            elif start_flag and depth < low_cutoff:
                buffer.end = loci
            elif start_flag and depth > low_cutoff:
                print(buffer.get_line(), end='')
                buffer = None
                start_flag = False
            prev_depth = depth
            prev_loci = loci
            prev_chr = chr_id


def breakpoint_screen2(bam=None, add_name='F'):
    """
    %s test.bam > test.bam.breakpoint.txt
    STDOUT: depth of reads at each loci
    add_name: [T|F]
    Output: bam + '.summary' sum of split-reads at 10-flanking of mosaic region boundaries
    :return:
    """
    import pysam
    samfile = pysam.AlignmentFile(bam, "r")
    reads_all = samfile.fetch()
    buf = defaultdict(int)
    for read in reads_all:
        ###pysam's coordinate [0-based, 0-based), like bed, so have the following modifications
        # if type(read.reference_id) == int:
        #     read.reference_id += 1
        read.reference_start += 1
        # pysam 0.8.3, Future changed to latest version with reference_name
        try:
            reference_name = read.reference_name
        except AttributeError:
            reference_name = samfile.references[read.reference_id]
        ###
        if read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
            print_buff = "{}\t{}\t{}".format(reference_name, read.reference_start, "Start")
            if add_name == 'T':
                print_buff += "\t{}".format(read.qname)
            print(print_buff)
            buf["{}\t{}\t{}".format(reference_name, read.reference_start, "Head")] += 1
        if read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5:
            print_buff = "{}\t{}\t{}".format(reference_name, read.reference_end, "End")
            if add_name == 'T':
                print_buff += "\t{}".format(read.qname)
            print(print_buff)
            buf["{}\t{}\t{}".format(reference_name, read.reference_end, "Tail")] += 1
    with open(bam + '.summary', 'w') as fh:
        for i in buf:
            fh.write("{}\t{}\n".format(i, buf[i]))


def add_depth_to_mosaic(mosaic_bedpe=None, bkptsum_l=None,
                        bkptsum_r=None):
    """
    All is 1-based
    :param mosaic_bedpe:
    :param bkptsum: generated by brekpoint_screen2
    :param depth:
    :param with_read_name: whether there is readname in bkptsum file
    :return:
    """
    mbe = BedPE(mosaic_bedpe)
    breakpoint_coverage_cutoff = 5
    bkptdb = defaultdict(int)
    bkptdb_right = defaultdict(int)
    with open(bkptsum_l) as fh:
        # 1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            (chrid, loci, croptype, coverage) = mylist
            if int(coverage) < breakpoint_coverage_cutoff:
                continue
            bkptdb["_".join([chrid, loci])] += int(coverage)
    with open(bkptsum_r) as fh:
        # 1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            try:
                (chrid, loci, croptype, coverage) = mylist
            except ValueError:
                logging.debug(line)
            if int(coverage) < breakpoint_coverage_cutoff:
                continue
            bkptdb_right["_".join([chrid, loci])] += int(coverage)
    breakpoint_flanking = 10
    for lchr in mbe.bedpe_db:
        logging.debug("left chromosome is {}".format(lchr))
        for lpe in mbe.bedpe_db[lchr]:
            lchr = re.sub(r'.*_', '', lpe.left.chr)
            rchr = re.sub(r'.*_', '', lpe.right.chr)
            # logging.debug('searching {}'.format(lpe.get_line()))
            etc = ''
            for left in (lpe.left.start, lpe.left.end):
                sum = 0
                for i in range(left - breakpoint_flanking, left + breakpoint_flanking):
                    try:
                        # logging.debug(bkptdb["_".join([lchr, i])])
                        # logging.debug(sum)
                        sum += bkptdb["_".join([lchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            for right in (lpe.right.start, lpe.right.end):
                sum = 0
                for i in range(right - breakpoint_flanking, right + breakpoint_flanking):
                    try:
                        sum += bkptdb_right["_".join([rchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            if etc != '':
                print(lpe.get_line().rstrip() + "\t" + etc.rstrip())


def add_rdname_to_mosaic(mosaic_bedpe=None, bkptsum_l=None,
                         bkptsum_r=None):
    """
    Add read name to mosaic file
    :param mosaic_bedpe:
    :param bkptsum: generated by brekpoint_screen2
    :param depth:
    :return:
    """
    mbe = BedPE(mosaic_bedpe)
    bkptdb = defaultdict(str)
    bkptdb_right = defaultdict(str)
    with open(bkptsum_l) as fh:
        # 1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            (chrid, loci, croptype, rdname) = mylist
            bkptdb["_".join([chrid, loci])] += rdname + ";"
    with open(bkptsum_r) as fh:
        # 1       37      Head    1
        for line in fh:
            mylist = line.rstrip().split()
            try:
                (chrid, loci, croptype, rdname) = mylist
            except ValueError:
                logging.debug(line)
            bkptdb_right["_".join([chrid, loci])] += rdname + ";"
    breakpoint_flanking = 10
    for lchr in mbe.bedpe_db:
        logging.debug("left chromosome is {}".format(lchr))
        for lpe in mbe.bedpe_db[lchr]:
            lchr = re.sub(r'.*_', '', lpe.left.chr)
            rchr = re.sub(r'.*_', '', lpe.right.chr)
            # logging.debug('searching {}'.format(lpe.get_line()))
            etc = ''
            for left in (lpe.left.start, lpe.left.end):
                sum = ''
                for i in range(left - breakpoint_flanking, left + breakpoint_flanking):
                    try:
                        # logging.debug(bkptdb["_".join([lchr, i])])
                        # logging.debug(sum)
                        sum += bkptdb["_".join([lchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            for right in (lpe.right.start, lpe.right.end):
                sum = ''
                for i in range(right - breakpoint_flanking, right + breakpoint_flanking):
                    try:
                        sum += bkptdb_right["_".join([rchr, str(i)])]
                    except KeyError:
                        pass
                etc += str(sum) + "\t"
            if etc != '':
                print(lpe.get_line().rstrip() + "\t" + etc.rstrip())


def bedpe_to_ggplot(bedpe=None):
    """
    %s bedpe > bedpe.ggplot_in
    :param bedpe: like mosaic.bedpe generated by synal_to_mosaic
    :return:
    """
    from random import randint
    read = BedPE(bedpe)
    output_template = "{}\t" * 6
    output_template = output_template.rstrip()
    rand_max = 12000000
    # visit dict
    for chr_id in sorted(read.bedpe_db.keys()):
        chr_id_format = re.sub('.*_', '', chr_id)
        iter = 1
        ymin = 0
        # visit list
        for lp in read.bedpe_db[chr_id]:
            if 'NOT' in lp.left.name:
                sample1_id = re.sub('_.*', '', lp.left.chr)
                sample2_id = re.sub('_.*', '', lp.right.chr)
                yoffset = randint(0, rand_max)
            else:
                sample1_id = "SYN"
                sample2_id = "SYN"
                yoffset = 0
            x1size = lp.left.get_size(bed_format=False)
            x2size = lp.right.get_size(bed_format=False)
            if x1size > 1:
                x1min = iter
                x1max = iter + x1size - 1
                y1max = x1size + yoffset
                y1min = ymin + yoffset
                print(output_template.format(chr_id_format, x1min, x1max, y1min, y1max, sample1_id))
            if x2size > 1:
                x2min = iter
                x2max = iter + x2size - 1
                y2max = x2size * - 1 - yoffset
                y2min = ymin - yoffset
                print(output_template.format(chr_id_format, x2min, x2max, y2min, y2max, sample2_id))
            iter += max(x1size, x2size)

    # import gzip
    # with gzip.open(depth_l) as fh:
    #     for line in fh:
    #         (chr_id, loci, depth) = line.decode().rstrip().split()


def filter_bam_by_reads(reads=None, bam=None):
    """
    %s reads_name_file bam > reads_name_file.sam
    :param reads: file contain reads name
    :param bam: the bam file
    Output: STDOUT
    0422.
    :return:
    """
    import pysam
    reads_dt = {}
    with open(reads) as fh:
        for line in fh:
            reads_dt[line.rstrip()] = 1
    samfile = pysam.AlignmentFile(bam, "r")
    reads_all = samfile.fetch()
    buf = defaultdict(int)
    for read in reads_all:
        ###pysam's coordinate [0-based, 0-based), like bed, so have the following modifications
        # if type(read.reference_id) == int:
        #     read.reference_id += 1
        if read.qname in reads_dt:
            # print(dir(read))
            try:
                print(read.to_string())
            except AttributeError:
                # Deprecated in newer pysam
                print(read.tostring(htsfile=''))


def mtei_union(TIP_table=None, difftag=''):
    """
    :param TIP_table:
    :return:
    Input eg:
    workdir_TIP_chr10.combine.long
    B73_10  272712  272940  Mo17_10 150291  173980  10_NOTSYNAL13   .       +       +       10_NOTSYNAL13   B73_10  260657  276914  DHH00001        0       +       |       Mo17_10 128780  223337  RLX27954        0
    """
    len_left = {}
    len_right = {}
    te_left = {}
    te_right = {}
    with open(TIP_table) as fh:
        for line in fh:
            mylist = line.split('\t')
            len_left[mylist[6]] = set(range(int(mylist[1]), int(mylist[2]) + 1))
            len_right[mylist[6]] = set(range(int(mylist[4]), int(mylist[5]) + 1))
            if mylist[12] == '':
                continue
            if not(difftag == '' or mylist[17] == difftag):
                continue
            if mylist[6] not in te_left:
                te_left[mylist[6]] = set(range(int(mylist[12]), int(mylist[13]) + 1))
            else:
                te_left[mylist[6]] = te_left[mylist[6]].union(set(range(int(mylist[12]), int(mylist[13]) + 1)))
            if len(mylist) < 19 or mylist[19] == '':
                continue
            if mylist[6] not in te_right:
                te_right[mylist[6]] = set(range(int(mylist[19]), int(mylist[20]) + 1))
            else:
                te_right[mylist[6]] = te_left[mylist[6]].union(set(range(int(mylist[19]), int(mylist[20]) + 1)))
    for k in len_left:
        a = len(len_left[k])
        a2 = len(len_right[k])
        if k in te_left:
            b = len(len_left[k].intersection(te_left[k]))
            try:
                c = b / a
            except ZeroDivisionError:
                c = 0
        else:
            b = 0
            c = 0
        if k in te_right:
            b2 = len(len_right[k].intersection(te_right[k]))
            try:
                c2 = b2 / a2
            except ZeroDivisionError:
                c2 = 0
        else:
            b2 = 0
            c2 = 0
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(k, a, b, c, a2, b2, c2))


# 0 species [Rice|Maize|others]
# 1 cds
# 2 GENOME
# 3 threads
edta_sh = """
~/bin/EDTA/EDTA.pl  --species {0} --cds {1} --genome {2} --anno 1 --threads {3}
# --curatedlib maizeTE02052020
"""


def edta(genome=None, cds=None, species='others', threads=40):
    """
    EDTA wrapper
    :param cds: cds.fasta
    :param species: default others, could also be Rice or Maize
    :param threads: threads, default 40 running, 10 submitting
    :param genome: genome.fasta
    :return:
    """
    cmd = "conda activate EDTA"
    cmd += edta_sh.format(species, cds, genome, threads)
    qsub(cmd, cpus=5, name='EDTA')


# The following three function is run sequentially
#
# CDS=
# python -m iga.project.a188 calcKs_OF OrthoFinder/Results_Jun10/Single_Copy_Orthologue_Sequences $CDS
# python -m iga.project.a188 collect_calcKs_OF OrthoFinder/Results_Jun10/Single_Copy_Orthologue_Sequences > single_copy_cds.aln"
# python -m iga.project.a188 get_dist single_copy_cds.aln
#
# Output Dist is single_copy_cds.aln.dist[012345]

calcKs_OF_sh = r"""
cd {0}
# export MAX_N_PID_4_TCOFFEE=1000
for g in OG*
do
    #t_coffee ${g}.fa -mode fmcoffee  > {g}.pep.aln &&  pal2nal.pl ${g}.pep.aln ${g}.cds   >${g}.paml_aln
    muscle -clw -in ${g}.fa -out {g}.pep.aln &&  pal2nal.pl ${g}.pep.aln ${g}.cds   >${g}.paml_aln
done
"""


def calcKs_OF(Single_Copy_Orthologue_Sequences=None, total_cds=None, debug='F', submit='T'):
    """
    :param Single_Copy_Orthologue_Sequences: The directory of Single_Copy_Orthologue_Sequences (in orthofinder2)
    :return: group.paml_aln which is coding sequences alignment guided by protein alignment
    """
    count = 0
    cds_dict = SeqIO.to_dict(SeqIO.parse(total_cds, "fasta"))
    os.chdir(Single_Copy_Orthologue_Sequences)
    for g in os.listdir('.'):
        if 'OG' not in g or not g.endswith('.fa'):
            continue
        if debug == 'T':
            logging.debug(g)
        if not os.path.exists(g+'.cds'):
            pep_dict = SeqIO.to_dict(SeqIO.parse(g, "fasta"))
            with open(g + '.cds', 'w') as fh:
                for p in pep_dict:
                    fh.write(cds_dict[p].format('fasta'))
        cmd = "t_coffee {0}  > {0}.aln &&  pal2nal.pl {0}.aln {0}.cds   >{0}.paml_aln".format(g)
        # cmd = "muscle -clw -in {0} -out {0}.aln && sed -i 's/MUSCLE (3.7)/CLUSTAL W/' {0}.aln &&  pal2nal.pl {0}.aln {0}.cds   >{0}.paml_aln".format(g)
        if submit == 'T':
            qsub(cmd, name=g, normal='T')
        else:
            sh(cmd)
        if debug == 'T':
            count += 1
            if count > 5:
                break
    return 0


collect_calcKS_OF = r"""
cat {0}/*paml_aln |grep -v CLUSTAL | sed "s/OsR498.*\s/indica    /; s/OS0.*\s/japonica  /">$TREE_DIR/total.aln

#sed -i '1iCLUSTAL W multiple sequence alignment' $TREE_DIR/total.aln

trimal -in total.aln > total.trimal.aln

distmat -nucmethod 4  -sequence $TREE_DIR/total.trimal.aln -outfile $TREE_DIR/total.aln.dist
"""


def collect_calcKs_OF(Single_Copy_Orthologue_Sequences=None, species_name_col=1):
    """
    :param Single_Copy_Orthologue_Sequences: directory of that directory
    :param species_name: if zemay_A188_oGxxx, then species name A188 is 1 (split('_')[1])
    :return: STDOUT, need to redirect to a file
    """
    species_name_col = int(species_name_col)
    print('CLUSTAL W multiple sequence alignment')
    for g in os.listdir(Single_Copy_Orthologue_Sequences):
        if g.endswith('paml_aln'):
            g_abs_path = os.path.join(Single_Copy_Orthologue_Sequences, g)
            if os.path.getsize(g_abs_path) == 0:
                continue
            sh('trimal -nogaps -in {0} -out {0}.format'.format(g_abs_path))
            if not os.path.exists(g_abs_path + '.format'):
                #there is a bug in pal2nal
                #sometime will output
                #
                # solyc_BGV006775_Solyc07g150117.1.1     ATCTAATGATTCTGGACGTAATCCTGGACGTGAAGA
                # solyc_BGV006865_Solyc07g150117.1.1     CTAATGATTCTGGACGTAATCCTGGACGTGAAG---AAT
                # solyc_Brandywine_Solyc07g150117.1.1    ATCTAATGATTCTGGACGTAATCCTGGACGTGAAGA
                # solyc_Floradade_Solyc07g150117.1.1     ATC4TAATGATTCTGGACGTAATCCTGGACGTGAAGA
                # solyc_M82_Solyc07g150117.1.1           ATCTAATGATTCTGGACGTAATCCTGGACGTGAAGA
                #
                # uder/Results_Jun10/Single_Copy_Orthologue_Sequences/OG0019286.fa.paml_aln (END)
                continue
            with open(g_abs_path + '.format') as fh:
                fh.readline()
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(' ') or line == '':
                        print('')
                        continue
                    else:
                        try:
                            (id, seq) = line.split()
                        except ValueError:
                            logging.debug(line)
                            continue
                        id_list = re.split(r'[-_]', id)
                        newid = id_list[species_name_col]
                        print("{: <30}{}".format(newid, seq))


get_dist_sh = """
trimal -nogaps -in {0} > {0}.trim
for i in `seq 0 5`
do
    distmat -nucmethod $i  -sequence {0}.trim -outfile {0}.trim.dist$i
done
"""


def get_dist(aln=None, qsub='F'):
    """
    :param aln:
    :return:
    """
    cmd = get_dist_sh.format(aln)
    if qsub == 'T':
        qsub(cmd, cpus=1, normal='T')
    else:
        sh(cmd)


if __name__ == "__main__":
    emain()
