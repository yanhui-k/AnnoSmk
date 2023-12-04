"""
GFF relevant utils
"""
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import chain
from statistics import mean

from parse import parse

from iga.apps.base import emain, sh

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


class Feat:
    r"""
    The feat data structure that is needed by GFF class, support:
    1. Get parent
    2. Get childs
    #chr01   .       gene    12132486        12138762        .       -       .       ID=CORNEG00007591;
    #Name=CORNE00007591-t5;Alias=maker-000023F|arrow_np1212-snap-gene-26.41
    """

    def __init__(self, gff_line=''):
        self.childs = []
        self.content = ''
        self.seqid = ''
        self.source = '.'
        self.type = ''
        self.start = None
        self.end = None
        self.score = '.'
        self.strand = None
        self.phase = '.'
        self.len = None
        self.attributes = None
        self.attr_dict = OrderedDict()
        self.parent = None
        self.ID = None

        if gff_line != "":
            self.content = gff_line.rstrip()
            mylist = self.content.split('\t')
            # Assign gff values by
            # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
            self.seqid = mylist[0]
            self.source = mylist[1]
            self.type = mylist[2]
            self.start = mylist[3]
            self.end = mylist[4]
            self.score = mylist[5]
            self.strand = mylist[6]
            self.phase = mylist[7]
            # Refactor affected vairables
            self.len = abs(int(self.end) - int(self.start)) + 1
            self.attributes = mylist[8]
            attr_list = self.attributes.rstrip(';').split(';')
            for a in attr_list:
                try:
                    (attr_key, attr_value) = parse("{}={}", a)
                except TypeError:
                    logger.error("Type Error on line {}, list {} and attribute: {}".format(self.attributes, attr_list, a))
                    continue
                self.attr_dict[attr_key] = attr_value
            if 'Parent' in self.attr_dict:
                self.parent = self.attr_dict['Parent']
            try:
                self.ID = self.attr_dict['ID']
            except KeyError:
                self.ID = None

    def get_parent(self):
        r"""
        Return the name of the parent
        :return:
        """
        return self.parent

    def add_child(self, child):
        r"""
        child is feat type
        :param child:
        :return:
        """
        if type(child) != Feat:
            raise TypeError
        self.childs.append(child)

    def get_all_child_feats(self, type=''):
        r"""
        return all lines that are descendants of this feat
        :param type:
        TODO: There i a bug here, will print gene mRNA CDS if use get_all_child_feats(CDS)
        :return:
        """
        if len(self.childs) == 0:
            if type == '' or self.type == type:
                # return when type is wild card or self.type equals specified type
                result = self.content + "\n"
        else:
            result = self.content + "\n"
            for i in self.childs:
                result += i.get_all_child_feats(type)
        return result

    def get_all_child_feats_obj(self, type=''):
        r"""
        return all descendants feats of this feat
        :param type:
        :return:
        """
        if type == '' or self.type == type:
            # return when type is wild card or self.type equals specified type
            result = [self]
        else:
            result = []
        if len(self.childs) == 0:
            pass
        else:
            for i in self.childs:
                result += i.get_all_child_feats_obj(type)
#        if any(isinstance(i, list) for i in result):
#            result = list(chain.from_iterable(result))
        return result

    def update_tag(self, tag, value):
        r"""
        Add or append a tag to this feat
        :param tag:
        :param value:
        :return:
        """
        self.attr_dict[tag] = value
        self.__refactor__()

    def delete_tag(self, tag):
        r"""
        Remove a tag
        :param tag:
        :return:
        """
        if tag == "ID":
            logger.error("Can't delete ID item")
            return 1
        if tag in self.attr_dict:
            del self.attr_dict[tag]
            self.__refactor__()
        return 0

    def __refactor__(self):
        r"""
        Update all fields
        :return:
        """
        self.attributes = ''
        for i in self.attr_dict:
            self.attributes += "{}={};".format(i, self.attr_dict[i])
        self.ID = self.attr_dict['ID']
        self.content = "\t".join([self.seqid, self.source, self.type, str(self.start), str(self.end),
                                  str(self.score), self.strand, str(self.phase), self.attributes])
        if 'Parent' in self.attr_dict:
            self.parent = self.attr_dict['Parent']
        self.len = abs(int(self.end) - int(self.start)) + 1
        return 0

    def new(self, chr, start, end, strand, name, type='CDS', score='.', source='.'):
        self.start = start
        self.end = end
        self.seqid = chr
        self.strand = strand
        self.type = type
        self.score = score
        self.source = source
        self.update_tag("ID", name)
        # logging.debug(name)
        # logging.debug(self.ID)
        # logging.debug(self.content)
        # exit(1)

    def print_all_childs(self):
        result = self.get_all_child_feats()
        print(result, end="")


class GFF:
    r"""
    GFF class that support:
    1. reading a GFF into memory
    2. adding feat by gene ID or transcript ID
    3. print out
    ...3. extracting feat by level or by transcript ID into tab delimited file
    """

    def __init__(self, filename=''):
        # Top level, which has no parent, usually gene type
        self.top_level_list = []
        # A dictionary store all type
        self.GFF_dict = OrderedDict()
        # cds name is same
        count_cds = defaultdict(int)
        if filename != '':
            with open(filename) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    feat = Feat(line)
                    if "Parent" not in feat.content:
                        # Top level
                        self.top_level_list.append(feat.attr_dict['ID'])
                        self.GFF_dict[feat.ID] = feat
                    else:
                        if feat.type == "CDS":
                            # Manage duplicate CDS
                            original_cds_name = feat.attr_dict['ID']
                            count_cds[original_cds_name] += 1
                            if count_cds[original_cds_name] > 1:
                                new_name = "{}:cds{}".format(original_cds_name, count_cds[original_cds_name])
                                feat.update_tag("ID", new_name)
                                if original_cds_name in self.GFF_dict:
                                    new_first_cds = "{}:cds{}".format(original_cds_name, 1)
                                    self.GFF_dict[new_first_cds] = self.GFF_dict[original_cds_name]
                                    self.GFF_dict[new_first_cds].update_tag("ID", new_first_cds)
                                    del self.GFF_dict[original_cds_name]
                                self.GFF_dict[new_name] = feat
                            else:
                                self.GFF_dict[feat.ID] = feat
                        else:
                            self.GFF_dict[feat.ID] = feat
                        parent = feat.parent
                        self.GFF_dict[parent].add_child(feat)

    def print_out(self):
        r"""
        print out gff to screen
        :return:
        """
        total_result = ''
        for k in self.top_level_list:
            result = self.GFF_dict[k].get_all_child_feats()
            print(result.rstrip())
        return 0

    def to_str(self):
        r"""
        print out gff to screen
        :return:
        """
        total_result = ''
        for k in self.top_level_list:
            result = self.GFF_dict[k].get_all_child_feats()
            total_result += result.rstrip() + "\n"
        return total_result

    def get_attr(self, attr=''):
        r"""
        Return attr as a dict, like extracting _AED from maker GFF
        :return:
        """
        result = OrderedDict()
        for k in self.GFF_dict:
            if attr in self.GFF_dict[k].attr_dict:
                result[k] = self.GFF_dict[k].attr_dict[attr]
        return result

    def longest_mRNA(self):
        """
        :return: longest_table and longest_gff
        """
        longest_table = ""
        longest_gff = ""
        for k in self.top_level_list:
            mRNA_list = self.GFF_dict[k].get_all_child_feats_obj('mRNA') + \
                        self.GFF_dict[k].get_all_child_feats_obj('transcript')
            longest = 0
            #such as tRNA genes
            non_gene_flag = False
            if mRNA_list == []:
                # This is not a protein coding gene
                continue
            else:
                for mRNA in mRNA_list:
                    # logger.debug(mRNA.ID)
                    mRNA_len = 0
                    CDS_list = mRNA.get_all_child_feats_obj('CDS')
                    if not CDS_list:
                        logger.error("ERROR: mRNA {} do not have CDS type".format(mRNA))
                        non_gene_flag = True
                        continue
                    else:
                        for CDS in CDS_list:
                            mRNA_len += CDS.len
                    mRNA.abs_len = mRNA_len
                    # logger.debug(mRNA.abs_len)
                    if longest < mRNA.abs_len:
                        longest = mRNA.abs_len
                        self.GFF_dict[k].longest = mRNA.ID
            if non_gene_flag:
                logger.error("ERROR: gene {} do not have CDS type".format(k))
                continue
            else:
                longest_table += "{}\t{}".format(k, self.GFF_dict[k].longest) + "\n"
                longest_gff += self.GFF_dict[k].content + "\n"
                longest_gff += self.GFF_dict[self.GFF_dict[k].longest].get_all_child_feats()
        return [longest_table, longest_gff]

    def append(self, chr, start, end, strand, name, type='CDS', score='.', source='.'):
        """
        Insert new gff record
        :param chr:
        :param start:
        :param end:
        :param strand:
        :param name:
        :return:
        """
        feat = Feat()
        feat.new(chr=chr, start=start, end=end, strand=strand, name=name, type=type, score=score, source=source)
        self.top_level_list.append(feat.attr_dict['ID'])
        self.GFF_dict[feat.ID] = feat


def longest_mRNA(gff=None):
    r"""
    print mRNA's longest to gff.longest.table and gff.longest.gff
    :param gff:
    :return:
    """
    gff_obj = GFF(gff)
    (longest_table, longest_gff) = gff_obj.longest_mRNA()
    with open(gff + "longest.table", 'w') as fh:
        fh.write(longest_table)
    with open(gff + "longest.gff", 'w') as fh:
        fh.write(longest_gff)
    return 0


def extract_gff_tag(gff=None, tag=None):
    r"""
    Extract specific attribute of gene or mRNA or any Item by the attribute name
    like python %s maker.gff _AED
    :param GFF:
    :param tag:
    :return:
    """
    gff_db = GFF(gff)
    tag_dict = gff_db.get_attr(tag)
    for k in tag_dict:
        print("{}\t{}".format(k, tag_dict[k]))


def fix_gt_gff(gff=None):
    r"""
    Fix the resulting gff files from gt gtf2gff3
    :param gff:
    :return:
    """
    count = defaultdict(int)
    with open(gff) as fh:
        for line in fh:
            if(line.startswith('#')):
                print(line, end='')
                continue
            feat = Feat(line)
            if feat.type == "gene":
                feat.update_tag("ID",
                                feat.attr_dict['gene_id'].replace('gene:', ''))
            elif feat.type == "mRNA":
                feat.update_tag("Parent",
                                feat.attr_dict['gene_id'].replace('gene:', ''))
                feat.update_tag("ID",
                                feat.attr_dict['transcript_id'])
            else:
                feat.parent = feat.attr_dict['transcript_id']
                prefix = "{}:{}".format(feat.parent, feat.type)
                count[prefix] += 1
                feat.update_tag("ID",
                                "{}-{}".format(prefix, count[prefix]))
            print(feat.content)


def gff2bed(GFF=None):
    """
    Convert gff to bed with jcvi
    :param GFF:
    :return: Output GFF.bed
    """
    cmd = "python -m jcvi.formats.gff bed --type=mRNA --key=ID {0} -o {0}.bed".format(GFF)
    sh(cmd)
    return GFF + ".bed"


#Bed relavent utils
class Loci:
    """
    Loci object which could also be looked as bed object
    """

    def __init__(self, chr='.', start='-1', end='-1', name='.', score='.', strand='.'):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = score
        self.strand = strand

    def get_size(self, bed_format=True):
        result = int(self.end) - int(self.start)
        if not bed_format:
            result += 1
        return result

    def get_line(self):
        return "\t".join([self.chr, str(self.start), str(self.end), self.name, self.score, self.strand]) + "\n"


class Bed:
    """
    A BED class that support read bed files, store them into a dict
    """

    def __init__(self, bed):
        # The bed list
        self.bed_list = []
        # The dict with gene name as keys
        self.bed_dict = {}
        self.load(bed)

    def load(self, bed):
        with open(bed) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                (chr, start, end, name, score, strand) = ['.'] * 6
                mylist = line.rstrip().split('\t')
                if len(mylist) >= 3:
                    (chr, start, end) = mylist[:3]
                    for i in range(3):
                        mylist.pop(0)
                else:
                    logger.error("Wrong format of BED, is it tab delmited with at least three field?")
                    exit(1)
                if len(mylist) >= 1:
                    name = mylist.pop(0)
                if len(mylist) >= 1:
                    score = mylist.pop(0)
                if len(mylist) >= 1:
                    strand = mylist.pop(0)
                loci = Loci(chr, start, end, name, score, strand)
                self.bed_list.append(loci)
                self.bed_dict[name] = loci

    def select_name(self, name, format='loci'):
        """
        Select bed by the name, returning the loci object
        :param name:
        :param format: loci|bed, loci will return the Loci object while bed will return text formated bed lines
        :return:
        """
        return_list = []
        if '\n' in name:
            name_list = name.splitlines()
        elif type(name) == str:
            name_list = [name]
        else:
            name_list = name
        for this_name in name_list:
            this_loci = self.bed_dict[this_name]
            return_list.append(this_loci)
        if format == 'loci':
            if len(return_list) == 1:
                return_list = return_list[0]
            return return_list
        if format == 'bed':
            return_text = ''
            for r in return_list:
                return_text += r.get_line() + "\n"
            return return_text

    def sum_size(self):
        """
        Total size of bed
        :return:
        """
        sum = 0
        for i in self.bed_list:
            sum += i.get_size()
        return sum

    def stat(self, short='F'):
        """
        stat bed size
        :param short: [T/F], whether to use short format, T is use
        :return:
        """
        size_list = []
        # threshold = 50
        for i in self.bed_list:
            size_list.append(i.get_size())
        if short != 'T':
            print("""Max: {}
            Min: {}
            Mean: {}
            Count: {}
            Sum: {}
            """.format(max(size_list), min(size_list), mean(size_list), len(size_list), sum(size_list)))
        else:
            print(sum(size_list))
        return 0

    def write(self, output):
        """
        write bed to output
        :return:
        """
        buffer = ''
        for i in self.bed_list:
            buffer += i.get_line()
        with open(output, 'w') as fh:
            fh.write(buffer)

    def change_offset(self, seq_id, offset):
        """
        :param seq_id:
        :param offset:
        :return:
        """
        offset = int(offset)
        for i in self.bed_list:
            if i.chr == seq_id:
                i.start = max(0, i.start - offset)
                i.end = max(0, i.end - offset)


def select_bed_by_name(gene_list_file=None, gene_bed=None):
    """
    :param gene_list_file: Input gene_list, just make sure the first column is gene name
    :param gene_bed: Input gene bed files
    :return: selected bed file name, usually gene_list_file.bed
    """
    gene_bed_obj = Bed(gene_bed)
    gene_list = []
    with open(gene_list_file) as fh:
        for line in fh:
            gene_list.append(line.split()[0])
    select_bed_text = gene_bed_obj.select_name(gene_list, format='bed')
    select_bed_file = gene_list_file + '.bed'
    with open(select_bed_file, 'w') as fh:
        fh.write(select_bed_text)
    return select_bed_file


def sum_bed(bed_file=None):
    """
    Input:
        chrUN   0       100000  99
        chrUN   0       100000  99
        chrUN   0       100000  99
    Output:
        chrUN   0       100000  297
    :param bed_file:
    :return:
    """
    dict_abc = defaultdict(int)
    if os.path.exists(bed_file):
        with open(bed_file) as fh:
            buffer = fh.read()
    elif bed_file == '-':
        buffer = sys.stdin.read()
    else:
        logger.error("INPUT error: {}\nEither '-' or a real file input is supported".format(bed_file))
        return 1
    for line in buffer.splitlines():
        mylist = line.strip().split()
        this_key = "\t".join(mylist[0:3])
        this_value = int(mylist[3])
        dict_abc[this_key] += this_value
    for k in dict_abc:
        print("{}\t{}".format(k, dict_abc[k]))


if __name__ == "__main__":
    emain()

