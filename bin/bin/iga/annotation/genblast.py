"""
genblast relevant utils
"""

import re

from iga.apps.base import emain, sh
from Bio import SeqIO

import logging


class Feat:
    def __init__(self, gene_id, content, rank=-1, score=0, length=0, exon_num=0):
        self.gene_id = gene_id
        self.content = content
        self.exon_num = exon_num
        self.length = length
        rank_search = re.search(r'(.*)-R(\d+)', gene_id)
        if rank_search is not None:
            self.core_gene_id = rank_search[1]
            self.rank = int(rank_search[2])
        if score == '0.0' or score == '0' or '-' in score or score == '.':
            self.score = 1
        else:
            self.score = float(score)

    def append(self, line):
        self.content += line
        self.exon_num += 1


class Best:
    def __init__(self, score, rank, gene_id):
        self.score = score
        self.rank = rank
        self.gene_id = gene_id


class GFF:
    def __init__(self, filename, transcript_key="transcript"):
        self.GFF_dict = {}
        self.best = {}
        with open(filename) as fh:
            transcript_id = ""
            score = 0
            count = 0
            for line in fh:
                count += 1
                if line.startswith('#') or line.rstrip() == "":
                    continue
                mylist = line.rstrip().split()
                if len(mylist) < 6 or len(mylist) < 3:
                    logging.warning("Error: Array element number is {} on \
                            line {}".format(len(mylist), count))
                    exit()
                if mylist[2] == transcript_key:
                    transcript_id = self.get_ID_from_last_feat(mylist[-1])
                    score = mylist[5]
                    gene_len = int(mylist[4]) - int(mylist[3])
                    self.GFF_dict[transcript_id] = \
                        Feat(gene_id=transcript_id, score=mylist[5], content=line, length=gene_len)

                    core_gene_id = self.GFF_dict[transcript_id].core_gene_id
                    this_score = self.GFF_dict[transcript_id].score
                    this_rank = self.GFF_dict[transcript_id].rank

                    if (core_gene_id not in self.best or
                            this_score > self.best[core_gene_id].score):
                        self.best[core_gene_id] = Best(this_score, this_rank, transcript_id)
                else:
                    #                    print(transcript_id)
                    #                    print(line)
                    self.GFF_dict[transcript_id].append(line)

    @staticmethod
    def get_ID_from_last_feat(lastfeat):
        lastfeat = re.sub(r";.*", "", lastfeat.replace('ID=', ''))
        return lastfeat

    def filter(self, div_threshold=0.9, rank_keep=2, exon_num=2, len_cutoff=50000):
        my_new_dict = {}
        for gene_id in self.GFF_dict:
            core_gene_id = self.GFF_dict[gene_id].core_gene_id
            try:
                div_result = self.GFF_dict[gene_id].score / self.best[core_gene_id].score
            except ZeroDivisionError as e:
                print("\t".join([gene_id, str(self.GFF_dict[gene_id].score), str(self.best[core_gene_id].score)]))

            if (self.GFF_dict[gene_id].rank > rank_keep or
                    self.GFF_dict[gene_id].score / self.best[core_gene_id].score <
                    div_threshold or
                    self.GFF_dict[gene_id].length > len_cutoff and
                    self.GFF_dict[gene_id].exon_num == exon_num or
                    self.GFF_dict[gene_id].length > len_cutoff and
                    self.GFF_dict[gene_id].rank > rank_keep):
                pass
            else:
                my_new_dict[gene_id] = self.GFF_dict[gene_id]

        self.GFF_dict = my_new_dict.copy()

    def print_out(self):
        for k in self.GFF_dict:
            print(self.GFF_dict[k].content, end='')


#def run_genblast():


def filter_genblast(genblast_gff=None):
    """
    Filter genblast result.
    :param genblast_gff:
    :return:
    """
    #   000028F|arrow_np1212    genBlastG       transcript      2151627 2152061 23.3019 -       .       ID=sp_A7M942_PSAC_CUSGR
    #   000028F|arrow_np1212    genBlastG       coding_exon     2151627 2152061 .       -       .       ID=sp_A7M942_PSAC_CUSGR
    #   000005F|arrow_np1212    genBlastG       transcript      3483885 3484160 23.2295 -       .       ID=sp_A7M942_PSAC_CUSGR

    gff_read = GFF(genblast_gff)
    gff_read.filter()
    gff_read.print_out()


def filter_early_stop(fasta=None):
    """
    Filter fasta contaning early stops
    :param fasta: pep fasta
    :return:
    """

    raw_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    genome_dict = {}
    # rename genome_dict
    for old_key in raw_dict.keys():
        if re.search(r'[\*|U].', raw_dict[old_key].seq.__str__()):
            #logging.debug('{} with seq {} contains gap'.format(old_key, raw_dict[old_key].seq))
            logging.debug(re.search(r'\*|U[ATCG]', raw_dict[old_key].seq.__str__()))
            continue
        else:
            genome_dict[old_key] = raw_dict[old_key]

    with open(fasta + "noStop.fa", "w") as output_handle:
        for k in genome_dict.keys():
            SeqIO.write(genome_dict[k], output_handle, "fasta")
            print(k)


def run(qry=None, ref=None, PREFIX=None, rank=3):
    """
    run genblast
    :param qry:
    :param ref:
    :return:
    """
    cmd2 = "genblast -p genblastg -q {} -t {} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r {} -c 0.5 -s 0 -i 15 \
-x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {}.genblast"
    cmd1 = "genblast -p genblastg -q {} -t {} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r {} -c 0.5 -s 0 -i 15 \
-x 20 -n 20 -v 2 -h 1 -j 0 -norepair -gff -cdna -pro -o {}.genblast"
    cmd = cmd2.format(qry, ref, rank, PREFIX)
    job = sh(cmd)
    if type(job) == int:
        logging.error("Genblast (h2) failed with return code {}".format(job))
        sh("rm {}.genblast_1.1c_2.3_s2_tdshift2_tddis0_tcls0.0_m2_score_i0_d16_0*".format(PREFIX))
        cmd = cmd1.format(qry, ref, rank, PREFIX)
        job = sh(cmd)
        if type(job) == int:
            logging.error("Genblast (h1) failed with return code {}".format(job))
            exit(1)


if __name__ == "__main__":
    emain()
