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
    score =r"""-o:-11
-e:-1
AA:4
AR:-1
AN:-2
AD:-2
AC:0
AQ:-1
AE:-1
AG:0
AH:-2
AI:-1
AL:-1
AK:-1
AM:-1
AF:-2
AP:-1
AS:1
AT:0
AW:-3
AY:-2
AV:0
RR:5
RN:0
RD:-2
RC:-3
RQ:1
RE:0
RG:-2
RH:0
RI:-3
RL:-2
RK:2
RM:-1
RF:-3
RP:-2
RS:-1
RT:-1
RW:-3
RY:-2
RV:-3
NN:6
ND:1
NC:-3
NQ:0
NE:0
NG:0
NH:1
NI:-3
NL:-3
NK:0
NM:-2
NF:-3
NP:-2
NS:1
NT:0
NW:-4
NY:-2
NV:-3
DD:6
DC:-3
DQ:0
DE:2
DG:-1
DH:-1
DI:-3
DL:-4
DK:-1
DM:-3
DF:-3
DP:-1
DS:0
DT:-1
DW:-4
DY:-3
DV:-3
CC:9
CQ:-3
CE:-4
CG:-3
CH:-3
CI:-1
CL:-1
CK:-3
CM:-1
CF:-2
CP:-3
CS:-1
CT:-1
CW:-2
CY:-2
CV:-1
QQ:5
QE:2
QG:-2
QH:0
QI:-3
QL:-2
QK:1
QM:0
QF:-3
QP:-1
QS:0
QT:-1
QW:-2
QY:-1
QV:-2
EE:5
EG:-2
EH:0
EI:-3
EL:-3
EK:1
EM:-2
EF:-3
EP:-1
ES:0
ET:-1
EW:-3
EY:-2
EV:-2
GG:6
GH:-2
GI:-4
GL:-4
GK:-2
GM:-3
GF:-3
GP:-2
GS:0
GT:-2
GW:-2
GY:-3
GV:-3
HH:8
HI:-3
HL:-3
HK:-1
HM:-2
HF:-1
HP:-2
HS:-1
HT:-2
HW:-2
HY:2
HV:-3
II:4
IL:2
IK:-3
IM:1
IF:0
IP:-3
IS:-2
IT:-1
IW:-3
IY:-1
IV:3
LL:4
LK:-2
LM:2
LF:0
LP:-3
LS:-2
LT:-1
LW:-2
LY:-1
LV:1
KK:5
KM:-1
KF:-3
KP:-1
KS:0
KT:-1
KW:-3
KY:-2
KV:-2
MM:5
MF:0
MP:-2
MS:-1
MT:-1
MW:-1
MY:-1
MV:1
FF:6
FP:-4
FS:-2
FT:-2
FW:1
FY:3
FV:-1
PP:7
PS:-1
PT:-1
PW:-4
PY:-3
PV:-2
SS:4
ST:1
SW:-3
SY:-2
SV:-2
TT:5
TW:-2
TY:-2
TV:0
WW:11
WY:2
WV:-3
YY:7
YV:-1
VV:4
AB:-2
AJ:-1
AZ:-1
RB:-1
RJ:-2
RZ:0
NB:4
NJ:-3
NZ:0
DB:4
DJ:-3
DZ:1
CB:-3
CJ:-1
CZ:-3
QB:0
QJ:-2
QZ:4
EB:1
EJ:-3
EZ:4
GB:-1
GJ:-4
GZ:-2
HB:0
HJ:-3
HZ:0
IB:-3
IJ:3
IZ:-3
LB:-4
LJ:3
LZ:-3
KB:0
KJ:-3
KZ:1
MB:-3
MJ:2
MZ:-1
FB:-3
FJ:0
FZ:-3
PB:-2
PJ:-3
PZ:-1
SB:0
SJ:-2
SZ:0
TB:-1
TJ:-1
TZ:-1
WB:-4
WJ:-2
WZ:-2
YB:-3
YJ:-1
YZ:-2
VB:-3
VJ:2
VZ:-2
BB:4
BJ:-3
BZ:0
JJ:3
JZ:-3
ZZ:4
Xx:-1
*x:-4
**:1
"""
    cmd0 = "echo {} > alignscore.txt; ln -s `which formatdb`; ln -s `which blastall`;".format(score)
    cmd2 = "genblastG -p genblastg -q {} -t {} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r {} -c 0.5 -s 0 -i 15 \
-x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {}.genblast"
    cmd1 = "genblastG -p genblastg -q {} -t {} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r {} -c 0.5 -s 0 -i 15 \
-x 20 -n 20 -v 2 -h 1 -j 0 -norepair -gff -cdna -pro -o {}.genblast"
    cmd = cmd0 + cmd2.format(qry, ref, rank, PREFIX)
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
