configfile:"config/config.yaml"

import os

PREFIX=config["PREFIX"]
REF=config["REF"]
PEP=config["PEP"]

include:"rules/rule.evidence.smk",
include:"rules/rule.makerfor2.smk"

wildcard_constraints:
    PREFIX=PREFIX

rule all:
    input:
        expand("result/{PREFIX}/evidence/rnaseq.fasta",PREFIX=PREFIX),
        expand("result/{PREFIX}/evidence/total_est.gff",PREFIX=PREFIX),
        expand("result/{PREFIX}/evidence/repeat.gff",PREFIX=PREFIX),
        expand("result/{PREFIX}/evidence/genblast.gff",PREFIX=PREFIX),
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=3),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=3),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=3)
