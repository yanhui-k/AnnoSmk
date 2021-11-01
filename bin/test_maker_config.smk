# -*- coding = utf-8 -*-
# @Time : 2021/10/29 10:22
# @Author : 严慧
# @File : maker_config.py
# @Software : PyCharm

configfile: "config/config.yaml"

import os
from configparser import ConfigParser

PREFIX = config["PREFIX"]
REF = config["REF"]
FLNCESTGFF = os.path.abspath(config["FLNCESTGFF"])
ESTGFF = config["ESTGFF"]
CDNAFASTA = os.path.abspath(config["CDNAFASTA"])
PEPGFF = os.path.abspath(config["PEPGFF"])
REPEATGFF = os.path.abspath(config["REPEATGFF"])


def prepare_opts(estgff=None, pepgff=None, rmgff=None, round=None,
                 snap_hmm="", augustus_species="", output_file=None):
    # find the abspath of estgff,pepgff,rmgff
    # if round=1,snap_hmm="",if round=2,R1/R1.hmm,if round=3,R2/R2.hmm
    # if round=1,augustus_species="",if round=2,R1,if round=3,R2
    # add the information to the opts.ctl
    if snap_hmm != "" and augustus_species != "" and round != "1":
        snap_hmm_dir = os.path.abspath(snap_hmm)
        augustus_species = augustus_species
        est2genome = "0"
        protein2genome = "0"
        alt_splice = "1"
        if round == "2":
            trna = "1"
        else:
            trna = "0"
    elif round == "1":
        snap_hmm_dir = ""
        augustus_species = ""
        est2genome = "1"
        protein2genome = "1"
        alt_splice = "0"
        trna = "0"
    else:
        exit(1)
    config = ConfigParser()
    config.read("maker_opts.ctl")
    estgff_dir = os.path.abspath(estgff)
    pepgff_dir = os.path.abspath(pepgff)
    rmgff_dir = os.path.abspath(rmgff)
    config.set("maker_opts", "est_gff", estgff_dir)
    config.set("maker_opts", "protein_gff", pepgff_dir)
    config.set("maker_opts", "rm_gff", rmgff_dir)
    config.set("maker_opts", "snaphmm", snap_hmm_dir)
    config.set("maker_opts", "augustus_species", augustus_species)
    config.set("maker_opts", "est2genome", est2genome)
    config.set("maker_opts", "protein2genome", protein2genome)
    config.set("maker_opts", "alt_splice", alt_splice)
    config.set("maker_opts", "trna", trna)
    config.set("maker_opts", "model_org", "")

    output_dir = os.path.abspath(output_file)
    with open("opts.yaml", "w", encoding="utf-8") as file:
        config.write(file)
    lines = open("opts.yaml").readlines()
    file = open(output_dir, "w")
    for s in lines:
        s = s.replace(" =", "=")
        s = s.replace("aed_threshold", "AED_threshold")
        file.write(s.replace("tmp=", "TMP="))
    file.close()


rule all:
    input:
        expand("result/{PREFIX}/R{round}/maker_opts{round}_{pren}.ctl", PREFIX=PREFIX, round=1, pren=0),
        expand("result/{PREFIX}/R{round}/maker_opts{round}_{pren}.ctl", PREFIX=PREFIX, round=2, pren=1),
        expand("result/{PREFIX}/R{round}/maker_opts{round}_{pren}.ctl", PREFIX=PREFIX, round=3, pren=2)

rule pre_pre_estgff:
    input:
        expand("{FLNCESTGFF}", FLNCESTGFF=FLNCESTGFF)
    output:
        "total_est.gff"
    shell:
        "cp {input} {output}"

rule pre_pre_pepgff:
    input:
        expand("{PEPGFF}", PEPGFF=PEPGFF)
    output:
        "total_pep.gff"
    shell:
        "cp {input} {output}"

rule pre_pre_rmgff:
    input:
        expand("{REPEATGFF}", REPEATGFF=REPEATGFF)
    output:
        "rm.gff"
    shell:
        "cp {input} {output}"

rule pre_opts:
    input:
        snap_hmm = "result/{PREFIX}/R{pren}/{PREFIX}.genome.contig.fa.masked.fa_R{pren}.hmm",
        estgff = "total_est.gff",
        pepgff = "total_pep.gff",
        rmgff = "rm.gff"
    params:
        round = "{round}",
        augustus_species = "{PREFIX}.genome.contig.fa.masked.fa_R{pren}_direct"
    output:
        opts_file = "result/{PREFIX}/R{round}/maker_opts{round}_{pren}.ctl"
    run:
        prepare_opts(estgff="data_test/genblast.gff", pepgff=input.pepgff,
                    rmgff="data_test/genblast.gff", round=params.round,
                    snap_hmm=input.snap_hmm,
                    augustus_species=params.augustus_species,
                    output_file=output.opts_file)
