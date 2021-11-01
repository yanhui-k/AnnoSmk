import pathlib
import random

import glob
import os
from configparser import ConfigParser

rule make_fasta:
    output:
        fasta = expand("{REF}",REF=REF)
    shell:
        "touch {output.fasta}"

checkpoint split_fasta:
    input:
        fasta = rules.make_fasta.output.fasta
    output:
        lane_dir = directory("result/{PREFIX}/sample/")
    script:
        "../bin/split_fasta1.py"

rule mv:
    input:
        fasta = "result/{PREFIX}/sample/{lane_number}.fa"
    output:
        fa = "result/{PREFIX}/R{round}/{lane_number}.fa"
    shell:
        "cp {input} {output}"

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
    config.read("config/maker_opts.ctl")
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

rule pre_pre_opts:
    output:
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",round=0,PREFIX=PREFIX)

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

def get_opts(wildcards):
    opts_file = checkpoints.pre_exe.get(**wildcards).output[0]
    round = glob_wildcards(f"result/{wildcards.PREFIX}/R{{round}}/maker_exe.ctl").round
    for num in round:
        pren = int(num) -1
    opts_file = expand(rules.pre_opts.output,**wildcards,pren=pren)
    return opts_file

checkpoint pre_exe:
    output:
        "result/{PREFIX}/R{round}/maker_exe.ctl",
        "result/{PREFIX}/R{round}/maker_bopts.ctl"
    shell:
        '''
        maker -CTL
        mv maker_exe.ctl result/{wildcards.PREFIX}/R{wildcards.round}/
        mv maker_bopts.ctl result/{wildcards.PREFIX}/R{wildcards.round}/
        '''

# rule pre_bopts/exe:
# 现在用上面的方法
#     input:
#         "config/bopts.txt"
#     output:
#         "result/{PREFIX}/R{round}/maker_bopts.ctl"
#     shell:
#         '''
#         cat {input} > {output}
#         '''

rule run_maker:
    input:
        g="result/{PREFIX}/R{round}/{lane_number}.fa",
        opts=get_opts,
        bopts="result/{PREFIX}/R{round}/maker_bopts.ctl",
        exe="result/{PREFIX}/R{round}/maker_exe.ctl"
    output:
        log="result/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_master_datastore_index.log"
    shell:
        '''
        maker -genome {input.g} {input.opts} {input.bopts} {input.exe}
        wait
        cp -rf {wildcards.lane_number}.maker.output result/{wildcards.PREFIX}/R{wildcards.round}/
        rm -rf {wildcards.lane_number}.maker.output
        '''

rule alt_log:
    input:
        rules.run_maker.output
    output:
        "result/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_total_master_datastore_index.log"
    shell:
        '''
        cat {input} |sed "s/\t/\t{wildcards.lane_number}.maker.output\//" \
        > {output}
        '''

def get_log(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"result/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    log = expand(rules.alt_log.output, **wildcards, lane_number=lane_numbers)
    return log

rule merge_log:
    input:
        get_log
    output:
        "result/{PREFIX}/R{round}/total_master_datastore_index.log"
    shell:
        '''
        cat {input} > {output}
        '''

rule gff3_merge:
    input:
        rules.merge_log.output
    output:
        all_gff="result/{PREFIX}/R{round}/genome.all.gff",
        noseq_gff="result/{PREFIX}/R{round}/genome.all.noseq.gff",
        all_fasta="result/{PREFIX}/R{round}/total.all.maker.proteins.fasta"
    shell:
        '''
        fasta_merge -d {input}
        gff3_merge -o {output.all_gff} -d {input}
        gff3_merge -n -o {output.noseq_gff} -d {input}
        mv total.all.maker.proteins.fasta {output.all_fasta}
        '''

rule get_genome_maker_gff:
    input:
        rules.gff3_merge.output.noseq_gff
    output:
        "result/{PREFIX}/R{round}/genome.maker.gff"
    shell:
        '''
        awk '$2=="maker"' {input} > {output}
        '''

def get_fa(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"result/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    fa = expand(rules.mv.output, **wildcards, lane_number=lane_numbers)
    return fa

rule get_ref_fa:
    input:
        get_fa
    output:
        "result/{PREFIX}/R{round}/ref.fa"
    shell:
        "cat {input} > {output}"

rule maker2zff:
    input:
        "result/{PREFIX}/R{round}/genome.all.gff"
    output:
        ann="result/{PREFIX}/R{round}/genome.ann",
        dna="result/{PREFIX}/R{round}/genome.dna"
    shell:
        '''
        maker2zff -x 0.25 -l 50 {input}
        mv genome.ann result/{wildcards.PREFIX}/R{wildcards.round}/.
        mv genome.dna result/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom1:
    input:
        ann="result/{PREFIX}/R{round}/genome.ann",
        dna="result/{PREFIX}/R{round}/genome.dna"
    output:
        "result/{PREFIX}/R{round}/gene-stats.log"
    shell:
        "fathom -gene-stats {input.ann} {input.dna} >{output} 2>&1"

rule fathom2:
    input:
        ann="result/{PREFIX}/R{round}/genome.ann",
        dna="result/{PREFIX}/R{round}/genome.dna"
    output:
        "result/{PREFIX}/R{round}/validate.log"
    shell:
        "fathom -validate {input.ann} {input.dna} >{output} 2>&1"

rule fathom3:
    input:
        ann="result/{PREFIX}/R{round}/genome.ann",
        dna="result/{PREFIX}/R{round}/genome.dna"
    output:
        uann="result/{PREFIX}/R{round}/uni.ann",
        udna="result/{PREFIX}/R{round}/uni.dna"
    shell:
        '''
        fathom -categorize 1000 {input.ann} {input.dna}
        mv uni.ann result/{wildcards.PREFIX}/R{wildcards.round}/.
        mv uni.dna result/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom4:
    input:
        uann="result/{PREFIX}/R{round}/uni.ann",
        udna="result/{PREFIX}/R{round}/uni.dna"
    output:
        exann="result/{PREFIX}/R{round}/export.ann",
        exdna="result/{PREFIX}/R{round}/export.dna"
    shell:
        '''
        fathom -export 1000 -plus {input.uann} {input.udna}
        mv export.ann result/{wildcards.PREFIX}/R{wildcards.round}/.
        mv export.dna result/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule forge:
    input:
        exann="result/{PREFIX}/R{round}/export.ann",
        exdna="result/{PREFIX}/R{round}/export.dna"
    output:
        "result/{PREFIX}/R{round}/forge.log"
    shell:
        "forge {input.exann} {input.exdna} >{output} 2>&1"

rule hmm_assembler:
    input:
        files=rules.forge.output,
        dir="result/{PREFIX}/R{round}/"
    output:
        "result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm"
    shell:
        '''
        hmm-assembler.pl snap_trained {input.dir} > {output}
        '''

rule fathom_to_genbank:
    input:
        uann="result/{PREFIX}/R{round}/uni.ann",
        udna="result/{PREFIX}/R{round}/uni.dna"
    output:
        "result/{PREFIX}/R{round}/augustus.gb"
    shell:
        '''
        fathom_to_genbank.pl --annotation_file {input.uann} --dna_file {input.udna}  --genbank_file {output} --number 500
        '''
#fathom_to_genbank.pl文件需要修改perl的路径

rule perl_cat:
    input:
        "result/{PREFIX}/R{round}/augustus.gb"
    output:
        "result/{PREFIX}/R{round}/genbank_gene_list.txt"
    script:
        "../bin/cat.py"

rule get_subset_of_fastas:
    input:
        txt="result/{PREFIX}/R{round}/genbank_gene_list.txt",
        udna="result/{PREFIX}/R{round}/uni.dna"
    output:
        "result/{PREFIX}/R{round}/genbank_gene_seqs.fasta"
    shell:
        '''
        get_subset_of_fastas.pl -l {input.txt} -f {input.udna} -o {output}
        '''

rule randomSplit:
    input:
        "result/{PREFIX}/R{round}/augustus.gb"
    output:
        "result/{PREFIX}/R{round}/augustus.gb.test"
    shell:
        '''
        randomSplit.pl {input} 250
        '''

rule autoAugA:
    input:
        fasta="result/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gb="result/{PREFIX}/R{round}/augustus.gb",
        cdna=expand("{CDNAFASTA}",CDNAFASTA=CDNAFASTA)
    output:
        "result/{PREFIX}/R{round}/autoAug/autoAugPred_abinitio/shells/aug1",
        "result/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    shell:
        '''
        autoAug.pl --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --trainingset={input.gb} --cdna={input.cdna} --noutr
        cp -r autoAug result/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

#如果有autoAug文件夹或者在/nfs/yanhui/.conda/envs/repeat/config/species下
#有species=altra.genome.contig.fa.masked.fa_R1_direct的文件夹都会报错

rule autoAugB:
    input:
        fasta="result/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gff="result/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    output:
        directory("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells")
    shell:
        '''
        autoAug.pl --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --useexisting --hints={input.gff} \
        -v -v -v  --index=1
        mv -f autoAug/autoAugPred_hints result/{wildcards.PREFIX}/R{wildcards.round}/autoAug/.
        rm -r autoAug
        '''

rule busco:
    input:
        rules.gff3_merge.output.all_fasta
    output:
        "result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta"
    params:
        dir_busco="total.all.maker.proteins.fasta.busco.embryophyta"
    conda:
        "../env/busco.yaml"
    shell:
        """
        cd result/{wildcards.PREFIX}/R{wildcards.round}/
        busco -f -c 64 -m prot -i ../../../{input} -o {params.dir_busco} -l embryophyta_odb10
        """
























