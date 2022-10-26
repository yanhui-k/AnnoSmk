import pathlib
import random

import glob
import os
from configparser import ConfigParser

def replace_list(list,a,b):
    new_list = []
    while list:
        new_item = list.pop().replace(a,b,1)
        new_list.append(new_item)
    return new_list

folder_flnc = os.path.join(PREFIX,"*_subreads.fastq.gz")
samples_flnc=glob.glob(folder_flnc)
prefix_flnc=replace_list(samples_flnc,"_subreads.fastq.gz","")
prefix_flnc=replace_list(prefix_flnc,PREFIX,"")
prefix_flnc=replace_list(prefix_flnc,"/","")

if samples_flnc == []:
    prefix_est=["rnaseq"]
    augustus_cdna="annotation_smk/{PREFIX}/evidence/rnaseq.fasta"
elif samples_flnc != []:
    prefix_est=["rnaseq","flnc"]
    augustus_cdna="annotation_smk/{PREFIX}/evidence/flnc.fasta"

localrules:
    all,
    pre_exe,
    pre_opts,
    merge_log,
    fathom_to_genbank,
    maker_gff_only

rule make_fasta:
    output:
        fasta = expand("{REF}",REF=REF)
    shell:
        '''
        touch {output.fasta}
        '''

checkpoint split_fasta:
    input:
        fasta = rules.make_fasta.output.fasta
    output:
        lane_dir = directory("annotation_smk/{PREFIX}/sample/")
    params:
        size = 4000000
    shell:
        "split_fasta.py {input} {output} {params.size}"

rule cp:
    input:
        fasta = "annotation_smk/{PREFIX}/sample/{lane_number}.fa"
    output:
        fa = "annotation_smk/{PREFIX}/R{round}/{lane_number}.fa"
    shell:
        "cp {input} {output}"

rule pre_0_opts:
    output:
        "annotation_smk/{PREFIX}/maker_opts_0.ctl"
    shell:
        '''
        maker_opts.sh > {output}
        '''

def prepare_opts(input_opt=None, estgff=None, pepgff=None, rmgff=None, round=None,
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
    config.read(input_opt)
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

#Modifications are required based on the evidence provided

rule pre_pre_hmm:
    output:
        hmm = expand("annotation_smk/{PREFIX}/pre_R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",round=0,PREFIX=PREFIX),
        aug_dir = expand("annotation_smk/{PREFIX}/pre_R{round}/autoAug/autoAugPred_hints/shells",round=0,PREFIX=PREFIX)
    shell:
        '''
        touch {output.hmm}
        touch {output.aug_dir}
        '''

rule pre_pre_pepgff:
    input:
        expand("annotation_smk/{PREFIX}/evidence/genblast.gff",PREFIX=PREFIX)
    output:
        "total_pep.gff"
    shell:
        "sort -n {input} | uniq > {output}"

rule pre_pre_rmgff:
    input:
        expand("annotation_smk/{PREFIX}/evidence/repeat.gff",PREFIX=PREFIX)
    output:
        "repeat.gff"
    shell:
        "cp {input} {output}"

rule pre_pre_estgff:
    input:
        expand("annotation_smk/{PREFIX}/evidence/total_est.gff",PREFIX=PREFIX)
    output:
        "total_est.gff"
    shell:
        "sort -n {input} | uniq > {output}"

hmm_dict = {"1":f"annotation_smk/{PREFIX}/pre_R0/{PREFIX}.genome.contig.fa.masked.fa_R0.hmm",
            "2":f"annotation_smk/{PREFIX}/R1/{PREFIX}.genome.contig.fa.masked.fa_R1.hmm",
            "3":f"annotation_smk/{PREFIX}/R2/{PREFIX}.genome.contig.fa.masked.fa_R2.hmm"}

def get_hmm(wildcards):
    round = wildcards.round
    hmm_file = hmm_dict[round]
    return hmm_file

augustus_species_dict = {"1":"","2":f"{PREFIX}.genome.contig.fa.masked.fa_R1_direct",
                         "3":f"{PREFIX}.genome.contig.fa.masked.fa_R2_direct"}

def get_augustus_species(wildcards):
    round = wildcards.round
    augustus_species = augustus_species_dict[round]
    return augustus_species

augustus_dict = {"1":f"annotation_smk/{PREFIX}/pre_R0/autoAug/autoAugPred_hints/shells",
                 "2":f"annotation_smk/{PREFIX}/R1/autoAug/autoAugPred_hints/shells",
                 "3":f"annotation_smk/{PREFIX}/R2/autoAug/autoAugPred_hints/shells"}

def get_augustus_dir(wildcards):
    round = wildcards.round
    augustus_dir = augustus_dict[round]
    return augustus_dir

rule pre_opts:
    input:
        input_opt = "annotation_smk/{PREFIX}/maker_opts_0.ctl",
        snap_hmm = get_hmm,
        estgff = "total_est.gff",
        pepgff = "total_pep.gff",
        rmgff = "repeat.gff",
        augustus_dir = get_augustus_dir
    params:
        round = "{round}",
        augustus_species = get_augustus_species
    output:
        opts_file = "annotation_smk/{PREFIX}/R{round}/maker_opts{round}.ctl"
    run:
        prepare_opts(input_opt=input.input_opt, estgff=input.estgff, pepgff=input.pepgff, 
                    rmgff=input.rmgff, round=params.round,
                    snap_hmm=input.snap_hmm,
                    augustus_species=params.augustus_species,
                    output_file=output.opts_file)

# def get_opts(wildcards):
#     opts_file = checkpoints.pre_exe.get(**wildcards).output[0]
#     round = glob_wildcards(f"annotation_smk/{wildcards.PREFIX}/R{{round}}/maker_exe.ctl").round
#     # for num in round:
#     #     pren = int(num) - 1
#     opts_file = expand(rules.pre_opts.output,**wildcards,pren=lambda w: w-1)
#     return opts_file

checkpoint pre_exe:
    output:
        "annotation_smk/{PREFIX}/R{round}/maker_exe.ctl",
        "annotation_smk/{PREFIX}/R{round}/maker_bopts.ctl"
    shell:
        '''
        maker -CTL 2>/dev/null
        cp maker_exe.ctl annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/
        cp maker_bopts.ctl annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/
        '''

rule run_maker:
    input:
        g="annotation_smk/{PREFIX}/R{round}/{lane_number}.fa",
        opts="annotation_smk/{PREFIX}/R{round}/maker_opts{round}.ctl",
        bopts="annotation_smk/{PREFIX}/R{round}/maker_bopts.ctl",
        exe="annotation_smk/{PREFIX}/R{round}/maker_exe.ctl"
    output:
        log="annotation_smk/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_master_datastore_index.log"
    threads:
        THREADS
    shell:
        '''
        maker -c {threads} -genome {input.g} {input.opts} {input.bopts} {input.exe} -fix_nucleotides
        wait
        cp -rf {wildcards.lane_number}.maker.output annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/
        rm -rf {wildcards.lane_number}.maker.output
        '''

rule alt_log:
    input:
        rules.run_maker.output.log
    output:
        "annotation_smk/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_total_master_datastore_index.log"
    shell:
        '''
        cat {input} |sed "s/\t/\t{wildcards.lane_number}.maker.output\//" \
        > {output}
        '''

def get_log(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"annotation_smk/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    log = expand(rules.alt_log.output, **wildcards, lane_number=lane_numbers)
    return log

rule merge_log:
    input:
        get_log
    output:
        "annotation_smk/{PREFIX}/R{round}/total_master_datastore_index.log"
    shell:
        '''
        cat {input} > {output}
        '''

rule gff3_merge:
    input:
        rules.merge_log.output
    output:
        all_gff="annotation_smk/{PREFIX}/R{round}/genome.all.gff",
        noseq_gff="annotation_smk/{PREFIX}/R{round}/genome.all.noseq.gff",
        all_fasta="annotation_smk/{PREFIX}/R{round}/total.all.maker.proteins.fasta"
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
        "annotation_smk/{PREFIX}/R{round}/genome.maker.gff"
    shell:
        '''
        awk '$2=="maker"' {input} > {output}
        '''

def get_fa(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"annotation_smk/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    fa = expand(rules.cp.output.fa, **wildcards, lane_number=lane_numbers)
    return fa

rule get_ref_fa:
    input:
        get_fa
    output:
        "annotation_smk/{PREFIX}/R{round}/ref.fa"
    shell:
        "cat {input} > {output}"

rule maker2zff:
    input:
        "annotation_smk/{PREFIX}/R{round}/genome.all.gff"
    output:
        ann="annotation_smk/{PREFIX}/R{round}/genome.ann",
        dna="annotation_smk/{PREFIX}/R{round}/genome.dna"
    shell:
        '''
        maker2zff -x 0.25 -l 50 {input}
        mv genome.ann annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv genome.dna annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom1:
    input:
        ann="annotation_smk/{PREFIX}/R{round}/genome.ann",
        dna="annotation_smk/{PREFIX}/R{round}/genome.dna"
    output:
        "annotation_smk/{PREFIX}/R{round}/gene-stats.log"
    shell:
        "fathom -gene-stats {input.ann} {input.dna} >{output} 2>&1"

rule fathom2:
    input:
        ann="annotation_smk/{PREFIX}/R{round}/genome.ann",
        dna="annotation_smk/{PREFIX}/R{round}/genome.dna"
    output:
        "annotation_smk/{PREFIX}/R{round}/validate.log"
    shell:
        "fathom -validate {input.ann} {input.dna} >{output} 2>&1"

rule fathom3:
    input:
        ann="annotation_smk/{PREFIX}/R{round}/genome.ann",
        dna="annotation_smk/{PREFIX}/R{round}/genome.dna"
    output:
        uann="annotation_smk/{PREFIX}/R{round}/uni.ann",
        udna="annotation_smk/{PREFIX}/R{round}/uni.dna"
    shell:
        '''
        fathom -categorize 1000 {input.ann} {input.dna}
        mv uni.ann annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv uni.dna annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom4:
    input:
        uann="annotation_smk/{PREFIX}/R{round}/uni.ann",
        udna="annotation_smk/{PREFIX}/R{round}/uni.dna"
    output:
        exann="annotation_smk/{PREFIX}/R{round}/export.ann",
        exdna="annotation_smk/{PREFIX}/R{round}/export.dna"
    shell:
        '''
        fathom -export 1000 -plus {input.uann} {input.udna}
        mv export.ann annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv export.dna annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule forge_hmm_assembler:
    input:
        exann="annotation_smk/{PREFIX}/R{round}/export.ann",
        exdna="annotation_smk/{PREFIX}/R{round}/export.dna"
    output:
        "annotation_smk/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm"
    params:
        fo="{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm"
    shell:
        '''
        cd annotation_smk/{wildcards.PREFIX}/R{wildcards.round}
        forge export.ann export.dna >forge.log 2>&1
        hmm-assembler.pl snap_trained . > {params.fo}
        cd -
        '''

rule fathom_to_genbank:
    input:
        uann="annotation_smk/{PREFIX}/R{round}/uni.ann",
        udna="annotation_smk/{PREFIX}/R{round}/uni.dna"
    output:
        "annotation_smk/{PREFIX}/R{round}/augustus.gb"
    params:
        path=expand("{PATH}",PATH=PATH)
    shell:
        '''
        perl {params.path}/bin/fathom_to_genbank.pl --annotation_file {input.uann} --dna_file {input.udna}  --genbank_file {output} --number 500
        '''
#fathom_to_genbank.pl文件需要修改perl的路径

rule perl_cat:
    input:
        "annotation_smk/{PREFIX}/R{round}/augustus.gb"
    output:
        "annotation_smk/{PREFIX}/R{round}/genbank_gene_list.txt"
    shell:
        "cat.py {input} {output}"

rule get_subset_of_fastas:
    input:
        txt="annotation_smk/{PREFIX}/R{round}/genbank_gene_list.txt",
        udna="annotation_smk/{PREFIX}/R{round}/uni.dna"
    output:
        "annotation_smk/{PREFIX}/R{round}/genbank_gene_seqs.fasta"
    params:
        path=expand("{PATH}",PATH=PATH)
    shell:
        '''
        perl {params.path}/bin/get_subset_of_fastas.pl -l {input.txt} -f {input.udna} -o {output}
        '''

# rule randomSplit:
#     input:
#         "annotation_smk/{PREFIX}/R{round}/augustus.gb"
#     output:
#         "annotation_smk/{PREFIX}/R{round}/augustus.gb.test"
#     shell:
#         '''
#         randomSplit.pl {input} 250
#         '''


rule autoAugA:
    input:
        fasta="annotation_smk/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gb="annotation_smk/{PREFIX}/R{round}/augustus.gb",
        rna=augustus_cdna
    output:
        "annotation_smk/{PREFIX}/R{round}/autoAug/autoAugPred_abinitio/shells/aug1",
        "annotation_smk/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    message:
        "If this step reports an error, you can delete autoAug and ~/.conda/envs/annotation/config/species/{wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct, and try again"
    threads: THREADS
    shell:
        '''
        autoAug.pl --cpus={threads} --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --trainingset={input.gb} --cdna={input.rna} --noutr
        cd autoAug/autoAugPred_abinitio/shells
        ./aug1
        cd ../../../
        cp -r autoAug annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

#如果有autoAug文件夹或者在/nfs/yanhui/.conda/envs/annotation/config/species下
#有species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct的文件夹都会报错

rule autoAugB:
    input:
        fasta="annotation_smk/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gff="annotation_smk/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    output:
        directory("annotation_smk/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells")
    threads: THREADS
    shell:
        '''
        autoAug.pl --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --useexisting --hints={input.gff} \
        -v -v -v  --index=1
        cd autoAug/autoAugPred_hints/shells/
        augustus --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct --UTR=off \
        --hintsfile=../../hints/hints.E.gff \
        --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg --exonnames=on  \
        --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
        ../../seq/split/genome_clean.split.1.fa > aug1.out
        cd ../../../
        cp -r autoAug/autoAugPred_hints annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/autoAug/.
        rm -r autoAug
        '''

rule busco:
    input:
        rules.gff3_merge.output.all_fasta
    output:
        directory("annotation_smk/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta")
    params:
        dir_busco="total.all.maker.proteins.fasta.busco.embryophyta"
    container:"docker://ezlabgva/busco:v5.4.3_cv1"
    shell:
        """
        cd annotation_smk/{wildcards.PREFIX}/R{wildcards.round}/
        busco -f -c 64 -m prot -i ../../../{input} -o {params.dir_busco} -l embryophyta_odb10
        """
rule AED:
    input:
        rules.gff3_merge.output.all_gff
    output:
        "annotation_smk/{PREFIX}/R{round}/AED.csv"
    shell:
        "AED_cdf_generator.pl -b 0.025 {input} > {output}"


rule maker_gff_only:
    input:
        expand(rules.gff3_merge.output.all_gff,PREFIX=PREFIX, round=3)
    output:
        gff=expand("annotation_smk/{PREFIX}/maker_gff_only.gff",PREFIX=PREFIX),
        fasta=expand("annotation_smk/{PREFIX}/maker_fasta_only.fa",PREFIX=PREFIX)
    shell:
        """
        awk -v line=$(awk '/##FASTA/{{print NR}}' {input}) '{{if(NR<line){{print $0}}}}' {input} > {output.gff}
        awk -v line=$(awk '/##FASTA/{{print NR}}' {input}) '{{if(NR>line){{print $0}}}}' {input} > {output.fasta}
        """


rule format_maker_gff:
    input:
        expand(rules.maker_gff_only.output.gff,PREFIX=PREFIX, round=3)
    output:
        expand("annotation_smk/{PREFIX}.gff",PREFIX=PREFIX)
    shell:
        "maker_gff.py {input} {output}"




















