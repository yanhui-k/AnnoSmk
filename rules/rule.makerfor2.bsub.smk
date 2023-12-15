import pathlib
import random

import glob
import os
from configparser import ConfigParser

if config["FLNC"] == None:
    prefix_est=["rnaseq"]
    augustus_cdna="rnaseq.fasta"
elif config["FLNC"] != None:
    prefix_est=["rnaseq","flnc"]
    augustus_cdna="flnc.fasta"

localrules:
    all,
    pre_exe,
    pre_opts,
    merge_log,
    fathom_to_genbank,
    maker_gff_only


checkpoint split_fasta:
    input:
        fasta="AnnoSmk/{PREFIX}/{PREFIX}.fa"
    output:
        lane_dir=directory("AnnoSmk/{PREFIX}/sample/")
    params:
        size=1000000
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_split_fasta.tsv"
    shell:
        "split_fasta.py {input} {output} {params.size}"

rule cp:
    input:
        fasta="AnnoSmk/{PREFIX}/sample/{lane_number}.fa"
    output:
        fa="AnnoSmk/{PREFIX}/R{round}/{lane_number}.fa"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_cp_R{round}_{lane_number}.tsv"
    shell:
        "cp {input} {output}"

rule pre_0_opts:
    output:
        "AnnoSmk/{PREFIX}/maker_opts_0.ctl"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_maker_opts_0.tsv"
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
            trna = "0"
        else:
            trna = "1"
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
        hmm=expand("AnnoSmk/{PREFIX}/pre_R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",round=0,PREFIX=PREFIX),
        aug_dir=expand("AnnoSmk/{PREFIX}/pre_R{round}/autoAug/autoAugPred_hints/shells",round=0,PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/pre_pre_hmm_R0.tsv"
    shell:
        '''
        touch {output.hmm}
        touch {output.aug_dir}
        '''

rule pre_pre_pepgff:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/genblast.gff",PREFIX=PREFIX)
    output:
        "total_pep.gff"
    benchmark:
        "AnnoSmk/benchmark/pre_pre_pepgff.tsv"
    shell:
        "cp {input} {output}"

rule pre_pre_rmgff:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/repeat.gff",PREFIX=PREFIX)
    output:
        "repeat.gff"
    benchmark:
        "AnnoSmk/benchmark/pre_pre_rmgff.tsv"
    shell:
        "cp {input} {output}"

rule pre_pre_estgff:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/total_est.gff",PREFIX=PREFIX)
    output:
        "total_est.gff"
    benchmark:
        "AnnoSmk/benchmark/pre_pre_estgff.tsv"
    shell:
        "cp {input} {output}"

hmm_dict = {"1":f"AnnoSmk/{PREFIX}/pre_R0/{PREFIX}.genome.contig.fa.masked.fa_R0.hmm",
            "2":f"AnnoSmk/{PREFIX}/R1/{PREFIX}.genome.contig.fa.masked.fa_R1.hmm",
            "3":f"AnnoSmk/{PREFIX}/R2/{PREFIX}.genome.contig.fa.masked.fa_R2.hmm"}

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

augustus_dict = {"1":f"AnnoSmk/{PREFIX}/pre_R0/autoAug/autoAugPred_hints/shells",
                 "2":f"AnnoSmk/{PREFIX}/R1/autoAug/autoAugPred_hints/shells",
                 "3":f"AnnoSmk/{PREFIX}/R2/autoAug/autoAugPred_hints/shells"}

def get_augustus_dir(wildcards):
    round = wildcards.round
    augustus_dir = augustus_dict[round]
    return augustus_dir

rule pre_opts:
    input:
        input_opt="AnnoSmk/{PREFIX}/maker_opts_0.ctl",
        snap_hmm=get_hmm,
        estgff="total_est.gff",
        pepgff="total_pep.gff",
        rmgff="repeat.gff",
        augustus_dir=get_augustus_dir
    params:
        round="{round}",
        augustus_species=get_augustus_species
    output:
        opts_file="AnnoSmk/{PREFIX}/R{round}/maker_opts{round}.ctl"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_pre_opts_R{round}.tsv"
    run:
        prepare_opts(input_opt=input.input_opt, estgff=input.estgff, pepgff=input.pepgff, 
                    rmgff=input.rmgff, round=params.round,
                    snap_hmm=input.snap_hmm,
                    augustus_species=params.augustus_species,
                    output_file=output.opts_file)

# def get_opts(wildcards):
#     opts_file = checkpoints.pre_exe.get(**wildcards).output[0]
#     round = glob_wildcards(f"AnnoSmk/{wildcards.PREFIX}/R{{round}}/maker_exe.ctl").round
#     # for num in round:
#     #     pren = int(num) - 1
#     opts_file = expand(rules.pre_opts.output,**wildcards,pren=lambda w: w-1)
#     return opts_file

checkpoint pre_exe:
    output:
        "AnnoSmk/{PREFIX}/R{round}/maker_exe.ctl",
        "AnnoSmk/{PREFIX}/R{round}/maker_bopts.ctl"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_pre_exe_R{round}.tsv"
    shell:
        '''
        maker -CTL 2>/dev/null
        cp maker_exe.ctl AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/
        cp maker_bopts.ctl AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/
        '''


rule run_maker:
    input:
        g="AnnoSmk/{PREFIX}/R{round}/{lane_number}.fa",
        opts="AnnoSmk/{PREFIX}/R{round}/maker_opts{round}.ctl",
        bopts="AnnoSmk/{PREFIX}/R{round}/maker_bopts.ctl",
        exe="AnnoSmk/{PREFIX}/R{round}/maker_exe.ctl"        
    output:
        log="AnnoSmk/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_master_datastore_index.log"
    threads:
        THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_run_maker_R{round}_{lane_number}.tsv"
    shell:
        '''
        bjobs -l $LSB_JOBID 
        /bin/rm -rf {wildcards.lane_number}.maker.output 2>> AnnoSmk/gossypium/R{wildcards.round}/{wildcards.lane_number}.log
        until [[ ! -e {wildcards.lane_number}.maker.output/{wildcards.lane_number}.db ]]
        do
            sleep 60
            /bin/rm -rf {wildcards.lane_number}.maker.output 2>> AnnoSmk/gossypium/R{wildcards.round}/{wildcards.lane_number}.log
        done
        echo "------START-------$(date)" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log
        echo "1.rm *.maker.output" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log               
        echo "2.run maker" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log 
        mpiexec -n {threads} maker -base {wildcards.lane_number} -genome {input.g} {input.opts} {input.bopts} {input.exe} -fix_nucleotides \
        > AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.error.maker 2>&1 & 
        sleep 300 
        until [[ -f "{wildcards.lane_number}.maker.output/{wildcards.lane_number}_master_datastore_index.log" ]] && \
        [[ `grep "Maker is now finished" AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.error.maker` ]]
        do
#            if [[ -f "{wildcards.lane_number}.maker.output/{wildcards.lane_number}_master_datastore_index.log" ]];then                
#                echo "$(date) but unfinish" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log 
#                sleep 60
            if [[ -f "{wildcards.lane_number}.maker.output/{wildcards.lane_number}_master_datastore_index.log" ]] && \
            [[ `grep "FAILED" {wildcards.lane_number}.maker.output/{wildcards.lane_number}_master_datastore_index.log` ]]; then 
                echo "$(date) get FAIL message. END now" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log 
#                /bin/rm -rf {wildcards.lane_number}.maker.output
#                wait
#                sleep 600
#                /bin/rm -rf {wildcards.lane_number}.maker.output
#                wait                                               
                bkill $LSB_JOBID
                exit 1 
            else
                echo "$(date) but unfinish" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log
                sleep 60                               
            fi
        done
        echo "$(date) finish maker" >> AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log 
        cp -rf {wildcards.lane_number}.maker.output AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/
        /bin/rm -rf {wildcards.lane_number}.maker.output
        /bin/rm -rf AnnoSmk/{PREFIX}/R{wildcards.round}/{wildcards.lane_number}.log 
        '''
#         maker -c {threads} -genome {input.g} {input.opts} {input.bopts} {input.exe} -fix_nucleotides

rule alt_log:
    input:
        "AnnoSmk/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_master_datastore_index.log"
    output:
        "AnnoSmk/{PREFIX}/R{round}/{lane_number}.maker.output/{lane_number}_total_master_datastore_index.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_alt_log_R{round}_{lane_number}.tsv"
    shell:
        '''
        cat {input} |sed "s/\t/\t{wildcards.lane_number}.maker.output\//" \
        > {output}
        '''

def get_log(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"AnnoSmk/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    log = expand(rules.alt_log.output, **wildcards, lane_number=lane_numbers)
    return log

rule merge_log:
    input:
        get_log
    output:
        "AnnoSmk/{PREFIX}/R{round}/total_master_datastore_index.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_merge_log_R{round}.tsv"
    shell:
        '''
        cat {input} > {output}
        '''

rule gff3_merge:
    input:
        rules.merge_log.output
    output:
        all_gff="AnnoSmk/{PREFIX}/R{round}/genome.all.gff",
        noseq_gff="AnnoSmk/{PREFIX}/R{round}/genome.all.noseq.gff",
        all_fasta="AnnoSmk/{PREFIX}/R{round}/total.all.maker.transcripts.fasta"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_gff3_merge_R{round}.tsv"
    shell:
        '''
        fasta_merge -d {input}
        gff3_merge -o {output.all_gff} -d {input}
        gff3_merge -n -o {output.noseq_gff} -d {input}
        cp total.all.maker.transcripts.fasta AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}
        '''

rule get_genome_maker_gff:
    input:
        rules.gff3_merge.output.noseq_gff
    output:
        "AnnoSmk/{PREFIX}/R{round}/genome.maker.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_get_genome_maker_gff_R{round}.tsv"
    shell:
        '''
        awk '$2=="maker"' {input} > {output}
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/*.maker.output
        '''

def get_fa(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"AnnoSmk/{wildcards.PREFIX}/sample/{{lane_number}}.fa").lane_number
    fa = expand(rules.cp.output.fa, **wildcards, lane_number=lane_numbers)
    return fa

rule get_ref_fa:
    input:
        get_fa
    output:
        "AnnoSmk/{PREFIX}/R{round}/ref.fa"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_get_ref_fa_R{round}.tsv"
    shell:
        "cat {input} > {output}"

rule maker2zff:
    input:
        "AnnoSmk/{PREFIX}/R{round}/genome.all.gff"
    output:
        ann="AnnoSmk/{PREFIX}/R{round}/genome.ann",
        dna="AnnoSmk/{PREFIX}/R{round}/genome.dna"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_maker2zff_R{round}.tsv"
    shell:
        '''
        maker2zff -x 0.25 -l 50 {input}
        mv genome.ann AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv genome.dna AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom1:
    input:
        ann="AnnoSmk/{PREFIX}/R{round}/genome.ann",
        dna="AnnoSmk/{PREFIX}/R{round}/genome.dna"
    output:
        "AnnoSmk/{PREFIX}/R{round}/gene-stats.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fathom1_R{round}.tsv"
    shell:
        "fathom -gene-stats {input.ann} {input.dna} >{output} 2>&1"

rule fathom2:
    input:
        ann="AnnoSmk/{PREFIX}/R{round}/genome.ann",
        dna="AnnoSmk/{PREFIX}/R{round}/genome.dna"
    output:
        "AnnoSmk/{PREFIX}/R{round}/validate.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fathom2_R{round}.tsv"
    shell:
        "fathom -validate {input.ann} {input.dna} >{output} 2>&1"

rule fathom3:
    input:
        ann="AnnoSmk/{PREFIX}/R{round}/genome.ann",
        dna="AnnoSmk/{PREFIX}/R{round}/genome.dna"
    output:
        uann="AnnoSmk/{PREFIX}/R{round}/uni.ann",
        udna="AnnoSmk/{PREFIX}/R{round}/uni.dna"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fathom3_R{round}.tsv"
    shell:
        '''
        fathom -categorize 1000 {input.ann} {input.dna}
        mv uni.ann AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv uni.dna AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule fathom4:
    input:
        uann="AnnoSmk/{PREFIX}/R{round}/uni.ann",
        udna="AnnoSmk/{PREFIX}/R{round}/uni.dna"
    output:
        exann="AnnoSmk/{PREFIX}/R{round}/export.ann",
        exdna="AnnoSmk/{PREFIX}/R{round}/export.dna"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fathom4_R{round}.tsv"
    shell:
        '''
        fathom -export 1000 -plus {input.uann} {input.udna}
        mv export.ann AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        mv export.dna AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

rule forge_hmm_assembler:
    input:
        exann="AnnoSmk/{PREFIX}/R{round}/export.ann",
        exdna="AnnoSmk/{PREFIX}/R{round}/export.dna"
    output:
        "AnnoSmk/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm"
    params:
        fo="{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_forge_hmm_assembler_R{round}.tsv"
    shell:
        '''
        cd AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}
        forge export.ann export.dna >forge.log 2>&1
        hmm-assembler.pl snap_trained . > {params.fo}
        cd -
        '''

rule fathom_to_genbank:
    input:
        uann="AnnoSmk/{PREFIX}/R{round}/uni.ann",
        udna="AnnoSmk/{PREFIX}/R{round}/uni.dna"
    output:
        "AnnoSmk/{PREFIX}/R{round}/augustus.gb"
    params:
        path=expand("{PATH}",PATH=PATH)
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fathom_to_genbank_R{round}.tsv"
    shell:
        '''
        perl {params.path}/bin/fathom_to_genbank.pl --annotation_file {input.uann} --dna_file {input.udna}  --genbank_file {output} --number 500
        '''
#fathom_to_genbank.pl文件需要修改perl的路径

rule perl_cat:
    input:
        "AnnoSmk/{PREFIX}/R{round}/augustus.gb"
    output:
        "AnnoSmk/{PREFIX}/R{round}/genbank_gene_list.txt"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_perl_cat_R{round}.tsv"
    shell:
        "cat.py {input} {output}"

rule get_subset_of_fastas:
    input:
        txt="AnnoSmk/{PREFIX}/R{round}/genbank_gene_list.txt",
        udna="AnnoSmk/{PREFIX}/R{round}/uni.dna"
    output:
        "AnnoSmk/{PREFIX}/R{round}/genbank_gene_seqs.fasta"
    params:
        path=expand("{PATH}",PATH=PATH)
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_get_subset_of_fastas_R{round}.tsv"
    shell:
        '''
        perl {params.path}/bin/get_subset_of_fastas.pl -l {input.txt} -f {input.udna} -o {output}
        '''

# rule randomSplit:
#     input:
#         "AnnoSmk/{PREFIX}/R{round}/augustus.gb"
#     output:
#         "AnnoSmk/{PREFIX}/R{round}/augustus.gb.test"
#     shell:
#         '''
#         randomSplit.pl {input} 250
#         '''


rule autoAugA:
    input:
        fasta="AnnoSmk/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gb="AnnoSmk/{PREFIX}/R{round}/augustus.gb",
        rna=augustus_cdna
    output:
        "AnnoSmk/{PREFIX}/R{round}/autoAug/autoAugPred_abinitio/shells/aug1",
        "AnnoSmk/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    params:
        path=expand("{PATH}",PATH=PATH)
    message:
        "If this step reports an error, you can delete autoAug and ~/.conda/envs/annotation/config/species/{wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct, and try again"
    threads: THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_autoAugA_R{round}.tsv"
    shell:
        '''
        /bin/rm -rf $(split_path.py $(which augustus))/../config/species/{PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct 
        cp -r {params.path}/bin/helpMod.pm $(split_path.py $(which augustus))
        /bin/rm -rf autoAug
        autoAug.pl --cpus={threads} --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --trainingset={input.gb} --cdna={input.rna} --noutr
        cd autoAug/autoAugPred_abinitio/shells
        let x=1
        while [ -e ./aug${{x}} ]
        do
            ./aug${{x}}
            let x=x+1
        done
        cd -
        cp -r autoAug AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/.
        '''

#如果有autoAug文件夹或者在/nfs/yanhui/.conda/envs/annotation/config/species下
#有species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct的文件夹都会报错

rule autoAugB:
    input:
        fasta="AnnoSmk/{PREFIX}/R{round}/genbank_gene_seqs.fasta",
        gff="AnnoSmk/{PREFIX}/R{round}/autoAug/hints/hints.E.gff"
    output:
        directory("AnnoSmk/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells")
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_autoAugB_R{round}.tsv"
    shell:
        '''
        autoAug.pl --species={wildcards.PREFIX}.genome.contig.fa.masked.fa_R{wildcards.round}_direct \
        --genome={input.fasta} --useexisting --hints={input.gff} \
        -v -v -v  --index=1
        cd autoAug/autoAugPred_hints/shells/
        let x=1
        while [ -e ./aug${{x}} ]
        do
            sed -i "s/extrinsic.E.cfg/extrinsic.M.RM.E.W.cfg/g" aug${{x}}
            ./aug${{x}}
            let x=x+1
        done
        cd ../../../
        cp -r autoAug/autoAugPred_hints AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/autoAug/.
        /bin/rm -rf autoAug
        '''

rule clean_memory:
    input:
        rules.merge_log.output,
        rules.autoAugB.output
    output:
        "AnnoSmk/{PREFIX}/R{round}/clean_memory.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_clean_memory_R{round}.tsv"
    shell:
        '''
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/*.maker.output > {output}
        '''

rule busco:
    input:
        rules.gff3_merge.output.all_fasta
    output:
        directory("AnnoSmk/{PREFIX}/R{round}/total.all.maker.transcripts.fasta.busco.embryophyta")
    params:
        dir_busco="total.all.maker.transcripts.fasta.busco.embryophyta"
    container:"docker://ezlabgva/busco:v5.4.3_cv1"
    threads: THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_busco_R{round}.tsv"
    shell:
        """
        cd AnnoSmk/{wildcards.PREFIX}/R{wildcards.round}/
        busco -f -c {threads} -m tran -i ../../../{input} -o {params.dir_busco} -l embryophyta_odb10
        """
rule AED:
    input:
        rules.gff3_merge.output.all_gff
    output:
        "AnnoSmk/{PREFIX}/R{round}/AED.csv"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_AED_R{round}.tsv"
    shell:
        "AED_cdf_generator.pl -b 0.025 {input} > {output}"


rule maker_gff_only:
    input:
        expand(rules.gff3_merge.output.all_gff,PREFIX=PREFIX, round=int(config["run_maker_Round"]))
    output:
        gff=expand("AnnoSmk/{PREFIX}/maker_gff_only.gff",PREFIX=PREFIX),
        fasta=expand("AnnoSmk/{PREFIX}/maker_fasta_only.fa",PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/maker_gff_only.tsv"
    shell:
        """
        awk -v line=$(awk '/##FASTA/{{print NR}}' {input}) '{{if(NR<line){{print $0}}}}' {input} > {output.gff}
        awk -v line=$(awk '/##FASTA/{{print NR}}' {input}) '{{if(NR>line){{print $0}}}}' {input} > {output.fasta}
        """


rule format_maker_gff:
    input:
        expand(rules.maker_gff_only.output.gff,PREFIX=PREFIX)
    output:
        pre_gff=expand("AnnoSmk/{PREFIX}_pre.gff",PREFIX=PREFIX),
        clu_0=temp("AnnoSmk/cluster0.gff")
    benchmark:
        "AnnoSmk/benchmark/format_maker_gff.tsv"
    shell:
        '''
        maker_gff.py {input} {output.pre_gff}
        getExon.py {output.pre_gff} {output.clu_0} > AnnoSmk/getExon.log
        /bin/rm -rf $(split_path.py $(which augustus))/../config/species/{PREFIX}.genome.contig.fa.masked.fa_R*_direct
        /bin/rm -f total_pep.gff
        /bin/rm -f repeat.gff
        /bin/rm -f total_est.gff
        '''

rule get_gene_overlap:
    input:
        rules.format_maker_gff.output.clu_0
    output:
        overlap_gene="AnnoSmk/overlap_gene.txt"
    benchmark:
        "AnnoSmk/benchmark/get_gene_overlap.tsv"
    shell:
        '''
        cd AnnoSmk
        maker_clean.sh
        /bin/rm -f cluster*
        cd -
        '''

rule clean_maker_gff:
    input:
        pre_gff=expand(rules.format_maker_gff.output.pre_gff,PREFIX=PREFIX),
        overlap_gene="AnnoSmk/overlap_gene.txt"
    output:
        expand("AnnoSmk/{PREFIX}.gff",PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/clean_maker_gff.tsv"
    shell:
        '''
        grep -vf {input.overlap_gene} {input.pre_gff} > {output}
        '''

















