import pathlib
import random

import glob
import os

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

rule pre_opts:
    input:
        "config/opt.txt"
    output:
        "result/{PREFIX}/R{round}/maker_opts.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule pre_exe:
    input:
        "config/exe.txt"
    output:
        "result/{PREFIX}/R{round}/maker_exe.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule pre_bopts:
    input:
        "config/bopts.txt"
    output:
        "result/{PREFIX}/R{round}/maker_bopts.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule run_maker:
    input:
        g="result/{PREFIX}/R{round}/{lane_number}.fa",
        opts="result/{PREFIX}/R{round}/maker_opts.ctl",
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
        mv export.ann result/{wildcards.PREFIX}/R{wildcards.round}/.
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
        cp -r autoAug result/{wildcards.PREFIX}/R{wildcards.round}/autoAug
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
        mv -f autoAug/autoAugPred_hints R1_autoAug
        rm -r autoAug
        '''

rule busco:
    input:
        rules.gff3_merge.output.all_fasta
    output:
        touch("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta")
    params:
        dir_busco="total.all.maker.proteins.fasta.busco.embryophyta"
    shell:
        """
        cd result/{wildcards.PREFIX}/R{wildcards.round}/
        run_busco -f -c 64 -m prot -i ../{input} -o {params.dir_busco} -l embryophyta_odb10
        """
























