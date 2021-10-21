import pathlib
import random

import glob
import os

# rule all:
#     input:
#         "R1/total_master_datastore_index.log",
#         "R1/ref.fa",
#         "R1/genbank_gene_seqs.fasta",
#         "R1/augustus.gb.test",
#         "R1_autoAug/autoAugPred_hints/shells"

rule make_fasta:
    output:
        fasta = touch("altra.genome.fa.masked.fa")

checkpoint split_fasta:
    input:
        fasta = rules.make_fasta.output.fasta
    output:
        lane_dir = directory("sample/")
    shell:
        '''
        python -m iga.annotation.maker2 split_fasta1 {input} {output}    
        '''

rule mv:
    input:
        fasta = "sample/{lane_number}.fa"
    output:
        fa = "R1/{lane_number}.fa"
    shell:
        "cat {input} > {output}"

rule pre_opts:
    input:
        "opt.txt"
    output:
        "R1/maker_opts.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule pre_exe:
    input:
        "exe.txt"
    output:
        "R1/maker_exe.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule pre_bopts:
    input:
        "bopts.txt"
    output:
        "R1/maker_bopts.ctl"
    shell:
        '''
        cat {input} > {output}
        '''

rule run_maker:
    input:
        g="R1/{lane_number}.fa",
        opts="R1/maker_opts.ctl",
        bopts="R1/maker_bopts.ctl",
        exe="R1/maker_exe.ctl",
    output:
        "{lane_number}.maker.output/{lane_number}_master_datastore_index.log"
    shell:
        '''
        maker -genome {input.g} {input.opts} {input.bopts} {input.exe}
        '''

rule alt_log:
    input:
        rules.run_maker.output
    output:
        "{lane_number}.maker.output/{lane_number}_total_master_datastore_index.log"
    shell:
        '''
        cat {input} |sed "s/\t/\t{wildcards.lane_number}.maker.output\//" \
        > {output}
        '''

def get_log(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"sample/{{lane_number}}.fa").lane_number
    log = expand(rules.alt_log.output, **wildcards, lane_number=lane_numbers)
    return log

rule merge_log:
    input:
        get_log
    output:
        "R1/total_master_datastore_index.log"
    shell:
        "cat {input} > {output}"

rule gff3_merge:
    input:
        rules.merge_log.output
    output:
        all_gff="R1/genome.all.gff",
        noseq_gff="R1/genome.all.noseq.gff"
        #fasta="R1/total.all.maker.proteins.fasta"
    shell:
        '''
        fasta_merge -d {input}
        gff3_merge -o {output.all_gff} -d {input}
        gff3_merge -n -o {output.noseq_gff} -d {input}
        '''

rule get_genome_maker_gff:
    input:
        rules.gff3_merge.output.noseq_gff
    output:
        "R1/genome.maker.gff"
    shell:
        '''
        awk '$2=="maker"' {input} > {output}
        '''

def get_fa(wildcards):
    lane_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"sample/{{lane_number}}.fa").lane_number
    fa = expand(rules.mv.output, **wildcards, lane_number=lane_numbers)
    return fa

rule get_ref_fa:
    input:
        get_fa
    output:
        "R1/ref.fa"
    shell:
        "cat {input} > {output}"

rule maker2zff:
    input:
        "R1/genome.all.gff"
    output:
        ann="R1/genome.ann",
        dna="R1/genome.dna"
    shell:
        '''
        maker2zff -x 0.25 -l 50 {input}
        mv genome.ann R1/.
        mv genome.dna R1/.
        '''

rule fathom1:
    input:
        ann="R1/genome.ann",
        dna="R1/genome.dna"
    output:
        "R1/gene-stats.log"
    shell:
        "fathom -gene-stats {input.ann} {input.dna} >{output} 2>&1"

rule fathom2:
    input:
        ann="R1/genome.ann",
        dna="R1/genome.dna"
    output:
        "R1/validate.log"
    shell:
        "fathom -validate {input.ann} {input.dna} >{output} 2>&1"

rule fathom3:
    input:
        ann="R1/genome.ann",
        dna="R1/genome.dna"
    output:
        uann="R1/uni.ann",
        udna="R1/uni.dna"
    shell:
        '''
        fathom -categorize 1000 {input.ann} {input.dna}
        mv uni.ann R1/.
        mv uni.dna R1/.
        '''

rule fathom4:
    input:
        uann="R1/uni.ann",
        udna="R1/uni.dna"
    output:
        exann="R1/export.ann",
        exdna="R1/export.dna"
    shell:
        '''
        fathom -export 1000 -plus {input.uann} {input.udna}
        mv export.ann R1/.
        mv export.ann R1/.
        '''

rule forge:
    input:
        exann="R1/export.ann",
        exdna="R1/export.dna"
    output:
        "R1/forge.log"
    shell:
        "forge {input.exann} {input.exdna} >{output} 2>&1"

rule hmm_assembler:
    input:
        files=rules.forge.output,
        dir="R1"
    output:
        "R1/altra.genome.contig.fa.masked.fa_R1.hmm"
    shell:
        '''
        hmm-assembler.pl snap_trained {input.dir} > {output}
        '''

rule fathom_to_genbank:
    input:
        uann="R1/uni.ann",
        udna="R1/uni.dna"
    output:
        "R1/augustus.gb"
    shell:
        '''
        fathom_to_genbank.pl --annotation_file {input.uann} --dna_file {input.udna}  --genbank_file {output} --number 500
        '''
#fathom_to_genbank.pl文件需要修改perl的路径

rule perl_cat:
    input:
        "R1/augustus.gb"
    output:
        "R1/genbank_gene_list.txt"
    script:
        "cat.py"

rule get_subset_of_fastas:
    input:
        txt="R1/genbank_gene_list.txt",
        udna="R1/uni.dna"
    output:
        "R1/genbank_gene_seqs.fasta"
    shell:
        '''
        export PATH=/nfs/yanhui/anntest/bundle/lh_bin/:$PATH
        get_subset_of_fastas.pl -l {input.txt} -f {input.udna} -o {output}
        '''

rule randomSplit:
    input:
        "R1/augustus.gb"
    output:
        "R1/augustus.gb.test"
    shell:
        '''
        randomSplit.pl {input} 250
        '''

rule autoAugA:
    input:
        fasta="R1/genbank_gene_seqs.fasta",
        gb="R1/augustus.gb",
        cdna="/nfs/yanhui/anntest/02.maker/flnc.fasta"
    output:
        "R1_autoAug/autoAugPred_abinitio/shells/aug1",
        "R1_autoAug/hints/hints.E.gff"
    shell:
        '''
        autoAug.pl --species=altra.genome.contig.fa.masked.fa_R1_direct \
        --genome={input.fasta} --trainingset={input.gb} --cdna={input.cdna} --noutr
        cp autoAug R1_autoAug
        '''

#如果有autoAug文件夹或者在/nfs/yanhui/.conda/envs/repeat/config/species下
#有species=altra.genome.contig.fa.masked.fa_R1_direct的文件夹都会报错

rule autoAugB:
    input:
        fasta="R1/genbank_gene_seqs.fasta",
        gff="R1_autoAug/hints/hints.E.gff"
    output:
        "R1_autoAug/autoAugPred_hints/shells"
    shell:
        '''
        autoAug.pl --species=altra.genome.contig.fa.masked.fa_R1_direct \
        --genome={input.fasta} --useexisting --hints={input.gff} \
        -v -v -v  --index=1
        cp autoAug/autoAugPred_hints/shells {output}
        '''

# rule clean:
#     input:
#         file="autoAug/autoAugPred_hints/shells"
#     output:
#         "R1/autoAug"
#     params:
#         dir="autoAug"
#     shell:
#         "mv -b {params.dir} R1/."

#rule busco:
#    input:
#        "R1/total.all.maker.proteins.fasta"
#    output:
#        touch("R1/total.all.maker.proteins.fasta.busco.embryophyta.v4.1.2")
#    conda:
#        "/nfs/yanhui/busco.yml"
#    shell:
#        """
#        busco -f -c 64 -m prot -i {input} -o {output} -l embryophyta_odb10
#        """
























