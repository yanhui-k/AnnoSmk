import os
import glob

if config["FLNC"] == None:
    prefix_est=["rnaseq"]
    augustus_cdna="rnaseq.fasta"
elif config["FLNC"] != None:
    prefix_est=["rnaseq","flnc"]
    augustus_cdna="flnc.fasta"

localrules:
    cat_flnc_gff,
    mv,
    cat_est_rnaseq_d,
    cat_est_rnaseq_s,
    prep_genblast

number=[1,2,3,4,5,6,7,8]



rule hisat2build:
    input:
        fa=expand("AnnoSmk/{PREFIX}/{PREFIX}.fa",PREFIX=PREFIX)
    output:
        expand("AnnoSmk/{PREFIX}/evidence/{PREFIX}.{number}.ht2",PREFIX=PREFIX,number=number)
    benchmark:
        "AnnoSmk/benchmark/hisat2build.tsv"   
    shell:
        '''
        hisat2-build {input.fa} {PREFIX}
        mv {PREFIX}.*.ht2 AnnoSmk/{PREFIX}/evidence/
        '''

def get_all_rnaseq(wildcards):
    folder_d=os.path.join(PREFIX,"*_1.fastq")
    samples_d=glob.glob(folder_d)
    prefix_d=replace_list(samples_d,"_1.fastq","")
    prefix_d=replace_list(prefix_d,PREFIX,"")
    prefix_d=replace_list(prefix_d,"/","")
    prefix_d_out=[i+"_d" for i in prefix_d]
    folder_all=os.path.join(PREFIX,"*.fastq")
    samples_all=glob.glob(folder_all)
    prefix_all=replace_list(samples_all,".fastq","")
    prefix_all=replace_list(samples_all,"_1.fastq","")
    prefix_all=replace_list(samples_all,"_1.fastq","")
    prefix_all=replace_list(prefix_all,PREFIX,"")
    prefix_all=replace_list(prefix_all,"/","")
    prefix_s=list(set(prefix_all)-set(prefix_d))
    prefix_s_out=[i+"_s" for i in prefix_s]
    prefix_all=prefix_d_out+prefix_s_out
    all_fastq=expand("AnnoSmk/{PREFIX}/evidence/{sample}.fasta",PREFIX=PREFIX,sample=prefix_all)
    return all_fastq
     

rule fastp_d:
    input:
        r1="{PREFIX}/{sample}_1.fastq",
        r2="{PREFIX}/{sample}_2.fastq",
        check_download="AnnoSmk/download.log"
    output:
        r1="AnnoSmk/{PREFIX}/evidence/{sample}_1.clean.fq",
        r2="AnnoSmk/{PREFIX}/evidence/{sample}_2.clean.fq"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fastp_d_{sample}.tsv"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule hisat2_d:
    input:
        r1="AnnoSmk/{PREFIX}/evidence/{sample}_1.clean.fq",
        r2="AnnoSmk/{PREFIX}/evidence/{sample}_2.clean.fq",
        index=rules.hisat2build.output
    params:
        index="AnnoSmk/{PREFIX}/evidence/{PREFIX}"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq.bam"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_hisat2_d_{sample}.tsv"
    threads: THREADS
    shell:
        "hisat2 --mp 3,1 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} | samtools sort -@ {threads} -o {output}"

rule fastp_s:
    input:
        r1="{PREFIX}/{sample}_s.fastq",
        check_download="AnnoSmk/download.log"
    output:
        r1="AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_fastp_s_{sample}.tsv"
    shell:
        "fastp -i {input.r1} -o {output.r1}"

rule hisat2_s:
    input:
        r1="AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq",
        index=rules.hisat2build.output
    params:
        index="AnnoSmk/{PREFIX}/evidence/{PREFIX}"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq.bam"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_hisat2_s_{sample}.tsv"
    threads: THREADS
    shell:
        "hisat2 --mp 3,1 -p {threads} -x {params.index} -U {input.r1} | samtools sort -@ {threads} -o {output}"

rule trinity_d:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq.bam"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq_trinity/Trinity-GG.fasta"
    container:
        "docker://trinityrnaseq/trinityrnaseq"
    params:
        dir="AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq_trinity"
    threads: THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_trinity_d_{sample}.tsv"
    shell:
        '''
        Trinity --genome_guided_bam {input} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU {threads} \
          --output {params.dir}
        '''


rule trinity_s:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq.bam"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq_trinity/Trinity-GG.fasta"
    container:
        "docker://trinityrnaseq/trinityrnaseq"
    params:
        dir="AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq_trinity"
    threads: THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_trinity_d_{sample}.tsv"
    shell:
        '''
        Trinity --genome_guided_bam {input} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU {threads} \
          --output {params.dir}
        '''

#def get_bam1(wildcards):
#    samples_2=glob.glob("*_s.fastq.gz")
#    samples_2=[x.replace("_s.fastq.gz","") for x in samples_2]
#    bam1 = expand(rules.hisat2_s.output, sample=samples_2)
#    return bam1
#
#def get_bam2(wildcards):
#    samples_2=glob.glob("*_1.fastq.gz")
#    samples_2=[x.replace("_1.fastq.gz","") for x in samples_2]
#    bam2 = expand(rules.hisat2.output, sample=samples_2)
#    return bam2

rule mod_est_flnc:
    input:
        "{PREFIX}/{sample}_clean.fasta"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_clean.fasta"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_mod_est_flnc_{sample}.tsv"
    shell:
        '''
        cat {input} |sed "s/>/>{wildcards.sample}_/;s/\//_/; s/\./_/; s/\s.*//;" > {output}
        '''

rule cat_est_flnc:
    input:
        expand("{flnc}",flnc=config["FLNC"])
    output:
        expand("AnnoSmk/{PREFIX}/evidence/flnc.fasta",PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/cat_est_flnc.tsv"
    shell:
        "cat {input} >> {output}"

rule cp_flnc:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/flnc.fasta",PREFIX=PREFIX)
    output:
        "flnc.fasta"
    benchmark:
        "AnnoSmk/benchmark/cp_flnc.tsv"
    shell:
        "cp {input} {output}"

rule cat_est_rnaseq_d:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}_d.clean.fq_trinity/Trinity-GG.fasta"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_d.fasta"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_cat_est_rnaseq_d_{sample}.tsv"
    shell:
        '''
        cat {input} |sed "s/>/>{wildcards.sample}_/;s/\//_/; s/\./_/; s/\s.*//;" > {output}
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_1.clean.fq
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_2.clean.fq
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_d.clean.fq.bam
        find AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_d.clean.fq_trinity/ -type f -not -name "Trinity-GG.fasta" -delete
        '''

rule cat_est_rnaseq_s:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}_s.clean.fq_trinity/Trinity-GG.fasta"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}_s.fasta"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_cat_est_rnaseq_s_{sample}.tsv"
    shell:
        '''
        cat {input} |sed "s/>/>{wildcards.sample}_/;s/\//_/; s/\./_/; s/\s.*//;" > {output}
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_s.clean.fq
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_s.clean.fq.bam
        find AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.sample}_s.clean.fq_trinity/ -type f -not -name "Trinity-GG.fasta" -delete
        '''
        
rule merge_rnaseq:
    input:
        expand("{transcriptome}",transcriptome=config["TRAN"])
    output:
        expand("AnnoSmk/{PREFIX}/evidence/rnaseq.fasta", PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/merge_rnaseq.tsv"
    shell:
        '''
        cat {input} >> {output}
        /bin/rm -rf AnnoSmk/{PREFIX}/evidence/*.clean.fq_trinity
        '''

rule cp_rnaseq:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/rnaseq.fasta", PREFIX=PREFIX)
    output:
        "rnaseq.fasta"
    benchmark:
        "AnnoSmk/benchmark/cp_rnaseq.tsv"
    shell:
        '''
        cp {input} {output}
        '''

rule minimap2:
    input:
        fa1=expand("AnnoSmk/{PREFIX}/{PREFIX}.fa",PREFIX=PREFIX),
        fa2=expand("AnnoSmk/{PREFIX}/evidence/{{sample}}.fasta",PREFIX=PREFIX)
    output:
        expand("AnnoSmk/{PREFIX}/evidence/{{sample}}.fasta.bam",PREFIX=PREFIX)
    threads:
        THREADS
    shell:
        "minimap2 -t {threads} -C5 -ax splice {input.fa1} {input.fa2} |samtools view -F 256 -b > {output}"

rule bedtools:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.bam"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_bedtools_{sample}.tsv"
    shell:
        "bedtools bamtobed -split -i {input} > {output}"

rule awk:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.bed"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_awk_{sample}.tsv"
    shell:
        "awk '$3-$2 > 1' {input} > {output}"

rule bed2gff3:
    input:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.bed"
    output:
        "AnnoSmk/{PREFIX}/evidence/{sample}.fasta.rawgff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_bed2gff3_{sample}.tsv"
    shell:
        "gt bed_to_gff3 {input} | sort -k9,9 -k1,1 -k7,7 -k4,4n  > {output}"

checkpoint split_flnc_gff3:
    input:
        pep="AnnoSmk/{PREFIX}/evidence/flnc.fasta.rawgff"
    output:
        lane_dir = directory("AnnoSmk/{PREFIX}/evidence/split_flnc_gff3/")
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_split_flnc_gff3.tsv"
    shell:
        '''
        mkdir AnnoSmk/{wildcards.PREFIX}/evidence/split_flnc_gff3/
        split -l 300000 -d {input.pep} AnnoSmk/{wildcards.PREFIX}/evidence/split_flnc_gff3/flnc_ 
        '''

rule format_gt_flnc_gff_to_maker_gff:
    input:
        "AnnoSmk/{PREFIX}/evidence/split_flnc_gff3/flnc_{lane_number}"
    output:
        "AnnoSmk/{PREFIX}/evidence/split_flnc_gff3/flnc_{lane_number}.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_format_gt_flnc_gff_to_maker_gff_{lane_number}.tsv"
    shell:
        "maker.py format_gt_gff_to_maker_gff {input}"

rule format_gt_rnaseq_gff_to_maker_gff:
    input:
        "AnnoSmk/{PREFIX}/evidence/rnaseq.fasta.rawgff"
    output:
        "AnnoSmk/{PREFIX}/evidence/rnaseq_maker.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_format_gt_rnaseq_gff_to_maker_gff.tsv"
    shell:
        '''
        maker.py format_gt_gff_to_maker_gff {input}
        mv {input}.gff {output}
        '''

def get_flnc_gff3(wildcards):
    lane_dir = checkpoints.split_flnc_gff3.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"AnnoSmk/{wildcards.PREFIX}/evidence/split_flnc_gff3/flnc_{{lane_number}}").lane_number
    gff = expand(rules.format_gt_flnc_gff_to_maker_gff.output, **wildcards, lane_number=lane_numbers)
    return gff

rule cat_flnc_gff:
    input:
        get_flnc_gff3
    output:
        "AnnoSmk/{PREFIX}/evidence/flnc_maker.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_cat_flnc_gff.tsv"
    shell:
        "cat {input} >> {output}"

rule merge_gff:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/{sample}_maker.gff",sample=prefix_est, PREFIX=PREFIX)
    output:
        expand("AnnoSmk/{PREFIX}/evidence/total_est.gff", PREFIX=PREFIX)
    benchmark:
        "AnnoSmk/benchmark/merge_gff.tsv"
    shell:
        "cat {input} >> {output}"



