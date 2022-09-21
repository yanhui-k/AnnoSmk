import os
import glob

def replace_list(list,a,b):
    new_list = []
    while list:
        new_item = list.pop().replace(a,b,1)
        new_list.append(new_item)
    return new_list

folder_flnc = os.path.join(PREFIX,"*_subreads.fastq.gz")
samples_flnc=glob.glob(folder_flnc)
if samples_flnc == []:
    prefix_est=["rnaseq"]
    augustus_cdna="annotation_smk/{PREFIX}/evidence/rnaseq.fasta"
elif samples_flnc != []:
    prefix_est=["rnaseq","flnc"]
    augustus_cdna="annotation_smk/{PREFIX}/evidence/flnc.fasta"

localrules:
    cat_flnc_gff,
    mv

folder_d = os.path.join(PREFIX,"*_1.fastq.gz")
samples_d=glob.glob(folder_d)
prefix_d=replace_list(samples_d,"_1.fastq.gz","")
prefix_d=replace_list(prefix_d,PREFIX,"")
prefix_d=replace_list(prefix_d,"/","")
prefix_d=[i+"_d" for i in prefix_d]

folder_s = os.path.join(PREFIX,"*_s.fastq.gz")
samples_s=glob.glob(folder_s)
prefix_s=replace_list(samples_s,"_s.fastq.gz","")
prefix_s=replace_list(prefix_s,PREFIX,"")
prefix_s=replace_list(prefix_s,"/","")
prefix_s=[i+"_s" for i in prefix_s]
prefix_all=prefix_d+prefix_s
number = [1,2,3,4,5,6,7,8]

rule BuildDatabase:
    input:
        fa=expand("{REF}",REF=REF)
    output:
        "{PREFIX}.nsq",
        "{PREFIX}.nog",
        "{PREFIX}.nni",
        "{PREFIX}.nnd",
        "{PREFIX}.nhr",
        "{PREFIX}.translation"
    shell:
        "BuildDatabase -name {PREFIX} -engine ncbi {input.fa} "

rule RepeatModeler:
    input:
        fa=expand("{REF}",REF=REF),
        nsq="{PREFIX}.nsq",
        nog="{PREFIX}.nog"
    output:
        famfa="annotation_smk/{PREFIX}/evidence/{PREFIX}-families.fa",
        famstk="annotation_smk/{PREFIX}/evidence/{PREFIX}-families.stk"
    params:
        db="{PREFIX}"
    threads: THREADS
    shell:
        '''
        RepeatModeler -pa {threads} -engine ncbi -database {params.db}
        cp {wildcards.PREFIX}-families.fa annotation_smk/{wildcards.PREFIX}/evidence/
        cp {wildcards.PREFIX}-families.stk annotation_smk/{wildcards.PREFIX}/evidence/
        '''

rule RepeatMasker:
    input:
        fa=expand("{REF}",REF=REF),
        famfa=rules.RepeatModeler.output.famfa
    output:
        catgz="annotation_smk/{PREFIX}/evidence/{PREFIX}.cat.gz",
        out="annotation_smk/{PREFIX}/evidence/{PREFIX}.out"
    params:
        dir="annotation_smk/{PREFIX}/evidence/"
    threads: THREADS
    shell:
        '''
        RepeatMasker -lib {input.famfa} {input.fa} -pa {threads} -dir {params.dir}
        mv annotation_smk/{PREFIX}/evidence/*.cat.gz annotation_smk/{PREFIX}/evidence/{PREFIX}.cat.gz
        mv annotation_smk/{PREFIX}/evidence/*.out annotation_smk/{PREFIX}/evidence/{PREFIX}.out
        mv annotation_smk/{PREFIX}/evidence/*.masked annotation_smk/{PREFIX}/evidence/{PREFIX}.masked
        mv annotation_smk/{PREFIX}/evidence/*.tbl annotation_smk/{PREFIX}/evidence/{PREFIX}.tbl
        '''

#rule post_repeatmasker_gunzip:
#    input:
#        catgz="annotation_smk/{PREFIX}/evidence/{REF}.cat.gz"
#    output:
#        cat="annotation_smk/{PREFIX}/evidence/{REF}.cat"
#    shell:
#        "gunzip {input.catgz} "
#
#rule ProcessRepeats:
#    input:
#        cat=rules.post_repeatmasker_gunzip.output.cat
#    output:
#        out="annotation_smk/{PREFIX}/evidence/{REF}.out",
#        tbl="annotation_smk/{PREFIX}/evidence/{REF}.tbl"
#    shell:
#        "ProcessRepeats -species Viridiplantae {input.cat}"

rule rmOutToGFF3:
    input:
        expand("annotation_smk/{PREFIX}/evidence/{PREFIX}.out", PREFIX=PREFIX)
    output:
        "annotation_smk/{PREFIX}/evidence/full_mask.gff3"
    shell:
        "rmOutToGFF3.pl {input} > {output}"
#Can't locate CrossmatchSearchEngine.pm in @INC
#Can't locate SearchEngineI.pm in @INC
#Can't locate SearchResultCollection.pm in @INC
#Can't locate SearchResult.pm in @INC
#Can't locate Matrix.pm in @INC
#Can't locate ArrayList.pm in @INC
#Can't locate ArrayListIterator.pm
#download from https://github.com/rmhubley/RepeatMasker

rule isolate_complex_repeats:
    input:
        rules.rmOutToGFF3.output
    output:
        "annotation_smk/{PREFIX}/evidence/full_mask.complex.gff3"
    shell:
        'grep -v -e "Satellite" -e ")n" -e "-rich" {input} > {output}'

rule reformat_to_work_with_MAKER:
    input:
        rules.isolate_complex_repeats.output
    output:
        "annotation_smk/{PREFIX}/evidence/repeat.gff"
    shell:
        "reformat_to_work_with_MAKER.py {input} {output}"

rule bedtools1:
    input:
        fa=expand("{REF}",REF=REF),
        gff3=rules.reformat_to_work_with_MAKER.output
    output:
        masked_fa="annotation_smk/{PREFIX}/evidence/genome.fa.masked.fa"
    log:
        "annotation_smk/{PREFIX}/evidence/bedtools.log"
    shell:
        '''
        bedtools maskfasta -soft -fi {input.fa} -bed {input.gff3} -fo {output.masked_fa}
        '''

rule hisat2build:
    input:
        fa = expand("{REF}",REF=REF)
    output:
        expand("annotation_smk/{PREFIX}/evidence/{PREFIX}.{number}.ht2",PREFIX=PREFIX,number=number)
    shell:
        '''
        hisat2-build {input} {PREFIX}
        mv {PREFIX}.*.ht2 annotation_smk/{PREFIX}/evidence/
        '''
        

rule fastp:
    input:
        r1 = "{PREFIX}/{sample}_1.fastq.gz",
        r2 = "{PREFIX}/{sample}_2.fastq.gz"
    output:
        r1 = "annotation_smk/{PREFIX}/evidence/{sample}_1.clean.fq.gz",
        r2 = "annotation_smk/{PREFIX}/evidence/{sample}_2.clean.fq.gz"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule hisat2:
    input:
        r1="annotation_smk/{PREFIX}/evidence/{sample}_1.clean.fq.gz",
        r2="annotation_smk/{PREFIX}/evidence/{sample}_2.clean.fq.gz",
        index=rules.hisat2build.output
    params:
        index="annotation_smk/{PREFIX}/evidence/{PREFIX}"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}_d.clean.fq.gz.bam"
    threads: THREADS
    shell:
        "hisat2 --rna-strandness RF --mp 3,1 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} | samtools sort -@ 30 -o {output}"

rule fastp_s:
    input:
        r1 = "{PREFIX}/{sample}_s.fastq.gz"
    output:
        r1 = "annotation_smk/{PREFIX}/evidence/{sample}_s.clean.fq.gz"
    shell:
        "fastp -i {input.r1} -o {output.r1}"

rule hisat2_s:
    input:
        r1="annotation_smk/{PREFIX}/evidence/{sample}_s.clean.fq.gz",
        index=rules.hisat2build.output
    params:
        index="annotation_smk/{PREFIX}/evidence/{PREFIX}"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}_s.clean.fq.gz.bam"
    threads: THREADS
    shell:
        "hisat2 --rna-strandness RF --mp 3,1 -p 30 -x tair -U {input.r1} | samtools sort -@ 30 -o {output}"


rule trinity:
    input:
        "annotation_smk/{PREFIX}/evidence/{sample}.clean.fq.gz.bam"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}.clean.fq.gz_trinity/Trinity-GG.fasta"
    params:
        dir="annotation_smk/{PREFIX}/evidence/{sample}.clean.fq.gz_trinity"
    threads: THREADS
    shell:
        '''
        Trinity --SS_lib_type RF \
          --genome_guided_bam {input} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU {threads} \
          --output {params.dir}
        '''
#
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

rule cat_est_flnc:
    input:
        expand("{flnc}",flnc=samples_flnc)
    output:
        "annotation_smk/{PREFIX}/evidence/flnc.fasta"
    shell:
        "maker.py cat_est {input} > {output}"

rule cat_est_rnaseq:
    input:
        expand("annotation_smk/{PREFIX}/evidence/{sample}.clean.fq.gz_trinity/Trinity-GG.fasta", PREFIX=PREFIX, sample=prefix_all)
    output:
        expand("annotation_smk/{PREFIX}/evidence/rnaseq.fasta", PREFIX=PREFIX)
    shell:
        "maker.py cat_est {input} > {output}"

rule minimap2:
    input:
        fa1 = expand("{REF}",REF=REF),
        fa2 = "annotation_smk/{PREFIX}/evidence/{sample}.fasta"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.bam"
    threads:
        THREADS
    shell:
        "minimap2 -t {threads} -C5 -ax splice {input.fa1} {input.fa2} |samtools view -F 256 -b > {output}"

rule bedtools:
    input:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.bam"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    shell:
        "bedtools bamtobed -split -i {input} > {output}"

rule awk:
    input:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.bed"
    shell:
        "awk '$3-$2 > 1' {input} > {output}"

rule bed2gff3:
    input:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.bed"
    output:
        "annotation_smk/{PREFIX}/evidence/{sample}.fasta.rawgff"
    shell:
        "gt bed_to_gff3 {input} | sort -k9,9 -k1,1 -k7,7 -k4,4n  > {output}"

checkpoint split_flnc_gff3:
    input:
        pep = "annotation_smk/{PREFIX}/evidence/flnc.fasta.rawgff"
    output:
        lane_dir = directory("annotation_smk/{PREFIX}/evidence/split_flnc_gff3/")
    shell:
        '''
        mkdir annotation_smk/{wildcards.PREFIX}/evidence/split_flnc_gff3/
        split -l 300 -d {input} annotation_smk/{wildcards.PREFIX}/evidence/split_flnc_gff3/flnc_ 
        '''

rule format_gt_flnc_gff_to_maker_gff:
    input:
        "annotation_smk/{PREFIX}/evidence/split_flnc_gff3/flnc_{lane_number}"
    output:
        "annotation_smk/{PREFIX}/evidence/split_flnc_gff3/flnc_{lane_number}.gff"
    shell:
        "maker.py format_gt_gff_to_maker_gff {input}"

rule format_gt_rnaseq_gff_to_maker_gff:
    input:
        "annotation_smk/{PREFIX}/evidence/rnaseq.fasta.rawgff"
    output:
        "annotation_smk/{PREFIX}/evidence/rnaseq_maker.gff"
    shell:
        '''
        maker.py format_gt_gff_to_maker_gff {input}
        mv {input}.gff {output}
        '''

def get_flnc_gff3(wildcards):
    lane_dir = checkpoints.split_flnc_gff3.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"annotation_smk/{wildcards.PREFIX}/evidence/split_flnc_gff3/flnc_{{lane_number}}").lane_number
    gff = expand(rules.format_gt_flnc_gff_to_maker_gff.output, **wildcards, lane_number=lane_numbers)
    return gff

rule cat_flnc_gff:
    input:
        get_flnc_gff3
    output:
        "annotation_smk/{PREFIX}/evidence/flnc_maker.gff"
    shell:
        "cat {input} >> {output}"

rule merge_gff:
    input:
        expand("annotation_smk/{PREFIX}/evidence/{sample}_maker.gff",sample=prefix_est, PREFIX=PREFIX)
    output:
        expand("annotation_smk/{PREFIX}/evidence/total_est.gff", PREFIX=PREFIX)
    shell:
        "cat {input} >> {output}"

rule pre_pep:
    output:
        pep=expand("{PEP}",PEP=PEP)
    shell:
        '''
        touch {output.pep}
        '''
rule prepep:
    input:
        pep = rules.pre_pep.output
    output:
        "annotation_smk/{PREFIX}/evidence/0.fa"
    shell:
        "modfasta.py {input} {output}"

checkpoint split_pep:
    input:
        pep = "annotation_smk/{PREFIX}/evidence/0.fa"
    output:
        lane_dir = directory("annotation_smk/{PREFIX}/evidence/sample/")
    params:
        size=600000
    shell:
        "split_fasta.py {input} {output} {params.size}"

rule mv:
    input:
        fasta = "annotation_smk/{PREFIX}/evidence/sample/{lane_number}.fa"
    output:
        fa = "annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa"
    shell:
        "cp {input} {output}"

rule prep_genblast:
    output:
        "annotation_smk/{PREFIX}/evidence/pep/alignscore.txt"
    params:
        ref="genome.fa.masked.fa",
        dir="annotation_smk/{PREFIX}/evidence/pep"
    shell:
        '''
        sh {input} > {output} 
        cd {params.dir}
        ln -s `which formatdb` 
        ln -s `which blastall`
        ln -s ../{params.ref}
        cd -
        '''

rule genblast:
    input:
        QRY="annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa",
        REF="annotation_smk/{PREFIX}/evidence/genome.fa.masked.fa",
        alignscore="annotation_smk/{PREFIX}/evidence/pep/alignscore.txt"
    output:
        gff="annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.gff",
        pro="annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.pro"
    params:
        dir="annotation_smk/{PREFIX}/evidence/pep",
        ref="genome.fa.masked.fa",
        qry="{lane_number}.fa"
    shell:
        """
        cd {params.dir}
        genblastG -q {params.qry} -t {params.ref} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
            -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {params.qry}
        wait
        until [ -s {params.qry}*.DNA ]; do
            rm {params.qry}*.blast
            rm {params.qry}*.blast.report
            genblastG -q {params.qry} -t {params.ref} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
            -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {params.qry}
        done
        mv {params.qry}_*.gff {params.qry}.gff
        mv {params.qry}_*.pro {params.qry}.pro       
        cd -
        """

rule filter_genblast:
    input:
        rules.genblast.output.gff
    output:
        "annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.slim.genblast.gff"
    shell:
        "genblast2.py filter_genblast {input} > {output}"

rule filter_early_stop:
    input:
        rules.genblast.output.pro
    output:
        "annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.genblast.noearly_stop.id"
    shell:
        "genblast2.py filter_early_stop {input} > {output}"

rule selectGFF:
    input:
        stopid=rules.filter_early_stop.output,
        slimgff=rules.filter_genblast.output
    output:
        "annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.filter.genblast.gff"
    shell:
        "selectGFF.pl {input.stopid} {input.slimgff} > {output}"

rule sed_gff:
    input:
        rules.selectGFF.output
    output:
        "annotation_smk/{PREFIX}/evidence/pep/{lane_number}.fa.final.gff"
    shell:
        "sed 's/transcript/protein_match/; s/coding_exon/match_part/' {input} > {output}"

def get_genblast_gff(wildcards):
    lane_dir = checkpoints.split_pep.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"annotation_smk/{wildcards.PREFIX}/evidence/sample/{{lane_number}}.fa").lane_number
    gff = expand(rules.sed_gff.output, **wildcards, lane_number=lane_numbers)
    return gff

rule merge_genblast_gff:
    input:
        get_genblast_gff
    output:
        "annotation_smk/{PREFIX}/evidence/genblast.gff"
    shell:
        "cat {input} > {output}"

