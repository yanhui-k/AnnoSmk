import os
import glob

folder_flnc = os.path.join(PREFIX,"*_subreads.fastq.gz")
samples_flnc=glob.glob(folder_flnc)
#samples_minimap=["flnc","rnaseq"]
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
        famfa="result/{PREFIX}/evidence/{PREFIX}-families.fa",
        famstk="result/{PREFIX}/evidence/{PREFIX}-families.stk"
    params:
        db="{PREFIX}"
    threads: THREADS
    shell:
        '''
        RepeatModeler -pa {threads} -engine ncbi -database {params.db}
        cp {wildcards.PREFIX}-families.fa result/{wildcards.PREFIX}/evidence/
        cp {wildcards.PREFIX}-families.stk result/{wildcards.PREFIX}/evidence/
        '''

rule RepeatMasker:
    input:
        fa=expand("{REF}",REF=REF),
        famfa=rules.RepeatModeler.output.famfa
    output:
        catgz="result/{PREFIX}/evidence/{PREFIX}.cat.gz",
        out="result/{PREFIX}/evidence/{PREFIX}.out"
    params:
        dir="result/{PREFIX}/evidence/"
    threads: THREADS
    shell:
        '''
        RepeatMasker -lib {input.famfa} {input.fa} -pa {threads} -dir {params.dir}
        mv result/{PREFIX}/evidence/*.cat.gz result/{PREFIX}/evidence/{PREFIX}.cat.gz
        mv result/{PREFIX}/evidence/*.out result/{PREFIX}/evidence/{PREFIX}.out
        mv result/{PREFIX}/evidence/*.masked result/{PREFIX}/evidence/{PREFIX}.masked
        mv result/{PREFIX}/evidence/*.tbl result/{PREFIX}/evidence/{PREFIX}.tbl
        '''

#rule post_repeatmasker_gunzip:
#    input:
#        catgz="result/{PREFIX}/evidence/{REF}.cat.gz"
#    output:
#        cat="result/{PREFIX}/evidence/{REF}.cat"
#    shell:
#        "gunzip {input.catgz} "
#
#rule ProcessRepeats:
#    input:
#        cat=rules.post_repeatmasker_gunzip.output.cat
#    output:
#        out="result/{PREFIX}/evidence/{REF}.out",
#        tbl="result/{PREFIX}/evidence/{REF}.tbl"
#    shell:
#        "ProcessRepeats -species Viridiplantae {input.cat}"

rule rmOutToGFF3:
    input:
        expand("result/{PREFIX}/evidence/{PREFIX}.out", PREFIX=PREFIX)
    output:
        "result/{PREFIX}/evidence/full_mask.gff3"
    shell:
        "perl_module/rmOutToGFF3.pl {input} > {output}"
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
        "result/{PREFIX}/evidence/full_mask.complex.gff3"
    shell:
        'grep -v -e "Satellite" -e ")n" -e "-rich" {input} > {output}'

rule reformat_to_work_with_MAKER:
    input:
        rules.isolate_complex_repeats.output
    output:
        "result/{PREFIX}/evidence/repeat.gff"
    script:
        "../bin/reformat_to_work_with_MAKER.py"

rule bedtools1:
    input:
        fa=expand("{REF}",REF=REF),
        gff3=rules.reformat_to_work_with_MAKER.output
    output:
        masked_fa="result/{PREFIX}/evidence/genome.fa.masked.fa"
    log:
        "result/{PREFIX}/evidence/bedtools.log"
    shell:
        '''
        bedtools maskfasta -soft -fi {input.fa} -bed {input.gff3} -fo {output.masked_fa}
        '''

rule hisat2build:
    input:
        fa = expand("{REF}",REF=REF)
    output:
        expand("result/{PREFIX}/evidence/{PREFIX}.{number}.ht2",PREFIX=PREFIX,number=number)
    shell:
        '''
        hisat2-build {input} {PREFIX}
        mv {PREFIX}.*.ht2 result/{PREFIX}/evidence/
        '''
        

rule fastp:
    input:
        r1 = "{PREFIX}/{sample}_1.fastq.gz",
        r2 = "{PREFIX}/{sample}_2.fastq.gz"
    output:
        r1 = "result/{PREFIX}/evidence/{sample}_1.clean.fq.gz",
        r2 = "result/{PREFIX}/evidence/{sample}_2.clean.fq.gz"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule hisat2:
    input:
        r1="result/{PREFIX}/evidence/{sample}_1.clean.fq.gz",
        r2="result/{PREFIX}/evidence/{sample}_2.clean.fq.gz",
        index=rules.hisat2build.output
    params:
        index="result/{PREFIX}/evidence/{PREFIX}"
    output:
        "result/{PREFIX}/evidence/{sample}.clean.fq.gz.bam"
    threads: THREADS
    shell:
        "hisat2 --rna-strandness RF --mp 3,1 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} | samtools sort -@ 30 -o {output}"

rule trinity:
    input:
        "result/{PREFIX}/evidence/{sample}.clean.fq.gz.bam"
    output:
        ggfa="result/{PREFIX}/evidence/{sample}.clean.fq.gz_trinity/Trinity-GG.fasta"
    params:
        dir="result/{PREFIX}/evidence/{sample}.clean.fq.gz_trinity"
    threads: THREADS
    shell:
        '''
        Trinity --SS_lib_type RF \
          --genome_guided_bam {input} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU {threads} \
          --output {params.dir}
        '''

#def get_bam(wildcards):
#    folder_2 = os.path.join(PREFIX,"*_1.fastq.gz")
#    samples_2=glob.glob(folder_2)
#    samples_2=[x.replace("_1.fastq.gz","") for x in samples_2]
#    samples_2=[x.replace("/","") for x in samples_2]
#    samples_2=[x.replace(PREFIX,"") for x in samples_2]
#    bam = expand(rules.hisat2.output, PREFIX=PREFIX, sample=samples_2)
#    return bam

def get_ggfa(wildcards):
    samples_2=glob.glob("*_1.fastq.gz")
    samples_2=[x.replace("_1.fastq.gz","") for x in samples_2]
    ggfa = expand(rules.trinity.output.ggfa, sample=samples_2)
    return ggfa

#rule samtools_merge:
#    input:
#        get_bam
#    output:
#        "result/{PREFIX}/evidence/merged.bam"
#    shell:
#        "samtools merge -@ 40 {output} {input}"
#
#rule stringtie:
#    input:
#        expand(rules.samtools_merge.output, PREFIX=PREFIX)
#    output:
#        "result/{PREFIX}/evidence/merged.gtf"
#    shell:
#        "stringtie -p 64 -o {output} {input}"
#
#rule gffread:
#    input:
#        gtf = rules.stringtie.output,
#        REF = expand("{REF}",REF=REF)
#    output:
#        "result/{PREFIX}/evidence/rnaseq.fasta"
#    shell:
#        "gffread -w {output} -g {input.REF} {input.gtf}"

rule cat_est_flnc:
    input:
        expand("{flnc}",flnc=samples_flnc)
    output:
        "result/{PREFIX}/evidence/flnc.fasta"
    shell:
        "python -m iga.annotation.maker cat_est {input} > {output}"

rule cat_est_rnaseq:
    input:
        get_ggfa
    output:
        "result/{PREFIX}/evidence/rnaseq.fasta"
    shell:
        "python -m iga.annotation.maker cat_est {input} > {output}"

rule minimap2:
    input:
        fa1 = expand("{REF}",REF=REF),
        fa2 = "result/{PREFIX}/evidence/{sample}.fasta"
    output:
        "result/{PREFIX}/evidence/{sample}.fasta.bam"
    threads:
        THREADS
    shell:
        "minimap2 -t {threads} -C5 -ax splice {input.fa1} {input.fa2} |samtools view -F 256 -b > {output}"

rule bedtools:
    input:
        "result/{PREFIX}/evidence/{sample}.fasta.bam"
    output:
        "result/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    shell:
        "bedtools bamtobed -split -i {input} > {output}"

rule awk:
    input:
        "result/{PREFIX}/evidence/{sample}.fasta.raw.bed"
    output:
        "result/{PREFIX}/evidence/{sample}.fasta.bed"
    shell:
        "awk '$3-$2 > 1' {input} > {output}"

rule bed2gff3:
    input:
        "result/{PREFIX}/evidence/{sample}.fasta.bed"
    output:
        "result/{PREFIX}/evidence/{sample}.fasta.rawgff"
    shell:
        "gt bed_to_gff3 {input} | sort -k9,9 -k1,1 -k7,7 -k4,4n  > {output}"

rule format_gt_gff_to_maker_gff:
    input:
        rules.bed2gff3.output
    output:
        "result/{PREFIX}/evidence/{sample}.fasta.rawgff.gff"
    shell:
        "python -m iga.annotation.maker format_gt_gff_to_maker_gff {input}"

rule cat_gff:
    input:
        expand("result/{PREFIX}/evidence/{sample}.fasta.rawgff.gff",PREFIX=PREFIX,sample=["rnaseq","flnc"])
    output:
        expand("result/{PREFIX}/evidence/total_est.gff", PREFIX=PREFIX)
    shell:
        "cat {input} > {output}"

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
        "result/{PREFIX}/evidence/0.fa"
    script:
        "../bin/modfasta.py"

checkpoint split_pep:
    input:
        pep = "result/{PREFIX}/evidence/0.fa"
    output:
        lane_dir = directory("result/{PREFIX}/evidence/sample/")
    params:
        size=600000
    script:
        "../bin/split_fasta1.py"

rule mv:
    input:
        fasta = "result/{PREFIX}/evidence/sample/{lane_number}.fa"
    output:
        fa = "result/{PREFIX}/evidence/pep/{lane_number}.fa"
    shell:
        "cp {input} {output}"

rule prep_genblast:
    input:
        "config/alignscore.sh"
    output:
        "result/{PREFIX}/evidence/pep/alignscore.txt"
    params:
        ref="genome.fa.masked.fa",
        dir="result/{PREFIX}/evidence/pep"
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
        QRY="result/{PREFIX}/evidence/pep/{lane_number}.fa",
        REF="result/{PREFIX}/evidence/genome.fa.masked.fa",
        alignscore="result/{PREFIX}/evidence/pep/alignscore.txt"
    output:
        gff="result/{PREFIX}/evidence/pep/{lane_number}.fa.gff",
        pro="result/{PREFIX}/evidence/pep/{lane_number}.fa.pro"
    params:
        dir="result/{PREFIX}/evidence/pep",
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
        "result/{PREFIX}/evidence/pep/{lane_number}.fa.slim.genblast.gff"
    shell:
        "python -m iga.annotation.genblast2 filter_genblast {input} > {output}"

rule filter_early_stop:
    input:
        rules.genblast.output.pro
    output:
        "result/{PREFIX}/evidence/pep/{lane_number}.fa.genblast.noearly_stop.id"
    shell:
        "python -m iga.annotation.genblast2 filter_early_stop {input} > {output}"

rule selectGFF:
    input:
        stopid=rules.filter_early_stop.output,
        slimgff=rules.filter_genblast.output
    output:
        "result/{PREFIX}/evidence/pep/{lane_number}.fa.filter.genblast.gff"
    shell:
        "bin/selectGFF.pl {input.stopid} {input.slimgff} > {output}"

rule sed_gff:
    input:
        rules.selectGFF.output
    output:
        "result/{PREFIX}/evidence/pep/{lane_number}.fa.final.gff"
    shell:
        "sed 's/transcript/protein_match/; s/coding_exon/match_part/' {input} > {output}"

def get_genblast_gff(wildcards):
    lane_dir = checkpoints.split_pep.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"result/{wildcards.PREFIX}/evidence/sample/{{lane_number}}.fa").lane_number
    gff = expand(rules.sed_gff.output, **wildcards, lane_number=lane_numbers)
    return gff

rule merge_gff:
    input:
        get_genblast_gff
    output:
        "result/{PREFIX}/evidence/genblast.gff"
    shell:
        "cat {input} > {output}"

