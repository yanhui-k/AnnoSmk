import os
import glob

samples_trans=glob.glob(r"altra.*_1.clean.fq.gz")
samples_trans=[x.replace("_1.clean.fq.gz","") for x in samples_trans]
samples_flnc=glob.glob(r"workdir_isoseq_altra/*hq.fasta.gz")
samples_rnaseq=glob.glob(r"*trinity/Trinity-GG.fasta")
samples_minimap=["flnc","rnaseq"]
number = [1,2,3,4,5,6,7,8]
number_to_split_pep=10
PREFIX=range(1,number_to_split_pep+1)


rule all:
    input:
        "test.altra.genome.fa.masked.fa",
        "flnc.fasta",
        "rnaseq.fasta",
        "total_est.gff",
        "genblast.gff",
        "flnc.fasta.rawgff.gff",


rule BuildDatabase:
    input:
        fa="altra.genome.fa"
    output:
        "altra.genome.nsq","altra.genome.nog",
        "altra.genome.nni","altra.genome.nnd",
        "altra.genome.nhr","altra.genome.translation"
    shell:
        "BuildDatabase -name altra -engine ncbi {input.fa} "

rule RepeatModeler:
    input:
        fa="altra.genome.fa",
        nsq="altra.genome.nsq",
        nog="altra.genome.nog",
    output:
        famfa="altra.genome.fa-families.fa",
        famstk="altra.genome.fa-families.stk"
    shell:
        "RepeatModeler -engine ncbi -pa 100 -database {input.fa}"

rule RepeatMasker:
    input:
        fa="altra.genome.fa",
        famfa=rules.RepeatModeler.output.famfa
    output:
        dir="custom_lib.out/altra.genome.fa.cat.gz"
    shell:
        "RepeatMasker -lib {input.famfa} {input.fa}-pa 30 -dir {output.dir}"

rule post_repeatmasker_gunzip:
    input:
        catgz="custom_lib.out/altra.genome.fa.cat.gz"
    output:
        cat="custom_lib.out/altra.genome.fa.cat"
    shell:
        "gunzip {input.catgz} "

rule ProcessRepeats:
    input:
        cat=rules.post_repeatmasker_gunzip.output.cat
    output:
        out="custom_lib.out/altra.genome.fa.out",
        tbl="custom_lib.out/altra.genome.fa.tbl"
    shell:
        "ProcessRepeats -species $species {input.cat}"

rule rmOutToGFF3:
    input:
        rules.ProcessRepeats.output.out
    output:
        "full_mask.gff3"
    shell:
        "rmOutToGFF3.pl {input} > {output}"

rule isolate_complex_repeats:
    input:
        rules.rmOutToGFF3.output
    output:
        "full_mask.complex.gff3"
    shell:
        'grep -v -e "Satellite" -e ")n" -e "-rich" {input} > {output}'

rule reformat_to_work_with_MAKER:
    input:
        rules.isolate_complex_repeats.output
    output:
        "full_mask.complex.reformat.gff3"
    shell:
        '''
        cat {input} | \
            perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > {output}
        '''

rule bedtools1:
    input:
        fa="altra.genome.fa",
        gff3=rules.reformat_to_work_with_MAKER.output
    output:
        "test.altra.genome.fa.masked.fa"
    log:
        "logs/bedtools.log"
    shell:
        "bedtools maskfasta -soft -fi {input.fa} -bed {input.gff3} -fo {output}"

rule hisat2build:
    input:
        fa = "altra.genome.fa"
    output:
        expand("altra.genome.fa.{sample}.ht2",sample=number)
    shell:
        "hisat2-build {input} {input}"

rule hisat2:
    input:
        fa = "altra.genome.fa",
        r1 = "{sample}_1.clean.fq.gz",
        r2 = "{sample}_2.clean.fq.gz",
        index=rules.hisat2build.output
    output:
        "{sample}.clean.fq.gz.bam",
    shell:
        "hisat2 --rna-strandness RF --mp 3,1 -p 30 -x {input.fa} -1 {input.r1} -2 {input.r2} | samtools sort -@ 30 -o {output}"

rule trinity:
    input:
        "{sample}.clean.fq.gz.bam"
    output:
        dir="{sample}.clean.fq.gz_trinity",
        ggfa="{sample}.clean.fq.gz_trinity/Trinity-GG.fasta",
    shell:
        "Trinity --SS_lib_type RF \
          --genome_guided_bam {input} \
          --genome_guided_max_intron 15000 \
          --max_memory 20G --CPU 30 \
          --output {output.dir}"

rule cat_est_flnc:
    input:
        expand("{sample}",sample=samples_flnc)
    output:
        "flnc.fasta"
    shell:
        "python -m iga.annotation.maker cat_est {input} > {output}"

rule cat_est_rnaseq:
    input:
        expand("{sample}.clean.fq.gz_trinity/Trinity-GG.fasta",sample=samples_trans)
    output:
        "rnaseq.fasta"
    shell:
        "python -m iga.annotation.maker cat_est {input} > {output}"

rule minimap2:
    input:
        fa1 = "altra.genome.fa",
        fa2 = expand("{sample}.fasta",sample=samples_minimap)
    output:
        expand("{sample}.fasta.bam",sample=samples_minimap)
    threads:
        20
    shell:
        "minimap2 -t{threads} -C5 -ax splice {input.fa1} {input.fa2} |samtools view -F 256 -b > {output}"

rule bedtools:
    input:
        expand("{sample}.fasta.bam",sample=samples_minimap)
    output:
        expand("{sample}.fasta.raw.bed",sample=samples_minimap)
    shell:
        "bedtools bamtobed -split -i {input} > {output}"

rule awk:
    input:
        expand("{sample}.fasta.raw.bed",sample=samples_minimap)
    output:
        expand("{sample}.fasta.bed",sample=samples_minimap)
    shell:
        "awk '$3-$2 > 1' {input} > {output}"

rule bed2gff3:
    input:
        expand("{sample}.fasta.bed",sample=samples_minimap)
    output:
        expand("{sample}.fasta.rawgff",sample=samples_minimap)
    shell:
        "gt bed_to_gff3 {input} | sort -k9,9 -k1,1 -k7,7 -k4,4n  > {output}"

rule format_gt_gff_to_maker_gff:
    input:
        rules.bed2gff3.output
    output:
        expand("{sample}.fasta.rawgff.gff",sample=samples_minimap)
    shell:
        "python -m iga.annotation.maker format_gt_gff_to_maker_gff {input}"

rule cat_gff:
    input:
        rules.format_gt_gff_to_maker_gff.output
    output:
        "total_est.gff"
    shell:
        "cat {input} > {output}"

rule split_fastav3:
    input:
        pep="arath_med_sprot.pep"
    output:
        "arath_med_sprot.pep._/{sample}.fa"
    params:
        number=100
    shell:
        "split_fastav3.pl {input.pep} {params.number}"

rule prep_genblast:
    input:
        "alignscore.sh"
    output:
        "alignscore.txt"
    shell:
        "sh {input} > {output} && ln -s `which formatdb` && ln -s `which blastall`"

rule genblast:
    input:
        QRY="arath_med_sprot.pep._/{sample}.fa",
        REF="altra.genome.fa.masked.fa",
        alignscore="alignscore.txt"
    output:
        "arath_med_sprot.pep._/{sample}.fa.genblast_1.1c_2.3_s2_tdshift2_tddis0_tcls0.0_m2_score_i0_d16_0.gff",
        "arath_med_sprot.pep._/{sample}.fa.genblast_1.1c_2.3_s2_tdshift2_tddis0_tcls0.0_m2_score_i0_d16_0.pro"
    shell:
        """
        genblastG -p genblastg -q {input.QRY} -t {input.REF} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
            -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {output}
        """

rule filter_genblast:
    input:
        "arath_med_sprot.pep._/{sample}.fa.genblast_1.1c_2.3_s2_tdshift2_tddis0_tcls0.0_m2_score_i0_d16_0.gff"
    output:
        "arath_med_sprot.pep._/{sample}.fa.slim.genblast.gff"
    shell:
        "python -m iga.annotation.genblast2 filter_genblast {input} > {output}"

rule filter_early_stop:
    input:
        "arath_med_sprot.pep._/{sample}.fa.genblast_1.1c_2.3_s2_tdshift2_tddis0_tcls0.0_m2_score_i0_d16_0.pro"
    output:
        "arath_med_sprot.pep._/{sample}.fa.genblast.noearly_stop.id"
    shell:
        "python -m iga.annotation.genblast2 filter_early_stop {input} > {output}"

rule selectGFF:
    input:
        stopid=rules.filter_early_stop.output,
        slimgff=rules.filter_genblast.output
    output:
        "arath_med_sprot.pep._/{sample}.fa.filter.genblast.gff"
    shell:
        "selectGFF.pl {input.stopid} {input.slimgff} > {output}"

rule sed_gff:
    input:
        rules.selectGFF.output
    output:
        "arath_med_sprot.pep._/{sample}.fa.final.gff"
    shell:
        "sed 's/transcript/protein_match/; s/coding_exon/match_part/' {input} > {output}"

rule merge_gff:
    input:
        expand(rules.sed_gff.output,sample=PREFIX)
    output:
        "genblast.gff"
    shell:
        "cat {input} > {output}"
