localrules:
    split_by_Ns

rule split_by_Ns:
    input:
        fasta=expand("{REF}",REF=REF)
    output:
        pre_fa="AnnoSmk/{PREFIX}/{PREFIX}_pre.fa",
        pre_fai="AnnoSmk/{PREFIX}/{PREFIX}_pre.fa.fai",
        pre_dict="AnnoSmk/{PREFIX}/{PREFIX}_pre.fa.dict",
        pre_log="AnnoSmk/{PREFIX}/{PREFIX}_split.log",
        split_gff="AnnoSmk/{PREFIX}/{PREFIX}_split.gff",
        split_bed="AnnoSmk/{PREFIX}/{PREFIX}_split.bed",        
        split_fa="AnnoSmk/{PREFIX}/{PREFIX}_split.fa",
        split_txt="AnnoSmk/{PREFIX}/{PREFIX}_split.txt"
    params:
        minNs=5000
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_split_by_Ns.tsv"
    shell:
        """
        cp {input.fasta} {output.pre_fa}
        samtools faidx {output.pre_fa} --output {output.pre_fai}
        samtools dict {output.pre_fa} -o {output.pre_dict}
        picard ScatterIntervalsByNs -REFERENCE {output.pre_fa} -OUTPUT {output.split_gff} -MAX_TO_MERGE {params.minNs} -OUTPUT_TYPE ACGT &> {output.pre_log}
        mod_split_bed.py {output.split_gff} {output.split_bed}
        bedtools getfasta -fi {output.pre_fa} -bed {output.split_bed} -fo {output.split_fa}
        rename_gff_head.py {output.split_bed} {output.split_txt}
        """

rule maker_fasta:
    input:
        fa="AnnoSmk/{PREFIX}/{PREFIX}_split.fa"
    output:
        fasta="AnnoSmk/{PREFIX}/{PREFIX}.fa"
    shell:
        "samplify_fa_id.py {input.fa} {output.fasta}"


rule mod_gff:
    input:
        split_txt="AnnoSmk/{PREFIX}/{PREFIX}_split.txt",
        out_gff="AnnoSmk/{PREFIX}.gff"
    output:
        split_gff="AnnoSmk/{PREFIX}_fnl.gff"
    shell:
        """
        cp {input.out_gff} {output.split_gff}
        rename_gff.sh {input.split_txt} AnnoSmk/{PREFIX}_split.gff {output.split_gff}
        """