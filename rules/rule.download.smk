rule download:
    output:
        temp("{PREFIX}/{sample}.sra")
    benchmark:
        "AnnoSmk/benchmark/download_{PREFIX}_{sample}.tsv"
    shell:
        """
        wget -O {output} https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.sample}/{wildcards.sample}
        """

rule fastq_dump_d:
    input:
        "{PREFIX}/{sample}.sra"
    output:
        "{PREFIX}/{sample}_1.fastq",
        "{PREFIX}/{sample}_2.fastq"        
    benchmark:
        "AnnoSmk/benchmark/fastq_dump_{PREFIX}_{sample}.tsv"
    shell:
        '''
        fastq-dump --split-3 {input} -O {wildcards.PREFIX} 
        '''

rule fastq_dump_s:
    input:
        "{PREFIX}/{sample}.sra"
    output:
        "{PREFIX}/{sample}_s.fastq"        
    benchmark:
        "AnnoSmk/benchmark/fastq_dump_{PREFIX}_{sample}.tsv"
    shell:
        '''
        fastq-dump --split-3 {input} -O {wildcards.PREFIX} 
        mv {wildcards.PREFIX}/{wildcards.sample}.fastq {output}
        '''

#rule fastq_dump:
#    input:
#        "{PREFIX}/{sample}.sra"
#    output:
#        "{PREFIX}/{sample}.txt"        
#    benchmark:
#        "AnnoSmk/benchmark/fastq_dump_{PREFIX}_{sample}.tsv"
#    shell:
#        '''
#        fastq-dump --split-3 {input} -O {wildcards.PREFIX} &>> {output}
#        '''
#
#rule chech_download:
#    input:
#        expand("{PREFIX}/{sample}.txt",PREFIX=PREFIX, sample=config["DOWNLOAD_SRAID"])
#    output:
#        "AnnoSmk/download.log"
#    benchmark:
#        "AnnoSmk/benchmark/check_download.tsv"
#    shell:
#        '''
#        echo "\ndownload finish\n" > {output}
#        '''