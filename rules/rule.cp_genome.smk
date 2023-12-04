rule cp_genome:
    input:
        fasta=expand("{REF}",REF=REF)
    output:
        fasta="AnnoSmk/{PREFIX}/{PREFIX}.fa"
    shell:
        '''
        cp {input.fasta} {output.fasta}
        '''

rule cp_gff:
    input:
        gff="AnnoSmk/{PREFIX}.gff"
    output:
        gff="AnnoSmk/{PREFIX}_fnl.gff"
    shell:
        """
        cp {input.gff} {output.gff}
        """