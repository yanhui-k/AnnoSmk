rule BuildDatabase:
    input:
        "AnnoSmk/{PREFIX}/{PREFIX}.fa"
    output:
        "{PREFIX}.nsq",
        "{PREFIX}.nog",
        "{PREFIX}.nni",
        "{PREFIX}.nnd",
        "{PREFIX}.nhr",
        "{PREFIX}.translation"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_BuildDatabase.tsv"
    shell:
        '''
        BuildDatabase -name {wildcards.PREFIX} -engine ncbi AnnoSmk/{wildcards.PREFIX}/{wildcards.PREFIX}.fa
        '''

rule RepeatModeler:
    input:
        fa="AnnoSmk/{PREFIX}/{PREFIX}.fa",
        nsq="{PREFIX}.nsq",
        nog="{PREFIX}.nog"
    output:
        famfa="AnnoSmk/{PREFIX}/evidence/{PREFIX}-families.fa",
        famstk="AnnoSmk/{PREFIX}/evidence/{PREFIX}-families.stk"
    params:
        db="{PREFIX}"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_RepeatModeler.tsv"
    threads: THREADS
    shell:
        '''
        RepeatModeler -pa {threads} -engine ncbi -database {params.db}
        cp {wildcards.PREFIX}-families.fa AnnoSmk/{wildcards.PREFIX}/evidence/
        cp {wildcards.PREFIX}-families.stk AnnoSmk/{wildcards.PREFIX}/evidence/
        '''

rule RepeatMasker:
    input:
        fa="AnnoSmk/{PREFIX}/{PREFIX}.fa",
        famfa=rules.RepeatModeler.output.famfa
    output:
        out="AnnoSmk/{PREFIX}/evidence/{PREFIX}.out"
    params:
        dir="AnnoSmk/{PREFIX}/evidence/"
    threads: THREADS
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_RepeatMasker.tsv"
    shell:
        '''
        RepeatMasker -lib {input.famfa} {input.fa} -pa {threads} -dir {params.dir}
        mv AnnoSmk/{wildcards.PREFIX}/evidence/*.out AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.PREFIX}.out
        mv AnnoSmk/{wildcards.PREFIX}/evidence/*.masked AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.PREFIX}.masked
        mv AnnoSmk/{wildcards.PREFIX}/evidence/*.tbl AnnoSmk/{wildcards.PREFIX}/evidence/{wildcards.PREFIX}.tbl
        '''

#rule post_repeatmasker_gunzip:
#    input:
#        catgz="AnnoSmk/{PREFIX}/evidence/{REF}.cat.gz"
#    output:
#        cat="AnnoSmk/{PREFIX}/evidence/{REF}.cat"
#    shell:
#        "gunzip {input.catgz} "
#
#rule ProcessRepeats:
#    input:
#        cat=rules.post_repeatmasker_gunzip.output.cat
#    output:
#        out="AnnoSmk/{PREFIX}/evidence/{REF}.out",
#        tbl="AnnoSmk/{PREFIX}/evidence/{REF}.tbl"
#    shell:
#        "ProcessRepeats -species Viridiplantae {input.cat}"

rule rmOutToGFF3:
    input:
        expand("AnnoSmk/{PREFIX}/evidence/{PREFIX}.out", PREFIX=PREFIX)
    output:
        "AnnoSmk/{PREFIX}/evidence/full_mask.gff3"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_rmOutToGFF3.tsv"
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
        "AnnoSmk/{PREFIX}/evidence/full_mask.complex.gff3"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_isolate_complex_repeats.tsv"
    shell:
        'grep -v -e "Satellite" -e ")n" -e "-rich" {input} > {output}'

rule reformat_to_work_with_MAKER:
    input:
        rules.isolate_complex_repeats.output
    output:
        "AnnoSmk/{PREFIX}/evidence/repeat.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_reformat_to_work_with_MAKER.tsv"
    shell:
        "reformat_to_work_with_MAKER.py {input} {output}"

rule bedtools1:
    input:
        fa="AnnoSmk/{PREFIX}/{PREFIX}.fa",
        gff3=rules.reformat_to_work_with_MAKER.output
    output:
        masked_fa="AnnoSmk/{PREFIX}/evidence/genome.fa.masked.fa"
    log:
        "AnnoSmk/{PREFIX}/evidence/bedtools.log"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_bedtools1.tsv"
    shell:
        '''
        bedtools maskfasta -soft -fi {input.fa} -bed {input.gff3} -fo {output.masked_fa}
        '''