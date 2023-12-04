rule pre_pep:
    output:
        expand("{PEP}",PEP=PEP)
    benchmark:
        "AnnoSmk/benchmark/pre_pep.tsv"
    shell:
        '''
        touch {output}
        '''

rule prepep:
    input:
        pep=expand("{PEP}",PEP=PEP)
    output:
        "AnnoSmk/{PREFIX}/evidence/0.fa"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_prepep.tsv"
    shell:
        "modfasta.py {input.pep} {output}"

checkpoint split_pep:
    input:
        pep="AnnoSmk/{PREFIX}/evidence/0.fa"
    output:
        lane_dir=directory("AnnoSmk/{PREFIX}/evidence/sample/")
    params:
        size=600000
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_split_pep.tsv"
    shell:
        "split_fasta.py {input.pep} {output.lane_dir} {params.size}"

rule mv:
    input:
        fasta="AnnoSmk/{PREFIX}/evidence/sample/{lane_number}.fa"
    output:
        fa="AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_mv_{lane_number}.tsv"
    shell:
        "cp {input.fasta} {output.fa}"

rule prep_genblast:
    input:
        "AnnoSmk/{PREFIX}/evidence/genome.fa.masked.fa"
    output:
        "AnnoSmk/{PREFIX}/evidence/pep/alignscore.txt"
    params:
        ref="genome.fa.masked.fa",
        dir="AnnoSmk/{PREFIX}/evidence/pep"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_prep_genblast.tsv"
    shell:
        '''
        alignscore.sh > {output} 
        cd {params.dir}
        ln -sf `which formatdb` 
        ln -sf `which blastall`
        ln -sf ../{params.ref}
        cd -
        '''

rule genblast:
    input:
        QRY="AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa",
        REF="AnnoSmk/{PREFIX}/evidence/genome.fa.masked.fa",
        alignscore="AnnoSmk/{PREFIX}/evidence/pep/alignscore.txt"
    output:
        gff="AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.gff",
        pro="AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.pro"
    params:
        dir="AnnoSmk/{PREFIX}/evidence/pep",
        ref="genome.fa.masked.fa",
        qry="{lane_number}.fa"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_genblast_{lane_number}.tsv"
    shell:
        """
        cd {params.dir}
        awk '{{if ($0~/>/) name=$0; else seq[name]=seq[name]$0}} END {{for (i in seq) {{print i; print seq[i]}}}}' {params.qry} > {params.qry}_mod.fa
        genblastG -q {params.qry}_mod.fa -t {params.ref} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
            -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {params.qry}_mod.fa > {params.qry}_mod.fa.log
        wait
        mv {params.qry}_mod.fa_*.gff {params.qry}.gff
        mv {params.qry}_mod.fa_*.pro {params.qry}.pro
        if [ `grep -c 'stop codon not found between two HSPs that were originally separated by stop?' {params.qry}_mod.fa.log` -ne 0 ]; then
            rm {params.qry}_mod.fa*.blast
            rm {params.qry}_mod.fa*.blast.report
            seq=`grep 'current gene:' 1.fa.log | cut -d ' ' -f 3 | tail -n1`
            echo $seq
            var=`grep -n $seq {params.qry}_mod.fa`
            echo $var
            var1=`echo $var | cut -d ':' -f 1| awk '{{print int($0)+2}}'`
            var2=`echo $var | cut -d ':' -f 1| awk '{{print int($0)+3}}'`
            echo $var1
            echo $var2
            cp {params.qry}_mod.fa {params.qry}_mod.fa1
            sed -i '$var1,$var2d' {params.qry}_mod.fa1
            genblastG -q {params.qry}_mod.fa1 -t {params.ref} -e 1e-4 -g T -f F -a 0.5 -d 100000 -r 3 -c 0.5 -s 0 -i 15 \
            -x 20 -n 20 -v 2 -h 2 -j 0 -norepair -gff -cdna -pro -o {params.qry}_mod.fa1 > {params.qry}_mod.fa1.log
            mv -f {params.qry}_mod.fa1_*.gff {params.qry}.gff  
            mv -f {params.qry}_mod.fa1_*.pro {params.qry}.pro 
        fi                
        /bin/rm {params.qry}*.blast
        /bin/rm {params.qry}*.blast.report
        /bin/rm {params.qry}_*_0
        /bin/rm {params.qry}_*_0.DNA  
        cd -
        """

rule filter_genblast:
    input:
        rules.genblast.output.gff
    output:
        "AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.slim.genblast.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_filter_genblast_{lane_number}.tsv"
    shell:
        "genblast2.py filter_genblast {input} > {output}"

rule filter_early_stop:
    input:
        rules.genblast.output.pro
    output:
        "AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.genblast.noearly_stop.id"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_filter_early_stop_{lane_number}.tsv"
    shell:
        "genblast2.py filter_early_stop {input} > {output}"

rule selectGFF:
    input:
        stopid=rules.filter_early_stop.output,
        slimgff=rules.filter_genblast.output
    output:
        "AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.filter.genblast.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_selectGFF_{lane_number}.tsv"
    shell:
        "selectGFF.pl {input.stopid} {input.slimgff} > {output}"

rule sed_gff:
    input:
        rules.selectGFF.output
    output:
        "AnnoSmk/{PREFIX}/evidence/pep/{lane_number}.fa.final.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_sed_gff_{lane_number}.tsv"
    shell:
        "sed 's/transcript/protein_match/; s/coding_exon/match_part/' {input} > {output}"

def get_genblast_gff(wildcards):
    lane_dir = checkpoints.split_pep.get(**wildcards).output[0]
    lane_numbers = glob_wildcards(f"AnnoSmk/{wildcards.PREFIX}/evidence/sample/{{lane_number}}.fa").lane_number
    gff = expand(rules.sed_gff.output, **wildcards, lane_number=lane_numbers)
    return gff

rule merge_genblast_gff:
    input:
        get_genblast_gff
    output:
        "AnnoSmk/{PREFIX}/evidence/genblast.gff"
    benchmark:
        "AnnoSmk/benchmark/{PREFIX}_merge_genblast_gff.tsv"
    shell:
        '''
        cat {input} >> {output}
        /bin/rm -rf AnnoSmk/{wildcards.PREFIX}/evidence/pep
        '''