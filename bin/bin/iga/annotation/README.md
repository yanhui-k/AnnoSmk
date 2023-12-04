## Gene annotation with MAKER

### RepeatMasker

* Input
    - Genome
* Output
    - Repeat.gff (GFF format of repeat)
    - Repeat.tbl (Summary of repeat ratio)
    - Softmasked genome (complex repeat only, for maker and braker annotation)

```
# Step1. Run repeatmasker with both homology and de novo approach
python -m iga.annotation.repeat repeatmasker --species Viridiplantae elumb.contig.fa --threads 10
## Output:
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.cat.gz  
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.out
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.masked  
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.tbl

## The following step is slow  
python -m iga.annotation.repeat repeatmasker --denovo T elumb.contig.fa --threads 10
## Output:
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.cat.gz  
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.out
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.masked  
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.tbl

# Step2. Combine the two result 

python -m iga.annotation.repeat post_repeatmasker workdir_repeatmask_elumb.contig.fa  elumb.contig.fa

```

### Prepare RNA-seq evidence

```
cd reads/species
#$for i in *flcdna*bam; do python -m iga.annotation.isoseq isoseq_pb --threads 30 ${i} ${i}.IsoSeqPrimers.fasta; done
REF=sesep.genome.fa
mkdir -p annotation_evidence
# align and assemble RNA-Seq reads with Trinity. Make sure Hisat2 index is available
# This step may use 2h for each library
python -m iga.annotation.rnaseq reads_align_assembly  "sesep.ssRNA.CAS_Y5.leaf1_1.clean.fq.gz sesep.ssRNA.CAS_Y5.leaf1_2.clean.fq.gz"  $REF --threads 20
python -m iga.annotation.rnaseq reads_align_assembly  "sesep.ssRNA.CAS_Y6.leaf2_1.clean.fq.gz sesep.ssRNA.CAS_Y6.leaf2_2.clean.fq.gz"  $REF --threads 20
python -m iga.annotation.rnaseq reads_align_assembly  "sesep.ssRNA.CAS_Y7.leaf3_1.clean.fq.gz sesep.ssRNA.CAS_Y7.leaf3_2.clean.fq.gz"  $REF --threads 20

# waiting for job to finish
# Rename, and merge transcript.fasta in seperater folder into one single file. Local run for the following two lines
# This step is fasta
python -m iga.annotation.maker cat_est workdir_isoseq_*/*hq.fasta.gz > annotation_evidence/flnc_rna.fasta
python -m iga.annotation.maker cat_est *trinity/*GG.fasta > annotation_evidence/rnaseqGG.fasta 
cd annotation_evidence

# Submit to grid to execute. But this step is quite fasta, depending on the size of the est.fasta file.
python -m iga.annotation.maker fastq2gff flnc_rna.fasta ../${REF}
python -m iga.annotation.maker fastq2gff rnaseq.fasta ../${REF}


```


### Prepare Protein evidence

```
# One
python -m iga.annotation.maker prep_genblast elumb.contig.fa.masked.fa pep/elumb.homologous_pep.fa
* Input: 
    - rna_fasta.bam
    - protein_evidence.fasta
    - Genome (soft masked of complex region)
* OutputHMM: test_datisca/species/Sp_5/
```



### Maker
```


REF=elumb.contig.fa.masked.fa
FLNCESTGFF=flnc_rna.gff
ESTGFF=total_est.gff
CDNAFASTA=flnc_rna.fasta
PEPGFF=pep.gff
REPEATGFF=repeat.gff


#-+-First Round
ROUND=1
#python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${FLNCESTGFF} ${PEPGFF} ${REPEATGFF} --cpus 2
python -m iga.annotation.maker maker_check_resub ${REF}_R1
python -m iga.annotation.maker maker_collect     ${REF}_R1
python -m iga.annotation.maker maker_train       ${REF}_R1 --cdna_fasta ${CDNAFASTA}  --snap 'T' --augustus F
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Second Round
PREV=1
ROUND=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${FLNCESTGFF} ${PEPGFF} ${REPEATGFF} --round ${ROUND} --augustus_species ${REF}_R${PREV}_direct --snap_h
mm ${REF}_R${PREV} --update "alt_splice=1"
python -m iga.annotation.maker maker_check_resub ${REF}_R${ROUND}
python -m iga.annotation.maker maker_collect     ${REF}_R${ROUND}
python -m iga.annotation.maker maker_train       ${REF}_R${ROUND} --cdna_fasta ${CDNAFASTA} --augustus F
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Third Round #augustus is directly trained
ROUND=3
PREV=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${ESTGFF} ${PEPGFF} ${REPEATGFF} --round ${ROUND} --augustus_species ${REF}_R${PREV}_direct --snap_hmm $
{REF}_R${PREV} --update "trna=1;alt_splice=1"
#--queue Q64C1T_X4
python -m iga.annotation.maker maker_check_resub ${REF}_R${ROUND}
python -m iga.annotation.maker maker_collect     ${REF}_R${ROUND}
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..


# For Fourth Round, one can also use braker trained hmm. to see whether is better than round3
# Note one need to modify pep.fa and est.bam
python -m iga.annotation.maker braker --prefix elumb.contig.fa.masked.fa pep/elumb.homologous_pep.fa est/both_merge.bam
# Output trained hmm is located in /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/braker2/config/
# with species name is previously defined prefix

#-+-Fourth Round #augustus is directly trained
ROUND=4
PREV=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${ESTGFF} ${PEPGFF} ${REPEATGFF} --round ${ROUND} --augustus_species Sp_8 --snap_hmm ${REF}_R${P
#--queue Q64C1T_X4
python -m iga.annotation.maker maker_check_resub ${REF}_R${ROUND}
python -m iga.annotation.maker maker_collect     ${REF}_R${ROUND}
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

```

###  Refine, Rename and functional annotation

```
ROUND=3

# Genes will be renamed to ELUMB0001, transcripts will be renamed to ELUMB0001-t1
GenePrefix=ELUMB

cd ${{REF}}_R${{ROUND}}
python -m iga.annotation.maker pasa_refine ref.fa genome.maker.gff ../$CDNAFASTA
cp ref.fa.sqlite.gene_structures_post_PASA_updates.*.gff3 ref.fa.pasa.gff3
chmod -w ref.fa.sqlite.gene_structures_post_PASA_updates.*.gff3

bsub512 "python -m iga.annotation.maker maker_rename_gff ref.fa.pasa.gff3 --prefix ELUMB"
grep -i trna genome.maker.gff > trna.gff
cat ref.fa.pasa.gff3 trna.gff > ${{REF}}.gene_structure.gff3
gff_genome_to_genes.pl ${{REF}}.gene_structure.gff3 ref.fa > ${{REF}}.gene_structure.cds
cds2aa.pl ${{REF}}.gene_structure.cds > ${{REF}}.gene_structure.pep
python -m iga.assembly.assess busco --mode prot ${{REF}}.gene_structure.pep
python -m iga.annotation.maker func_anno pasa_raw.format.gff.pep
```