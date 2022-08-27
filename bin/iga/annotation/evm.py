""" EVM utile"""

evm_sh = r"""
# $head *gff
# ==> est.gff <==
# 000000F	StringTie	EST_match	11950	12646	1000	-	.	ID=MSTRG.1.1:exon-1;Name=MSTRG.1;Parent=MSTRG.1.1
# 000000F	StringTie	EST_match	12826	13501	1000	-	.	ID=MSTRG.1.1:exon-2;Name=MSTRG.1;Parent=MSTRG.1.1
# 000000F	StringTie	EST_match	14834	15508	1000	+	.	ID=MSTRG.2.1:exon-1;Name=MSTRG.2;Parent=MSTRG.2.1
# 000000F	StringTie	EST_match	16435	17476	1000	+	.	ID=MSTRG.2.1:exon-2;Name=MSTRG.2;Parent=MSTRG.2.1
# 000000F	StringTie	EST_match	17881	18718	1000	+	.	ID=MSTRG.2.1:exon-3;Name=MSTRG.2;Parent=MSTRG.2.1
# 000000F	StringTie	EST_match	19250	19564	1000	+	.	ID=MSTRG.3.1:exon-1;Name=MSTRG.3;Parent=MSTRG.3.1
# 000000F	StringTie	EST_match	19818	23918	1000	+	.	ID=MSTRG.3.1:exon-2;Name=MSTRG.3;Parent=MSTRG.3.1
# 000000F	StringTie	EST_match	24056	24985	1000	+	.	ID=MSTRG.3.1:exon-3;Name=MSTRG.3;Parent=MSTRG.3.1
# 000000F	StringTie	EST_match	19250	19564	1000	+	.	ID=MSTRG.3.2:exon-1;Name=MSTRG.3;Parent=MSTRG.3.2
# 000000F	StringTie	EST_match	19753	20014	1000	+	.	ID=MSTRG.3.2:exon-2;Name=MSTRG.3;Parent=MSTRG.3.2
# 
# ==> gene.gff <==
# 000042F	maker	gene	108408	110326	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0;Name=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0;score=1000
# 000042F	maker	mRNA	108408	110326	1137	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0;Name=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1;_AED=0.03;_eAED=0.03;_QI=141|1|1|1|0|0|4|0|331
# 000042F	maker	exon	108408	109087	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:1;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	exon	109221	109639	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:2;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	exon	110289	110326	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:3;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	five_prime_UTR	108408	108548	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:five_prime_utr;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	CDS	108549	109087	.	+	0	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:cds;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	CDS	109221	109639	.	+	1	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:cds;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	CDS	110289	110326	.	+	2	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1:cds;Parent=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.0-mRNA-1
# 000042F	maker	gene	110337	111540	.	+	.	ID=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.1;Name=maker-000042F|arrow_np1212-est_gff_StringTie-gene-1.1;score=1000
# 
# ==> pep.gff <==
# 000030F	genBlastG	nucleotide_to_protein_match	2650826	2652085	.	-	.	ID=AT3G18710.1R1-1-A1-E1;Parent=AT3G18710.1R1-1-A1
# 000006F	genBlastG	nucleotide_to_protein_match	2488347	2489606	.	+	.	ID=AT3G18710.1R2-2-A1-E1;Parent=AT3G18710.1R2-2-A1
# 000015F	genBlastG	nucleotide_to_protein_match	4391096	4392352	.	+	.	ID=AT3G18710.1R3-3-A1-E1;Parent=AT3G18710.1R3-3-A1
# 000062F	genBlastG	nucleotide_to_protein_match	85260	85504	.	+	.	ID=AT3G18710.1R4-4-A1-E1;Parent=AT3G18710.1R4-4-A1
# 000062F	genBlastG	nucleotide_to_protein_match	85598	85727	.	+	.	ID=AT3G18710.1R4-4-A1-E2;Parent=AT3G18710.1R4-4-A1
# 000062F	genBlastG	nucleotide_to_protein_match	85767	86522	.	+	.	ID=AT3G18710.1R4-4-A1-E3;Parent=AT3G18710.1R4-4-A1
# 000001F	genBlastG	nucleotide_to_protein_match	8917925	8919107	.	-	.	ID=AT3G18710.1R5-5-A1-E1;Parent=AT3G18710.1R5-5-A1
# 000001F	genBlastG	nucleotide_to_protein_match	8917794	8917867	.	-	.	ID=AT3G18710.1R5-5-A1-E2;Parent=AT3G18710.1R5-5-A1
# 000009F	genBlastG	nucleotide_to_protein_match	4184709	4184907	.	+	.	ID=sp_D5JWB3_SARED_ESCCAR1-1-A1-E1;Parent=sp_D5JWB3_SARED_ESCCAR1-1-A1
# 000009F	genBlastG	nucleotide_to_protein_match	4184926	4184972	.	+	.	ID=sp_D5JWB3_SARED_ESCCAR1-1-A1-E2;Parent=sp_D5JWB3_SARED_ESCCAR1-1-A1
#!/bin/bash
#BSUB -J evm_partition      # job name
#BSUB -n 1                   # number of tasks in job
#BSUB -q Q104C512G_X4              # queue
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

set -euxo pipefail

ROOT=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/input/genblast/falcon_genblast

#touch null.gff
#awk '$3=="coding_exon"' total.gff  | sed 's/coding_exon/nucleotide_to_protein_match/; s/-R/R/g' > nucleotide_to_protein_match.gff
#sort -k1,1 -k4,4n -k7,7 nucleotide_to_protein_match.gff > nucleotide_to_protein_match.sort.gff
#bsub -q Q104C512G_X4  -o output.%J -e error.%J "evidence_modeler.pl --genome falcon_v340_sgs_polish.fa --gene_predictions null.gff --protein_alignments n
ucleotide_to_protein_match.sort.gff --weights weights.txt >tmp.out 2>tmp.err"

#bsub -q Q104C512G_X4  -o output.%J -e error.%J "evidence_modeler.pl --genome test_evm.fa --gene_predictions  test_gene_format.gff --protein_alignments te
st_pep.gff --weights weights.txt >tmp.out 2>tmp.err"

#Prepare
#cat /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/rnaseq/liuhui/input_corfalcon/corfalcon.gff |grep "match_part" | sed 's/match_part/EST_match/' >
est.gff
#cat /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/rnaseq/liuhui/input_corfalcon/flcdna.gff |sed 's/StringTie/FLcdna/' | grep "match_part" | sed 's/
match_part/cDNA_match/' >> est.gff
#awk '$2=="maker"' gene.gff.link | sed 's/|arrow_np1212//' > gene.gff
#sed 's/|arrow_np1212//' est.gff.raw >est.gff
#sed 's/|arrow_np1212//' pep.gff.raw > pep.gff
#sed 's/|arrow_np1212//' falcon_v340_sgs_polish.fa.raw > falcon_v340_sgs_polish.fa
#exit

ROOT=$PWD
REF=falcon_v340_sgs_polish.fa
GENEGFF=gene.gff
PEPGFF=pep.gff
ESTGFF=est.gff

chunck_size=1000000
step_size=100000
#EVM pipeline
$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome ${REF} --gene_predictions ${GENEGFF} --protein_alignments ${PEPGFF} --transcript_alignments ${ESTGFF}
--segmentSize ${chunck_size} --overlapSize ${step_size} --partition_listing partitions_list.out


$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome ${REF} --weights `pwd`/weights.txt \
      --gene_predictions ${GENEGFF}  \
      --protein_alignments ${PEPGFF} \
      --transcript_alignments ${ESTGFF} \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list

#$EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

#$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
#
#$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ${REF}

cat */evm.out.gff3  > total_evm.out.gff3

#sed 's/|arrow_np1212//' falcon_v340_sgs_polish.fa.link > falcon_v340_sgs_polish.fa
#sed 's/|arrow_np1212//' total_format.gff > falcon_genblast_format.gff
"""