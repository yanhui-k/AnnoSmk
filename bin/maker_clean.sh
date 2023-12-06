#format_maker_gff_to_parseval_gff.py mpi_maker.gff mpi_maker_parseval.gff
#
##get exon line from maker result
#python getExon.py mpi_maker.gff cluster0.gff
touch cluster0_sort_overlap_uniq.txt

#cluster
bedtools cluster -i cluster0.gff -s > cluster1.txt
sort -k1,1 -k9,9 -k4,4n cluster1.txt > cluster1_sort.txt
#get potentially overlapping exons
awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster1_sort.txt > cluster1_sort_overlap.txt
sort -k1,1 -k9,9 -k4,4n cluster1_sort_overlap.txt | uniq > cluster1_sort_overlap_uniq.txt

i=0
j=1
until [[ $(awk 'END {print NR}' cluster${i}_sort_overlap_uniq.txt) -eq $(awk 'END {print NR}' cluster${j}_sort_overlap_uniq.txt) ]]
do
    i=`expr $i + 1`
    j=`expr $j + 1`
    cut -f 1-9 cluster${i}_sort_overlap_uniq.txt > cluster${j}.gff
    bedtools cluster -i cluster${j}.gff -s > cluster${j}.txt
    sort -k1,1 -k9,9 -k4,4n cluster${j}.txt > cluster${j}_sort.txt
    awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster${j}_sort.txt > cluster${j}_sort_overlap.txt
    sort -k1,1 -k9,9 -k4,4n cluster${j}_sort_overlap.txt | uniq > cluster${j}_sort_overlap_uniq.txt
    cut -f 9 cluster${j}_sort_overlap_uniq.txt | uniq > overlap_exon.txt
done


##cluster again
#cut -f 1-9 cluster1_sort_overlap_uniq.txt > cluster2.gff
#bedtools cluster -i cluster2.gff -s > cluster2.txt
#sort -k1,1 -k9,9 -k4,4n cluster2.txt > cluster2_sort.txt
#awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster2_sort.txt > cluster2_sort_overlap.txt
#sort -k1,1 -k9,9 -k4,4n cluster2_sort_overlap.txt | uniq > cluster2_sort_overlap_uniq.txt
#
##cluster again
#cut -f 1-9 cluster2_sort_overlap_uniq.txt > cluster3.gff
#bedtools cluster -i cluster3.gff -s > cluster3.txt
#sort -k1,1 -k9,9 -k4,4n cluster3.txt > cluster3_sort.txt
#awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster3_sort.txt > cluster3_sort_overlap.txt
#sort -k1,1 -k9,9 -k4,4n cluster3_sort_overlap.txt | uniq > cluster3_sort_overlap_uniq.txt
#
##cluster again
#cut -f 1-9 cluster3_sort_overlap_uniq.txt > cluster4.gff
#bedtools cluster -i cluster4.gff -s > cluster4.txt
#sort -k1,1 -k9,9 -k4,4n cluster4.txt > cluster4_sort.txt
#awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster4_sort.txt > cluster4_sort_overlap.txt
#sort -k1,1 -k9,9 -k4,4n cluster4_sort_overlap.txt | uniq > cluster4_sort_overlap_uniq.txt
#
##cluster again
#cut -f 1-9 cluster4_sort_overlap_uniq.txt > cluster5.gff
#bedtools cluster -i cluster5.gff -s > cluster5.txt
#sort -k1,1 -k9,9 -k4,4n cluster5.txt > cluster5_sort.txt
#awk -v OFS="\t" '{if(a1==$10 && a2==$9){print a3"\n"$0;a1=$10;a2=$9;a3=$0} else {a1=$10;a2=$9;a3=$0}}' cluster5_sort.txt > cluster5_sort_overlap.txt
#sort -k1,1 -k9,9 -k4,4n cluster5_sort_overlap.txt | uniq > cluster5_sort_overlap_uniq.txt
#
##clean
#cut -f 9 cluster5_sort_overlap_uniq.txt | uniq > overlap_exon.txt
awk -v FS="-" '{print $1"-"$2"-"$3"-"$4"-"$5}' overlap_exon.txt > overlap_gene_pre.txt
sort overlap_gene_pre.txt | uniq > overlap_gene.txt
#grep -vf overlap_gene.txt mpi_maker_parseval.gff > mpi_maker_parseval_filter.gff

#test
#singularity exec -e -B `pwd` /nfs/yanhui/bin/aegean.simg parseval ../../TAIR/tair_refr_parseval.gff mpi_maker_parseval_filter.gff -s -w -o tair-out.txt