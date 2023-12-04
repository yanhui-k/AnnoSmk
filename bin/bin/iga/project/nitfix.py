"""
scripts used in nitfix project
"""
import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


symbiosis_gene_fam_sh = r"""

RCG<-read.table("Orthogroups.GeneCount.symbiosis.tsv.count", header = T, row.names = 1, sep="\t")

library("pheatmap")
library("RColorBrewer")

# mycol =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)

RCGm <- as.matrix(RCG)

RCG <- ifelse(RCGm>5, 5, RCGm)

mycol = brewer.pal(n=6, name="PuBu")

#
# col_size = length(colnames(RCG))
#
# RCG <- RCG[rowSums(RCG[,c(1:col_size)]) > 5]
#
# RCG <- log2(RCG/rowMeans(RCG))

prefix="Orthogroups.GeneCount.symbiosis.tsv.count"

pdf(paste(prefix, "pdf", sep='.'), w=10, h=30)

a = pheatmap(RCG, show_rownames=T, col = mycol,
       main = "Orthogroups.GeneCount.symbiosis.tsv.count", cluster_rows=T, cluster_cols=F,
       fontsize_col = 20, angle_col ="45",border_color = 'white')

"""