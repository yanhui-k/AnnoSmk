"""
rna-seq relevant utils
"""
from collections import OrderedDict, defaultdict

from iga.apps.base import emain, get_prefix, logger, rscript
#import numpy as np


def remove_dup_table(table=None):
    tpm_dict = defaultdict(float)
    with open(table) as fh:
        for line in fh:
            mylist = line.rstrip().split('\t')
            try:
                tpm_float = float(mylist[-1])
            except ValueError:
                logger.warning(line)
                continue
            tpm_dict[mylist[0]] += tpm_float
    new_table = table + ".uniq"
    with open(new_table, 'w') as fh:
        fh.write("Gene ID\tTPM\n")
        for k in tpm_dict:
            fh.write("{}\t{}\n".format(k, tpm_dict[k]))
    return new_table

def merge_exp_table(tables=None):
    """
    Merge expression table produced by stringtie into one with TPM
    :param tables:
    :return:
    """
    import pandas as pd
    # sample_list = []
    # gene_dict = []
    z = pd.DataFrame()
    # 'Gene Name', 'Reference', 'Strand', 'Start', 'End',
    sub_list = ['Gene ID', 'TPM']
    for i, t in enumerate(tables):
        new_table = remove_dup_table(t)
        this_df = pd.read_table(new_table, sep='\t')
        logger.warning(this_df.columns)
        sub_df = this_df[list(sub_list)]
#        sub_df = sub_df.groupby('Gene ID').TPM.apply(lambda g: g.nlargest(2).sum())
        #Now change subdf's column name
        t_prefix = get_prefix(t)
        new_sublist = sub_list.copy()
        new_sublist[new_sublist.index('TPM')] = t_prefix
        sub_df.columns = new_sublist
        if i == 0:
            z = sub_df
        else:
            z = z.merge(sub_df, left_on=sub_list[0],
                    right_on=sub_list[0], how='outer')
        #sample_list.append(get_prefix(t))
    #z.to_csv('merged_expression1.txt', sep="\t", index=False, na_rep='0.0')
    prev_column = list(z.columns)
    prev_column[0] = ""
    z.columns = prev_column
    z = z.sort_index()
    z.to_csv('merged_expression.txt', sep="\t", index=False, na_rep='0.0')
    return z



pheatmap_sh = r"""
RCG<-read.table("{0}", header = T, row.names = 1, sep="\t")

library("pheatmap")
library("RColorBrewer")

# mycol =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)

RCGm <- as.data.frame(RCG)
# 
RCG <- ifelse(RCGm==0, 0.01, RCG)
# 
# col_size = length(colnames(RCG))
# 
# RCG <- RCG[rowSums(RCG[,c(1:col_size)]) > 5]
# 
# RCG <- log2(RCG/rowMeans(RCG))

RCG <- log2(RCG)

prefix="{0}"

pdf(paste(prefix, "pdf", sep='.'))

a = pheatmap(RCG, show_rownames=F, 
       main = "{1}", cluster_rows=T, cluster_cols=F,
       fontsize_col = 20, angle_col ="45",border_color = 'white')
"""


def pheatmap(table=None, main=''):
    if main == '':
        main = table
    main = "Expression Heatmap"
    cmd = pheatmap_sh.format(table, main)
    rscript(cmd)
    return 0


if __name__ == "__main__":
    emain()
