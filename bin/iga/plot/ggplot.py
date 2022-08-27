"""
ggplot wrappers
"""
from iga.apps.base import emain, rscript

import logging
import coloredlogs

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

theme_publication_r = r"""
theme_Publication <- function(base_size=14, base_family="sans") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))

}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506", 
      "#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
      "#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
"""
# 0 input_table
# 1 x=xx
# 2 y=xx
# 3 fill=xx
# 4 etc string like theme and others
barplot_r = r"""
library(ggplot2);
a=read.table("{0}", header=T, row.names=NULL); 
ggplot(a, aes(x={1}, y={2}, fill={3})) + geom_bar(stat="identity", position=position_dodge())  {4} 
ggsave("{0}.pdf", width = {5}, height = {6})
"""


def barplot(table=None, x='', y='', group='', theme='Publication', horizonal='F', pallette=""):
    r"""
    barplot with ggplot2
    :param table: Input table like
                Group   Number  Regulation
                Root-Nodule     1654    Down
                Root-Nodule     1496    Up
                Root-Leaf       4099    Down
    :param x: default is the 1st colomn, can be specified by header
    :param y: default is the 2nd colomn, can be specified by header
    :param group: default is the 3rd colomn, can be specified by header
    :param theme: available themes(Publication, minimal)
    :param horizonal: whether to plot horizonally. (T|F default F)
    :param pallette: Discrete use Set2, heatmap use RdBu, else use Spectral
    :return:
    """
    width = 6
    height = 4.5
    with open(table) as fh:
        header = fh.readline()
        (x1, y1, group1) = header.rstrip().split()
    if x == '':
        x = x1
    if y == '':
        y = y1
    if group == '':
        group = group1
    if pallette == "":
        etc = '+scale_fill_manual(values=c("#E69F00", "#56B4E9", "#8c55e6", "##e64d00"))'
    else:
        etc = '+scale_fill_brewer(palette="{}")'.format(pallette)
    if theme != "":
        etc += "+theme_" + theme + "()"
    if horizonal == 'T':
        etc += '+ coord_flip()'
        # (width, height) = (height, width)
    cmd = theme_publication_r + barplot_r.format(table, x, y, group, etc, width, height)
    rscript(cmd)


pheatmap_sh = r"""
RCG<-read.table("{0}", header = T, row.names = 1, sep="\t")

library("pheatmap")
library("RColorBrewer")

# mycol =colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)

RCGm <- as.matrix(RCG)

RCG <- ifelse(RCGm>5, 5, RCGm)

# 
# col_size = length(colnames(RCG))
# 
# RCG <- RCG[rowSums(RCG[,c(1:col_size)]) > 5]
# 
# RCG <- log2(RCG/rowMeans(RCG))

prefix="{0}"

pdf(paste(prefix, "pdf", sep='.'))

a = pheatmap(RCG, show_rownames=F, 
       main = "{1}", cluster_rows=T, cluster_cols=F,
       fontsize_col = 20, angle_col ="45",border_color = 'white')
"""


def pheatmap(table=None, main=''):
    if main == '':
        main = table
    cmd = pheatmap_sh.format(table, main)
    rscript(cmd)
    return 0


if __name__ == "__main__":
    emain()
