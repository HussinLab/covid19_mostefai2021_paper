library(ape)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(ggtree)

#code is the same for both phylogenetic tree, only the input file for the tree is changing
df_dates=read.table("../allSamples.random1kPerHap.infos.decimal.txt",h=T,sep="\t")
mytree_div=read.nexus("divergence_tree.nexus")
p <- ggtree(mytree_div)
p2 <- p %<+% df_dates + geom_tippoint(aes(color=haplotype))
pdf("allHaplotypes.1kperHap.divergenceTree.colorsHaps.pdf")
p2+scale_color_manual(values=colours_hap[as.factor(df_dates$haplotype)])+ theme_tree2()
dev.off()
