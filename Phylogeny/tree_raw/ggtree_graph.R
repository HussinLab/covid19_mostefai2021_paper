library(ape)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(ggtree)
colours_hap=c("#003FFF", "#28573F","#FF0000","#FF964F", "#FFCAE4", "#CB5590", "#A6D4FF", "#AD26FF", "#A9A9A9", "#B3F396", "#40E0D0", "#8F1C55","#03C0A4", "#FFD700","#000000","#009900","#E3B778")
haps=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","XVII")
names(colours_hap)=haps

#code is the same for both phylogenetic tree, only the input file for the tree is changing
df_dates=read.table("../allSamples.random1kPerHap.infos.decimal.tsv",h=T,sep="\t")
mytree_div=read.nexus("divergence_tree.nexus")
p <- ggtree(mytree_div)
p2 <- p %<+% df_dates + geom_tippoint(aes(color=haplotype))
pdf("allHaplotypes.1kperHap.divergenceTree.colorsHaps.pdf")
p2+scale_color_manual(values=colours_hap[as.factor(df_dates$haplotype)])+ theme_tree2()
dev.off()
