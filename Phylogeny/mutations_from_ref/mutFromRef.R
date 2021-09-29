library(ggplot2)

toRemove=read.table("outliersDates.treetime.txt",h=F)
df_moreInfos=read.table("allSamples.random1kPerHap.infos.decimal.withMoreInfos.txt",h=T)
df_moreInfos$haplotype = factor(df_moreInfos$haplotype,levels=haps)
df_moreInfos_sub=subset(df_moreInfos,!tip %in% toRemove$V1)
df_moreInfos_sub=subset(df_moreInfos,!name %in% toRemove$V1)
df_moreInfos_sub_sort=df_moreInfos_sub[order(df_moreInfos_sub$haplotype),]

pdf("MutFromRef.1kPerHap.noBadXV.pdf")
ggplot(df_moreInfos_sub,aes(x=as.Date(date),y=distance,color=as.factor(haplotype)))+geom_point(alpha=1)+scale_color_manual(values=colours_hap[as.factor(df_moreInfos_sub$haplotype)],name = "Haplotype")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
