library(ggplot2)

df_tempest=read.table("tempEst.rootToTip.txt",h=T)
df_dates=read.table("allSamples.random1kPerHap.infos.decimal.txt",h=T,sep="\t")
names(df_tempest)=c("tip","date_tempest","distance","residual")
df_tempest_merge=merge(df_tempest,df_dates,by.x="tip",by.y="name")
df_tempest_merge$haplotype = factor(df_tempest_merge$haplotype,levels=haps)

pdf("RootToTip.perHap.pdf")
ggplot(df_tempest_sub,aes(x=as.Date(date),y=distance,color=as.factor(haplotype)))+geom_point(alpha=1)+scale_color_manual(values=colours_hap[as.factor(df_tempest_sub$haplotype)],name = "Haplotype")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
