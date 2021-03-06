library(ape)
library('reshape2')
library('ggplot2')
library('pegas')
library(gridExtra)
library(reticulate)
library(igraph)
library(jsonlite)
library(RColorBrewer)
library('mapplots')
library('plyr')

print('Generate Fig S6')
print('Color by lineage')

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 


# Resize circle size
diviseur=10

# Which lineage to color
lign='Pangolin' #Nextstrain_clade		 Pangolin


pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')
df=read.table(paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),header=TRUE, sep='\t',quote="\"")
colnames(df)=c('ID','Seq','Continent','CountryToCheck','Region','Date','Pangolin_Lineage','Nextstrain_clade','hap')


lineage=read.table(paste(c(pathToOuputData,'info_',version,'.tsv'),collapse=''),header=TRUE, sep='\t',quote="\"")
colnames(lineage)=c('ID','Seq','Continent','CountryToCheck','Region','Date','Pangolin_Lineage','Nextstrain_clade','hap')

if(lign=='Pangolin'){

	# Only need those information
	lineage=lineage[,c('ID','Pangolin_Lineage')]
	lineage$Pangolin_Lineage=as.character(lineage$Pangolin_Lineage)

	# Count so that we can filter those too rare
	c=count(lineage$Pangolin_Lineage)
	c=c[order(c[,2]),]
	keep=c[c[,2]>5000,1]
	
	# To change the name (later) of those too few in number
	Sub=colsplit(as.character(lineage$Pangolin_Lineage),'\\.',c('A','B','C','D'))


	lineage=cbind(lineage,Sub)

	# Change the label of those without a lineage to None
	lineage[lineage[,'A']=='','A']='None'

	# Prepare the name of those who are more rare : only keep one number after the dot
	lineage[,'C']=apply( lineage[ , c('A','B') ] , 1 , paste , collapse = "." )

	# Still applied the filter for the new sub-lineage. Only keep the first letter if still not enought
	f=(count(lineage[,'C']))
	freq=f[f$freq>5000,'x']
	lineage[lineage[,'C']%in%freq,'A']=lineage[lineage[,'C']%in%freq,'C']

	# Format problematic name generated by the paste
	lineage[lineage[,'A']=='B. 1','A']='B.1'
	lineage[lineage[,'A']=='B.NA','A']='B'
	lineage[lineage[,'A']=='A.NA','A']='A'
	lineage[lineage[,'A']=='D.NA','A']='D'
	lineage[lineage[,'A']=='None.NA','A']='None'
	
	# If the filter was enought the first time, keep the full lineage
	lineage[lineage[,'Pangolin_Lineage']%in%keep,'A']=as.character(lineage[lineage[,'Pangolin_Lineage']%in%keep,'Pangolin_Lineage'])
	
}else{
	
	lineage=lineage[,c('ID','Nextstrain_clade')]
	
	lineage$A=as.character(lineage$Nextstrain_clade)
	lineage[is.na(lineage$A) | lineage$A=='','A']='None'


}


# Only one table to manipulate
df=merge(df,lineage,by=1)



pdf=paste(c("../figures/",version,"/"),collapse='')
pdf=paste(c(pdf,"db.NetworksLineage",lign,'Time.',version,'.pdf'),collapse='')
pdf(pdf, width=20,height=20)



net=readRDS(paste(c(pathToOuputData,'net.TimeCoherence.',version,'.',filter,'.rds'),collapse=''))
coordinate=readRDS(paste(c(pathToOuputData,'Coordinate.Time.',version,'.',filter,'.rds'),collapse=''))
d=readRDS(paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
h=readRDS(paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))




#Get matrix with association of each haplotype for each individual
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)


# Associate the color to the lineage
col=as.data.frame(c("#dd4733","#6bd545","#7942ce","#bfdb4d","#ca48bb","#72d589","#4d2773","#d7b942","#6075d1","#628631","#cb8ed8","#d58542","#729fc6","#963d2a","#77d7cc","#cc4072","#548871","#5c262d","#cdc78d","#34354d","#354427","#d58794","#8c6d45","#845c84"))
col=as.data.frame(col[1:length(unique(df$A)),1])
col[,1]=as.character(col[,1])
colnames(col)='Color'

# No lineage = grey
col$Lineage=c(sort(unique(df$A))[sort(unique(df$A))!='None'],'None')
col[col$Lineage%in%c('None','','Unknown'),'Color']='grey'

t=matrix(,nrow=nrow(ind.hap), ncol=0)


if(!file.exists(paste(c(pathToOuputData,'ind.hap.',lign,'.rds'),collapse=''))){

	# Count the number of sample by lineage for each haplotype 
	for (cont in col$Lineage) {
		o=1
		t2=c()

		while(o<=nrow(ind.hap)){ #number of haplotype
			Lineage=as.character(df[df[,'A']==cont,]$ID)
			t2=c(t2,sum(ind.hap[o,ifOfTheLineage]))
			o=o+1
		}
		#Bind lineage
		t=cbind(t,t2)
	}
	saveRDS(t,paste(c(pathToOuputData,'ind.hap.',lign,'.rds'),collapse='')) #Faster next time
}else{
	t=readRDS(paste(c(pathToOuputData,'ind.hap.',lign,'.rds'),collapse=''))
}


#Name haplotype and lineage
colnames(t)=col$Lineage
rownames(t)=rownames(ind.hap)


# Only show labels of haplotype with more than 300 samples
attr(net, "labels")[attr(net, "freq")<300]=''


plot(net, size=sqrt(attr(net, "freq"))/diviseur, scale.ratio=3, pie=t, fast=TRUE, labels=FALSE,threshold=c(0,5),show.mutation=0, bg=alpha(col$Color,1), cex=1.5, lwd=1.5, col='red')
par(new=TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')

legend('topleft', col$Lineage, fill=col$Color,  cex=2.5, bg='white', ncol=3,bty = "n")
legend.bubble(x='bottomleft',z=c(max(rowSums(t))*0.1,max(rowSums(t))*0.4,max(rowSums(t))*0.75,max(rowSums(t))),maxradius=sqrt(max(rowSums(t)))/(diviseur*2), bty='n',txt.cex=2)
# par(new=TRUE)
replot(coordinate)

plot(net, size=sqrt(attr(net, "freq"))/diviseur, scale.ratio=3, pie=t, fast=TRUE, labels=TRUE,threshold=c(0,5),show.mutation=0, bg=alpha(col$Color,1), cex=1.5, lwd=1.5, col='red')
par(new=TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
par(new=TRUE)
replot(coordinate)
