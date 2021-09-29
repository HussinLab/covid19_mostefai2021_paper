library(ape)
library('reshape2')
library('ggplot2')
library('pegas')
library('mapplots')

print('Generate Fig 4A and others')
print('Color for each mutation')

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 


# Resize circle size
diviseur=10



pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')
df=read.table(paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),header=TRUE, sep='\t',quote="\"")


# Position in the genome of the virus (can be obtained with a vcf)
position=as.data.frame(sort(c(241,313,1059,1163,3037,7540,8782,14408,14805,16647,18555,22227,22992,23063,23401,23403,25563,26144,28144,28881,28882,28883)))


pdf=paste(c("../figures/",version,"/"),collapse='')
pdf=paste(c(pdf,"db.NetworksMutationsTime.",version,".pdf"),collapse='')


net=readRDS(paste(c(pathToOuputData,'net.TimeCoherence.',version,'.',filter,'.rds'),collapse=''))
coordinate=readRDS(paste(c(pathToOuputData,'Coordinate.Time.',version,'.',filter,'.rds'),collapse=''))
d=readRDS(paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
h=readRDS(paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))


pdf(pdf, width=20,height=20)


# Length of haplotype
longueurHaplotype=nrow(position)



#Get matrix with association of each haplotype for each individual
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)



#One color by haplotype
diagonal=diag(nrow(ind.hap))
u=1
t=c()

#Calcul the number of sample for the haplotype
while(u<=nrow(ind.hap)){
	t=c(t,rownames(ind.hap)[u])
	diagonal[u,u]=sum(ind.hap[u,])
	u=u+1
}
diagonal=as.table(diagonal)
rownames(diagonal)=c(t)
colnames(diagonal)=c(t)

# Create a table for the color
color=as.data.frame(rownames(diagonal))
colnames(color)='hap'
color$col='white'

# Keep all haplotypes only once
haplotypeSequence=df[!duplicated(df$Seq),c('Seq','hap')]
colnames(haplotypeSequence)=c('pop','hap')


# Split for each mutation
# Keep the last 3 in a codon
sequence=colsplit(string=haplotypeSequence$pop,pattern='',names=seq(1, nrow(position)-2, 1))

rownames(sequence)=haplotypeSequence$hap

s=sequence
c=color
	

lh=longueurHaplotype-2

for (i in colnames(sequence)) {
	sequence=s
	color=c

	# Color for the nucleotide base
	myColorsS=c('yellow','blue','red','green')
	names(myColorsS)=c('T','A','C','G')

	# Extract those needed
	myColorsS=myColorsS[unique(sequence[,i])]

	if(i==lh){
		# For the codon, choose random color
		myColorsS=c('grey','#AD26FF','#00FFFFFF','#80FF00FF')
		names(myColorsS)=unique(sequence[,as.character(lh)])
	}

	# Order haplotype to be the same as the network
	sequence=sequence[match(color$hap,rownames(sequence)),]
	

	color=merge(color,sequence,by.x='hap',by.y=0,sort=FALSE)
	color[,as.character(i)]=as.character(color[,as.character(i)])

	# Attribute the good color for each haplotype
	for (j in names(myColorsS)) {
		color[color[,as.character(i)]==j,as.character(i)]=myColorsS[j]
	}

	attr(net, "labels")[attr(net, "freq")<300]=''

	
	plot(net,size=sqrt(attr(net, "freq"))/diviseur,  scale.ratio=3, pie=diagonal, fast=TRUE, labels=FALSE,threshold=0,show.mutation=0, main=paste(c('Identification of the mutation : ',position[i,1]),collapse=''), bg=alpha(color[,as.character(i)],0.8), col=color[,as.character(i)])
	par(new=TRUE)
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
	par(new=TRUE)
	replot(coordinate)
	legend('bottomright',names(myColorsS),fill=myColorsS, cex=2)

}
