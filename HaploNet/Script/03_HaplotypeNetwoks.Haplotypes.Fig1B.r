library(ape)
library('reshape2')
library('ggplot2')
library('pegas')
print('Generate figure 1B')

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 


# Resize circle size
diviseur=10

pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')
df=read.table(paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),header=TRUE, sep='\t',quote="\"")


pdf=paste(c("../figures/",version,"/"),collapse='')
dir.create(pdf)



# Get the network with time correction and its coordinate
net=readRDS(paste(c(pathToOuputData,'net.TimeCoherence.',version,'.',filter,'.rds'),collapse=''))
coordinate=readRDS(paste(c(pathToOuputData,'Coordinate.Time.',version,'.',filter,'.rds'),collapse=''))
d=readRDS(paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
h=readRDS(paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))


# db = start from the database
pdf=paste(c(pdf,"db.Networks.Time.Haplotypes.",version,".pdf"),collapse='')


pdf(pdf, width=20,height=20)



#Get matrix with association of each haplotype for each individual
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)


#To color by haplotype, will create a diagonal matrix
diagonal=diag(nrow(ind.hap))
u=1
t=c()

#Calcul the number of sample by haplotypes and put it in the diagonal
while(u<=nrow(ind.hap)){
	t=c(t,rownames(ind.hap)[u])
	diagonal[u,u]=sum(ind.hap[u,])
	u=u+1
}
diagonal=as.table(diagonal)
rownames(diagonal)=c(t)
colnames(diagonal)=c(t)

# Create the color
color=as.data.frame(rownames(diagonal))
colnames(color)='hap'
order <- attr(net,'labels')


# Put a color only for the most frequent haplotypes (the firsts ones), the others are white
b=c("#003FFF", "#28573F","#FF0000","#FF964F", "#FFCAE4", "#CB5590", "#A6D4FF", "#AD26FF", "#A9A9A9", "#B3F396", "#40E0D0", "#8F1C55","#03C0A4", "#FFD700","#000000","#009900","#E3B778")
a=as.data.frame(as.character(as.roman(seq(1, nrow(color), 1))))
colnames(a)='name'
a$hex='white'
for(y in seq(1,length(b),1)){
	a[y,'hex']=b[y]
}

color=(as.data.frame(attr(net,'labels')))
colnames(color)='name'
color=merge(color,a,by.x='name',by.y='name', all.x=TRUE, sort=FALSE)
color[is.na(color$hex),'hex']='white'
color=color[order(match(color[,1],attr(net,'labels'))),]


# show.mutation : 
# 0 - no information on the line
# 1 - one strike for each mutation between haplotype
# 2 - one dot for each mutation between haplotype
# 3 - number of mutation

#threshold : 0 - don't show alternatif link; c(0,X) show alternatif link who have at most X mutations

plot(net, size=sqrt(attr(net,'freq'))/diviseur, scale.ratio=3, pie=diagonal, fast=TRUE, labels=TRUE,threshold=0,show.mutation=0, main=1, bg=alpha(color$hex,1), col=color$hex)
par(new=TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
par(new=TRUE)
replot(coordinate)


plot(net, size=sqrt(attr(net,'freq'))/1diviseur0, scale.ratio=3, pie=diagonal, fast=TRUE, labels=FALSE,threshold=0,show.mutation=0, main=1, bg=alpha(color$hex,1), col=color$hex)
par(new=TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
par(new=TRUE)
replot(coordinate)
