library(ape)
library('reshape2')
library('ggplot2')
library('pegas')
library('mapplots')

print('Generate Fig 2A and B')

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 


# Resize circle size
diviseur=10

# Date to separe first and second wave
dateSeparate='2020-08-01'
month=months(as.Date(dateSeparate))


pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')
df=read.table(paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),header=TRUE, sep='\t',quote="\"")
# Remove samples without continental information (should not have any)
df=df[!is.na(df$Continent),]

# Get the network with time correction and its coordinate
net=readRDS(paste(c(pathToOuputData,'net.TimeCoherence.',version,'.',filter,'.rds'),collapse=''))
coordinate=readRDS(paste(c(pathToOuputData,'Coordinate.Time.',version,'.',filter,'.rds'),collapse=''))
d=readRDS(paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
h=readRDS(paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))


pdf=paste(c("../figures/",version,"/"),collapse='')
pdf=paste(c(pdf,"db.NetworksTime.Continents.Split.",version,".pdf"),collapse='')
pdf(pdf, width=20,height=20)


#Get matrix with association of each haplotype for each individual
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)


#Colorblind friendly
col=as.data.frame(c("#009E73", "#56B4E9", "#999999", "#D55E00", "#CC79A7", "#F0E442"))
colnames(col)='Color'
col$Continent=c('Asia','Europe','Africa','North America','South America','Oceania')

# If the information for the date is only the year, we can't associate with the first or second wave
df=df[df$Date!='2020' & !is.na(df$Date),]

s=colsplit(string=df$Date,pattern='-',names=c('Y','M','D'))
dateSeparate=colsplit(string=dateSeparate,pattern='-',names=c('Y','M','D'))

# Split first and second wave
before=df[(s$Y<dateSeparate$Y | (s$Y==dateSeparate$Y & s$M<dateSeparate$M)),]
after=df[!(s$Y<dateSeparate$Y | (s$Y==dateSeparate$Y & s$M<dateSeparate$M)),]



t=matrix(,nrow=nrow(ind.hap), ncol=0)

if(!file.exists(paste(c(pathToOuputData,'ind.hap.Continent.first.rds'),collapse=''))) {
	# Count the number of sample by haplotype for each continent
	for (cont in col$Continent) {
		o=1
		t2=c()

		while(o<=nrow(ind.hap)){ #number of haplotype
			idOfTheContinent=as.character(before[before[,'Continent']==cont,]$ID)

			if(sum(ind.hap[o,idOfTheContinent])!=0){
				t2=c(t2,sum(ind.hap[o,idOfTheContinent]))
			}else{
				# If no sample have this haplotype in the continent for the first wave, still need to have a number to have a circle
				t2=c(t2,0.01)
			}			

			o=o+1
		}
		#Bind all continents
		t=cbind(t,t2)
	}
	saveRDS(t,paste(c(pathToOuputData,'ind.hap.Continent.first.rds'),collapse='')) # Will be faster the second time
	
}else{
	t=readRDS(paste(c(pathToOuputData,'ind.hap.Continent.first.rds'),collapse=''))
}

#Name haplotype and continent
colnames(t)=col$Continent
rownames(t)=rownames(ind.hap)

# Keep all labels for the second wave
ta=attr(net, "labels")

# Remove labels for less frequent haplotypes
attr(net, "labels")[rowSums(t)<300]=''



plot(net, size=sqrt(rowSums(t))/diviseur, scale.ratio=3, pie=t, labels=TRUE, fast=TRUE, threshold=0,show.mutation=0, bg=alpha(col$Color,1), cex=1.5, lwd=1.5, col='red', main=paste(c('Before ',month,' 2020'),collapse=''))
par(new=TRUE)

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
par(new=TRUE)
replot(coordinate)
legend('topleft', col$Continent, fill=alpha(col$Color,1),  cex=2.5, bg='white',bty = "n")
legend.bubble(x='bottomleft',z=c(max(rowSums(t))*0.1,max(rowSums(t))*0.4,max(rowSums(t))*0.75,max(rowSums(t))),maxradius=sqrt(max(rowSums(t)))/(diviseur*2), bty='n',txt.cex=2)



# Put back all labels
attr(net, "labels")=ta

t=matrix(,nrow=nrow(ind.hap), ncol=0)

# Do the same thing as before, but for the second wave
if(!file.exists(paste(c(pathToOuputData,'ind.hap.Continent.second.rds'),collapse=''))){
	# Count the number of sample by haplotype for each continent
	for (cont in col$Continent) {
		o=1
		t2=c()

		while(o<=nrow(ind.hap)){ #number of haplotype


			idOfTheContinent=as.character(after[after[,'Continent']==cont,]$ID)
			
			if(sum(ind.hap[o,idOfTheContinent])!=0){
				t2=c(t2,sum(ind.hap[o,idOfTheContinent]))
			}else{
				# If no sample have this haplotype in the continent for the first wave, still need to have a number to have a circle
				t2=c(t2,0.01)
			}
			o=o+1
		}
		#Bind continent
		t=cbind(t,t2)
	}
	saveRDS(t,paste(c(pathToOuputData,'ind.hap.Continent.second.rds'),collapse=''))
}else{
	t=readRDS(paste(c(pathToOuputData,'ind.hap.Continent.second.rds'),collapse=''))
}


#Name haplotype and continent
colnames(t)=col$Continent
rownames(t)=rownames(ind.hap)


# Remove labels for less frequent haplotypes
attr(net, "labels")[rowSums(t)<300]=''



plot(net, size=sqrt(rowSums(t))/diviseur, scale.ratio=3, pie=t, labels=TRUE, fast=TRUE, threshold=0,show.mutation=0, bg=alpha(col$Color,1), cex=1.5, lwd=1.5, col='red', main=paste(c('After 1 ',month,' 2020'),collapse=''))
par(new=TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
par(new=TRUE)
replot(coordinate)
legend('topleft', col$Continent, fill=alpha(col$Color,1),  cex=2.5, bg='white',bty = "n")
legend.bubble(x='bottomleft',z=c(max(rowSums(t))*0.1,max(rowSums(t))*0.4,max(rowSums(t))*0.75,max(rowSums(t))),maxradius=sqrt(max(rowSums(t)))/(diviseur*2), bty='n',txt.cex=2)


# plot(net, size=sqrt(rowSums(t))/diviseur, scale.ratio=3, pie=t, fast=TRUE, labels=TRUE,threshold=0,show.mutation=0, bg=alpha(col$Color,1), cex=1.5, lwd=1.5, col='red', main=paste(c('After ',month,' 2020'),collapse=''))
# par(new=TRUE)
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white", border='white')
# legend('topleft', col$Continent, fill=alpha(col$Color,1),  cex=2.5, bg='white')

# legend.bubble(x='bottomleft',z=c(max(rowSums(t))*0.1,max(rowSums(t))*0.4,max(rowSums(t))*0.75,max(rowSums(t))),maxradius=sqrt(max(rowSums(t)))/(diviseur*2), bty='n',txt.cex=2)
# # par(new=TRUE)
# replot(coordinate)




