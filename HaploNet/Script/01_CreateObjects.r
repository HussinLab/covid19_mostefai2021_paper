library(data.table)
library(stringr)
library(pegas)
library('reshape2')
library('plyr')


freadf=function(...)return(as.data.frame(fread(...)))

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 

# Path to data
pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')

# Fasta containing all sequences
input=paste(c(pathToOuputData,version,'.fasta'),collapse='')

#########################################
#####								#####
#####		Information				#####
#####								#####
#########################################

# File containing the ID of the samples, the Sequence of the haplotype, the Continent, the country, the Region, the Data, the Lineage of Pangolin and Nextstrain clade, and the haplotype number
# In this order
# Could be obtained from the database
data=data=read.table(paste(c(pathToOuputData,'info_',version,'.tsv'),collapse=''), sep='\t',quote="\"", header=TRUE)

colnames(data)=c('ID','Seq','Continent','CountryToCheck','Region','Date','Pangolin_Lineage','Nextstrain_clade','hap')


# Facultatif : to avoid too many lineage
data$Pangolin_Lineage=colsplit(data$Pangolin_Lineage,'\\.',c('A','B'))$A
data$Nextstrain_clade=colsplit(data$Nextstrain_clade,'\\.',c('A','B'))$A

# Read the fasta in an easy manipulable format
sequence=read.table(input) 

# Put the name of the sample on the same line than the haplotype sequence
sequence=cbind(as.character(sequence[seq(1, nrow(sequence), 2),]),as.character(sequence[seq(2, nrow(sequence), 2),]))

# Count the number of sample having the haplotype and filter by the number determined above
c=count(sequence[,2])
c=as.character(c[c[,2]<=filter,1])

# Write a new fasta file containing only the sample having not a rare haplotype 
sequence2=sequence[!(sequence[,2] %in% c),]
write.table(sequence2,paste(c(input,2),collapse=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\n')
print(paste(c('Write : ',input,2),collapse=''))

# Filter the information table for the same filter
data=data[data$Seq %in% sequence2[,2],]

# Format the continent name 
data$Continent=as.character(data$Continent)
data$Continent=gsub('_',' ',data$Continent)


write.table(data,paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),quote=FALSE, row.names=FALSE,col.names=TRUE,sep='\t')
print(paste(c('Write : ',pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''))







##############################################
####							####
####		Create basic haplotype		####
####			network			####
####							####
##############################################
# Create the basic haplonet to have all haplotype labels when correcting for time (see 02)

# Read the new fasta file
print('Sequence read')
d=ape::read.dna(paste(c(input,2),collapse=''), format='fasta')

print('Get haplotype')
h <- pegas::haplotype(d)

print('Get ind.hap')
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)


ind.hap=melt(ind.hap,id=rownames(ind.hap)) #Format so each name will have the sequence next to it
ind.hap=ind.hap[order(-ind.hap$value),] #Sort by frequence

ind.hap=ind.hap[ind.hap$value!=0,] #Remove the non-associated (those with haplotype which did not pass the frequency filter)
ind.hap=ind.hap[!duplicated(ind.hap$hap),] #keep only one haplotype

# Format table to match haplotype sequence and name 
ind.hap=merge(ind.hap[,c('hap','pop')],data[,c('ID','hap')], by.x='pop',by.y='ID', sort=FALSE)
ind.hap=ind.hap[match(rownames(h),ind.hap$hap.x),]

# Name the haplotype with the labels in the network
rownames(h)=ind.hap$hap.y
net <- pegas::haploNet(h)

print('Writing objects for network in R')
saveRDS(h,paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))
print(paste(c(pathToOuputData,'haplotype.',version,'.',filter,'.rds'),collapse=''))

saveRDS(net,paste(c(pathToOuputData,'net.',version,'.',filter,'.rds'),collapse=''))
print(paste(c(pathToOuputData,'net.',version,'.',filter,'.rds'),collapse=''))
print(net)

saveRDS(d,paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
print(paste(c(pathToOuputData,'fasta.',version,'.',filter,'.rds'),collapse=''))
