library(ape)
library('plyr')
library('dplyr')
library('reshape2')
library('pegas')
library(igraph)
library(jsonlite)

# Number of samples having the haplotype to keep the haplotype
filter=50

# Date of the version
version='2021_03_12' 

# Split by half a month
dataToCheck=c(15,31)

pathToOuputData=paste(c('../Data/database/',version,'/'),collapse='')

df=read.table(paste(c(pathToOuputData,'NewInformationsTable.',version,'.',filter,'.txt'),collapse=''),header=TRUE, sep='\t',quote="\"")
df=df[,c('ID','Date','Seq','hap','Continent')] #Only need the date information

# So that we can filter
s=colsplit(string=df$Date,pattern='-',names=c('Y','M','D')) #Split to date

# Normally, after databaser, should not have problem with missing date
s[is.na(s$M) ,'M']=12 #If only the year known
s[(is.na(s$D) | s$D=='XX' | s$D=='') ,'D']=28 #If the month is know but not the day, it will be considering as an end of the month's sample

df[,'Date']=apply(s[,c('Y','M','D')],1,paste,collapse='-') #Modify the new date
df$Date2=as.Date(df$Date)


firstAppear=c()

#For each haplotype, will find the first time its appeared
for(i in unique(df$hap)){
	t=df[df[,'hap']==i,]
	t2=t[which.min(t$Date2)[1],]
	firstAppear=rbind(firstAppear,t2)
	
}

firstAppear=firstAppear[!is.na(firstAppear[,1]),] #Not supposed to remove something, but juste in case : If no date, so NA will appear


#Again, split the information
s=colsplit(string=firstAppear$Date2,pattern='-',names=c('Y','M','D'))
write.table(firstAppear,'First.txt',quote=FALSE, col.names=TRUE, sep='\t', row.names=FALSE)
# Can be used the next time : it will be faster

# Put either the first of the month or the 15 depending on which section the haplotype first appeared
s[s$D<=dataToCheck[1],'D']=1
s[s$D<=dataToCheck[2] & s$D>dataToCheck[1],'D']=(dataToCheck[1])


#Update the Date information
firstAppear$Date2=apply(s[,c('Y','M','D')],1,paste,collapse='-')
firstAppear$Date2=as.Date((firstAppear$Date2), format="%Y-%m-%d") #Convert in date


#Just to have one table to manipulate
firstAppear=cbind(firstAppear,s)
firstAppear$F='>' #To create temporary fasta files
firstAppear$F=apply(firstAppear[,c('F','hap')],1,paste,collapse='')

test=NULL #Will have all the keeped links
allnet=readRDS(paste(c(pathToOuputData,'net.',version,'.',filter,'.rds'),collapse='')) #To get all the information except the main link (size and order of the haplotype for graph), alter.link won't be good after modification

findCycle=function(mainLinkOrigin){ #This function find cycle after the merge. 
	
	# Rowname too long for my confort when verifying the table
	row.names(mainLinkOrigin)=NULL

	if(nrow(mainLinkOrigin)<3){
		# If less than 3 haplotypes, not supposed to have a cycle
		return(mainLinkOrigin)
	}

	# How this function search for a cycle :
	# 1. Cut a link between two haplotypes
	# 2. Enumerate all haplotypes linked to those two haplotypes (in variable side1 and side2)
	# 3. If there is at least one haplotype commun, this mean that those two haplotypes are in a cycle. Add them in the variable inCycle

	side1=c()
	side2=c()
	inCycle=c()
	hasCycle=FALSE

	# Starting from the first link in the network, we will verified all link
	link=1
	while(link<=nrow(mainLinkOrigin)){
		mainLink=mainLinkOrigin
		
		#Get the 2 haplotypes to split
		hapGroup1=mainLink[link,1]
		hapGroup2=mainLink[link,2]

		#Remove the link of the table
		mainLink=mainLink[-link,]
		group1=c(hapGroup1)
		group2=c(hapGroup2)
		

		last=-1 #Will stop when the script can't find new haplotype in either group
		while(length(group2)+length(group1)!=last){ #If no new haplotypes appeared, stop searching
		

			#Get the number of haplotype split in group (linked to the two haplotypes)
			last=length(group2)+length(group1) 


			#Will find all haplotypes who have a link with an haplotype of either of the group
			group1=c(group1,mainLink[mainLink[,1] %in% group1,2],mainLink[mainLink[,2] %in% group1,1])
			group2=c(group2,mainLink[mainLink[,1] %in% group2,2],mainLink[mainLink[,2] %in% group2,1])

			#Remove dublicated in a same group
			group1=group1[!duplicated(group1)]
			group2=group2[!duplicated(group2)]

			#If there is haplotype commun, still continue until all haplotypes in cycle are found
			if(sum(group1 %in% hapGroup2) || sum(group2 %in% hapGroup1)){
				inCycle=c(inCycle,hapGroup1,hapGroup2)
				hasCycle=TRUE
			}

			
		}
		link=link+1
	}
	if(!hasCycle){ # No cycle found
		return(mainLinkOrigin)
	}
	inCycle=unique(inCycle)
	# Cycle found : remove one
	return(removeCycle(inCycle,mainLinkOrigin))


}

removeCycle=function(inCycle,mainLinkOrigin){
	# This function will remove only one link

	# Only need haplotype and date 
	firstAppear2=firstAppear[,c('hap','Date2')]

	# Will extract haplotypes in cycle
	mainLinkCycle=mainLinkOrigin[mainLinkOrigin[,1] %in% inCycle & mainLinkOrigin[,2] %in% inCycle,]

	# Will extract haplotypes not in cycle 
	mainLinkOrigin=mainLinkOrigin[!(mainLinkOrigin[,1] %in% inCycle & mainLinkOrigin[,2] %in% inCycle),]

	# Will keep in memory the time difference between each haplotype in a cycle
	diff=c()
	for (i in seq(1, nrow(mainLinkCycle), 1)) {
		date1=firstAppear2[as.character(firstAppear2$hap)==as.character(mainLinkCycle[i,1]),'Date2']
		date2=firstAppear2[as.character(firstAppear2$hap)==as.character(mainLinkCycle[i,2]),'Date2']
		diff=c(diff,abs(as.numeric(date1-date2)))
	}

	# Specific case
	# In case the an haplotype has lots of mutations for two links, this is not probable, so we will remove one of those link (Australian haplotype has this problem). This section target specificaly the Australian haplotype and won't work for any other situation... (except if the situation is the same)
	# Can be removed if desire
	###
	diffStep=mainLinkCycle[,'step']
	if(sum(diffStep>=3)!=0 && sum(diffStep==max(diffStep))>1){
		# Found link with high number of mutation differences
		linkSame=c(mainLinkCycle[diffStep==max(diffStep),1],mainLinkCycle[diffStep==max(diffStep),2])
		diff2=diff
		# Will focus only on the highly differenciated haplotypes
		diff2[diffStep!=max(diffStep)]=0
		if(length(linkSame)!=length(unique(linkSame))){ #If there is an haplotype which has many link with large number of mutations
			# Remove the link which time difference is the largest
			mainLinkCycle=mainLinkCycle[-which.max(diff2),]
			diff=diff[-which.max(diff2)]
		}
	}
	###
	
	


	if(sum(diff==max(diff))>1){
		# If there is many links sharing the same time difference
		
		# Extract link sharing the highest time difference
		mainLinkCycle2=mainLinkCycle[which(diff==max(diff)),]

		# Remove them from the cycle table
		mainLinkCycle=mainLinkCycle[-which(diff==max(diff)),]

		# Extract the number of mutations of difference
		diff=mainLinkCycle2[,'step']

		# Remove the link having the higher time difference and mutations 
		mainLinkCycle2=mainLinkCycle2[-which.max(diff),]

		# Add the other links to the table
		mainLinkCycle=rbind(mainLinkCycle,mainLinkCycle2)

	}else{
		mainLinkCycle=mainLinkCycle[-which.max(diff),]
	}

	
	# Bind all links (from the cycle of not)
	mainLinkOrigin=rbind(mainLinkOrigin,mainLinkCycle)

	# Possible that we removed all cycles, but we need to check
	return(findCycle(mainLinkOrigin))

}



# Main function!
for (i in sort(unique(firstAppear$Y))) { #For each year in the dataset
	print(i)
	if(is.na(i)){ #Not suppose to be any, but just in case...
		next
	}	
	for (y in sort(unique(firstAppear[firstAppear$Y==i,]$M))) { #For each month in the year
		print(y)
		if(is.na(y)){ #Not suppose to be any, but just in case
			next
		}		

		for(date in dataToCheck){#For each date to check
			print(date)
			#Will get the line of the individual in the datafile who was sequenced before the date
			# The year before | the month before of this year | in the month until the date | date unknown but month before | month unknown 
			columnToKeep=c((firstAppear$Y<i | firstAppear$Y==i & firstAppear$M<y | (firstAppear$Y==i & firstAppear$M==y & firstAppear$D<=date & !is.na(firstAppear$D))  | (firstAppear$Y<=i & firstAppear$M<y & is.na(firstAppear$D))  | (firstAppear$Y<=i & is.na(firstAppear$M))))


			if(sum(columnToKeep)<=1){
				# If only one haplotype in this timeframe, skip
				next
			}
			
			# Extract sample from the time frame and write a fasta file 
			upTime=firstAppear[columnToKeep,]
			write.table(upTime[,c('F','Seq')],'Temp2.fasta',quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\n')


			d <- ape::read.dna('Temp2.fasta', format='fasta')
			sequence=read.table('Temp2.fasta') 

			# Get the name of the sample
			nam=sequence[grep('>',sequence[,1]),]

			# Get the sequence
			sequence=as.data.frame(sequence[!(sequence[,1] %in% nam),])
			
			names(d)=nam
			rownames(d)=nam

			h <- pegas::haplotype(d)

			ind.hap<-with(
				stack(setNames(attr(h, "index"), rownames(h))),
				table(hap=ind, pop=rownames(d)[values])
			)

			#Format so each name will have the sequence next to it
			ind.hap=melt(ind.hap,id=rownames(ind.hap))
			ind.hap=ind.hap[order(-ind.hap$value),]#Sort by frequence
			ind.hap=ind.hap[ind.hap$value!=0,] #Remove the non-associated (so that we keep only association between name and sequence)

			# Format name
			ind.hap$pop=gsub('>','',ind.hap$pop)
			ind.hap=ind.hap[match(rownames(h),ind.hap$hap),]
			rownames(h)=ind.hap$pop

			# Create the haplonet for the time frame
			net <- pegas::haploNet(h, getProb=FALSE)
			print(net)

			# Get main link
			net2=as.data.frame(net[,1:ncol(net)])
			
			# When too few haplotypes, I got a bug
			if(ncol(net2)==1){
				net2=t(net2)
				colnames(net2)=colnames(net[,1:ncol(net)])
			}
	
			# Since the labels is not the same when we generate a new haplonet, need to change the number by the labels so that we can associated the same haplotype between time frame
			for (l in seq(1, length(attr(net,'labels')), 1)) {
				net2[net2[,1]==l,1]=(attr(net,'labels')[l])
				net2[net2[,2]==l,2]=(attr(net,'labels')[l])
			}
			

			net2=as.matrix(net2, ncol=4)
			colnames(net2)=c('','','step','Prob')

			# Merge links of the two timeframe
			test=(rbind(test,net2))
			
			# Remove same link
			test=unique(test)
			colnames(test)=c('','','step','Prob')

			# Get all name-haplotypes and format so that the order of the haplotype is the same that in the non-corrected network
			allName=unique(c(test[,1],test[,2]))
			allName=allName[match(attr(allnet,'labels'),allName)]
			allName=allName[!is.na(allName)]



			u=1
			if(nrow(test)>1){
				# Find cycles after the merging
				test=findCycle(test)
			}

			
			# If we have already all samples in our graph, no need to continue
			if(sum(columnToKeep)==nrow(firstAppear)){
				break
			}
			
		}	

	}
}

# Change the labels to the numeric one
while(u<=length(allName)){
	test[test[,1]==allName[u],1]=as.numeric(u)
	test[test[,2]==allName[u],2]=as.numeric(u)
	u=u+1
}

# Recreate the haplonet object, use information from the non-corrected network
test=as.matrix(test,ncol=4)
mode(test)='numeric'

colnames(test)=c('','','step','Prob')

# Since the order is the same as the non-corrected network, get the frequencies by haplotype
attr(test, "freq") <- attr(allnet, "freq")
attr(test, "labels") <- allName
attr(test, "alter.links") <- NULL
class(test) <- "haploNet"


saveRDS(test,paste(c(pathToOuputData,'net.TimeCoherence.',version,'.',filter,'.rds'),collapse=''))

#Write main links
write.table(net[,c(1:4)], paste(c(pathToOuputData,'time.MainLinks.',version,'.',filter,'.txt'),collapse=''),quote=FALSE, row.names=FALSE,col.names=TRUE, sep='\t')




# To determine the position of each haplotypes in the graph (We give our coordinates) : 
print('Do this command please in the terminal')
print('R')
print('library(pegas)')
print(paste(c("net=readRDS('",pathToOuputData,"net.TimeCoherence.",version,'.',filter,".rds')"),collapse=''))
print(paste(c("coordinate=readRDS('",pathToOuputData,"Coordinate.Time.",version,'.',filter,".rds')"),collapse=''))
print("plot(net,labels=TRUE,threshold=0,show.mutation=0, size=sqrt(attr(net, 'freq'))/10, scale.ratio=3, fast=TRUE)")
print("replot(coordinate)")
print('o=replot()')
print('#Here, change the position of the circles in the graphic')
print(paste(c("saveRDS(o,'",pathToOuputData,"Coordinate.Time.",version,'.',filter,".rds')"),collapse=''))
print('dev.off()')
print('quit()')
