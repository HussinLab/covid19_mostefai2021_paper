#Extract from MSA positions that are not to mask
maskfile=problematic_sites_sarsCov2.2021_04_15.vcf

#extract only ref from alignment :
grep -A1 "EPI_ISL_402119" allSamples.random1kPerHap.fasta > ref.fasta
myref=ref.fasta

#get the positions on the MSA that match to the reference, so we can easily extract known positions
cat $myref | tail -1 | awk '{split($1,seq,"");n=1;for(i=1;i<=length($1);i++){if(seq[i]!="-"){print n,i;n++}}}' > allPosWithIndex_fromAWK.txt

#get the positions we want to mask
grep "mask" $maskfile | sed '/^#/d' | cut -f 2 > posToMask.2021_04_15.txt

#print the contrary, the ones we want to keep
#cat posToMask.2021_04_15.txt | awk 'BEGIN{s=1;e=29903;i=1}$1!=i{printf "%i-%i\n",i,$1-1;i=$1+1}$1==i{i++}' | tr '\n' ',' | sed 's/,$//' > toKeepFromCambridge.txt

#from the MSA, print the positions that matches the columns to remove
echo "posList=read.table(\"allPosWithIndex_fromAWK.txt\",h=F)
toExtract=read.table(\"posToMask.2021_04_15.txt\",h=F)
posList_kept=subset(posList,V1 %in% toExtract\$V1)
write.table(posList_kept\$V2,file=\"columnsToRemove.txt\",quote=F,row.names=F)" > extract.R
R CMD BATCH extract.R

#convert to bed format
sed '1d' columnsToRemove.txt | awk '{print "sarscov2\t"$1-1"\t"$1}' > columnsToRemove.bed

#get the size of the MSA
size=$(tail -n 1 ref.fasta | awk '{print length($0)}')


#print to bed format so we can get the reverse intersection
for f in $( seq 1 $size ) ; do echo "$f" >> allPos.txt ; done
awk '{print "sarscov2\t"$1-1"\t"$1}' allPos.txt > allPos.bed
rm allPos.txt

module load bedtools
bedtools subtract -a allPos.bed -b columnsToRemove.bed > columnsToKeep.bed
rm allPos.bed


#transform the bed file into intervals all onto one line
cut -f3 columnsToKeep.bed | awk 'NR==1{s=$1;i=$1+1}NR!=1 && $1!=i{printf "%i-%i\n",s,i-1;s=$1;i=$1+1}$1==i{i++}END{printf "%i-%i\n",s,i-1}' | tr '\n' ',' | sed 's/,$/\n/' > columnsToKeep.txt


#extract from the MSA the good columns
cat allSamples.random1kPerHap.fasta |awk 'substr($1,1,1)==">"{printf  "\n%s\n",$0}substr($1,1,1)!=">"{printf "%s",$1}' | awk 'NR!=1' |  awk 'NR%2==1' | tr ' ' '_' > onlyid.txt
cat allSamples.random1kPerHap.fasta  |awk 'substr($1,1,1)==">"{printf  "\n%s\n",$0}substr($1,1,1)!=">"{printf "%s",$1}' | awk 'NR!=1' | awk 'NR%2==0' | cut -c$(cat columnsToKeep.txt) > onlyseq.txt

paste onlyid.txt onlyseq.txt | awk '{printf "%s\n%s\n",$1,$2}' > allSamples.random1kPerHap.goodPos.fasta


