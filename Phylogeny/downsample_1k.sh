#For every haplotypes of interest, extract the samples information from a metadata file, look at the dates and make groups. For a downsampling of 1000 sample per haplotype, we would want, for the period of time we look at, 3 samples per day.
#so the sample_n function in dplyr will randomly extract 3 samples.

for f in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII
do
echo "library(dplyr)

df=read.table(\"infos/sorted.$f.infos.only2020.noOutlier.txt.match\",h=F)
df=df[-1,]

new_df <- df %>% group_by(V2) %>% sample_n(3,replace=T)
new_df=unique(new_df)
new_df=as.data.frame(new_df)

write.table(new_df\$V1,file=\"sorted.$f.infos.only2020.noOutlier.shuff3PerDate.txt\",row.names=F,col.names=F,quote=F)" > $f.shuf3.Rscript

R CMD BATCH $f.shuf3.Rscript

#take a look at the total number of samples extracted and randomly take the remaining samples from the ones that were not picked so we can have 1000 samples in total
i=$(wc -l sorted.$f.infos.only2020.noOutlier.shuff3PerDate.txt | cut -f 1 -d ' ')
j=$(expr 1000 - $i)
grep -Fwvf sorted.$f.infos.only2020.noOutlier.shuff3PerDate.txt infos/sorted.$f.infos.only2020.noOutlier.txt.match | shuf | head -n $j | cut -f 1 > $f.remaining$j.txt 

cat sorted.$f.infos.only2020.noOutlier.shuff3PerDate.txt $f.remaining$j.txt > $f.random1k.txt

done

rm *shuff3PerDate.txt
rm *remaining*
mv *.Rscript* scripts/

mv *.random1k.txt random_1k_fasta


#extract form the fasta file the selected samples
for f in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII
do
/lustre03/project/6005588/shared/bin/seqtk/seqtk subseq ../../sorted.$f.noOutliers.noReplicate.fasta random_1k_fasta/$f.random1k.txt > random_1k_fasta/$f.random1k.fasta
done
