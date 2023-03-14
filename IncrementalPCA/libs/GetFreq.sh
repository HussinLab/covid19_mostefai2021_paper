#!/bin/bash

tsvfile=$1
reference=$2

cat $tsvfile | cut -f1 | tr '[RYSWKMBDHVN\-]' '.'  | awk -v refseq=$(cat $reference  | awk 'substr($1,1,1)!=">"{printf "%s",$1}END{printf "\n"}') '
 BEGIN{
  l=length(refseq)
  split(refseq, ref, "")
  split("ATGC.", temp, "")
  for(i in temp){states[temp[i]]=0}
  for(s in states){for(i=1;i<=l;i++){count[s][i]=0}}
 }
 {
  split($1, seq, "");
  for(i=1;i<=l;i++){count[seq[i]][i]++}
 }
 END{
  printf "POS\tREF"
  for(s in states){printf "\t%s",s}
  printf "\n"
  for(i=1;i<=l;i++){
   printf "%i\t%s",i,ref[i];
   max=0
   for(s in states){
   	c=count[s][i]
   	printf "\t%s",c
   	if(s!="." && s!=ref[i] && max<c){
   	 max=c
   	}
   }
   printf "\n"
  }
 }'
