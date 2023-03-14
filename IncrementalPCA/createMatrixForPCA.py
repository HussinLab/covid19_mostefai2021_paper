#!/usr/bin/env python3


import os



import sys
from datetime import datetime
import timeit



import pickle as pk


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import hsv



import seaborn as sns
sns.set_theme(style="ticks", color_codes=True)


from sklearn.decomposition import PCA, IncrementalPCA

import csv
import pandas as pd
#import numba


print("libraries loaded")


########################################
########################################
########  SETTINGS #####################
########################################
########################################



datafile=sys.argv[1]
varfile=sys.argv[2]
minAlt=int(sys.argv[3])
print(varfile)
saved_snpnames=datafile+".snpNames_"+str(minAlt)+"_noref.pkl"
#saved_hash1=datafile+".hash1_"+str(minAlt)+"_noref.pkl"
#saved_hash2=datafile+".hash2_"+str(minAlt)+"_noref.pkl"
saved_df=datafile+".df_"+str(minAlt)+"_noref.pkl"

########################################
########################################
########  SELECTION OF SNPS ############
########################################
########################################


listnuc=["A","C","G","T"]

data=pd.read_csv(varfile, delim_whitespace=True)
data=data[data['POS']>55]
data=data[data['POS']<29804]

n_sample=int(data[1:2]["A"])+int(data[1:2]["C"])+int(data[1:2]["G"])+int(data[1:2]["T"])+int(data[1:2]["."])
too_much_missing=str(len(data[data['.']>int(n_sample/10)]))

print(too_much_missing+ " sites where remove because of more than "+str(n_sample/10)+" missing data (10%)")
data=data[data['.']<=int(n_sample/10)]


#data=data[data['REF']=="A"]
snp_names=[]
for ref in listnuc :
    alt_listnuc=[a for a in listnuc if a!=ref]
    for alt in alt_listnuc:
        d=list(data.loc[(data['REF']==ref) & (data[alt]> minAlt)]["POS"])
        snp_names=snp_names+[str(i)+"_"+ref+"_"+alt for i in d]
nb_snp=len(snp_names)
print(str(nb_snp)+" SNP kept with at least "+str(minAlt)+" alternative allele")
snp_names=sorted(snp_names, key=lambda i: int(i.split("_")[0]))
snplist=[(int(i.split("_")[0]),i.split("_")[2]) for i in snp_names]


#pk.dump(snplist, open(saved_snpnames,"wb"))

########################################
########################################
########  Hashtable 1 construction  ####
########################################
########################################


def list2bin(list_of_bool):
    return(int("".join(map(lambda x : str(int(x)),list_of_bool)), 2))


def transformseq(seq):
    l=[seq[snp[0]-1]==snp[1] for snp in snplist]
    return(list2bin(l))

n=0
hashTable1={}
start = timeit.default_timer()
with open(datafile, 'r') as read_obj:
    csv_reader = csv.reader(read_obj, delimiter='\t')
    for row in csv_reader:
        d=transformseq(row[0])
        if d in hashTable1:
            hashTable1[d]=hashTable1[d]+"|"+row[1]
        else:
            hashTable1[d]=row[1]
        n+=1
        if(n%10000==0):
            stop = timeit.default_timer()
            print(str(n)+" sequences read reaching "+str(len(hashTable1))+" unique sequences added in "+str(stop - start)+" seconds")
            start = timeit.default_timer()
print("step 1 DONE : "+str(n)+" sequences correspond to "+str(len(hashTable1))+" unique sequences")

#pk.dump(hashTable1, open(saved_hash1,"wb"))

########################################
########################################
########  Hashtable 2 construction  ####
########################################
########################################

hashTable2={}

def bin2list(element):
    return [x == '1' for x in format(element, str(nb_snp)+'b')]

kl=list(hashTable1.keys())
n=0
start = timeit.default_timer()
for k in kl:
    kid=hashTable1[k]
    hashTable2[kid]=bin2list(k)
    n+=1
    if(n%10000==0):
        stop = timeit.default_timer()
        print(str(n)+" over "+str(len(kl))+" sequence added in "+str(stop - start)+" seconds")
        start = timeit.default_timer()

print("step 2 DONE : hashtable inverted")



#pk.dump(hashTable2, open(saved_hash2,"wb"))

hashTable1.clear() #for memory

########################################
########################################
########  panda df construction  #######
########################################
########################################

 
data=pd.DataFrame.from_dict(hashTable2, orient='index',columns = snp_names)

print("step 3 DONE : hashtable to matrix")

print(data.iloc[0:10,0:10])
print(saved_df)
pk.dump(data, open(saved_df,"wb"))
