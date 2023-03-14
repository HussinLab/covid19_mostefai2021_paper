#!/usr/bin/env python3
import os
import sys
import csv

import pickle as pk

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, IncrementalPCA

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
sns.set_theme(style="ticks", color_codes=True)

#Not tested with another number of component
n_components = 20


############################
############################
### Function set global variables #######

#dictionnary to label each mutation
refseq = ""
def load_ref(file):
    global refseq
    refseq = "".join(open(file).read().split('\n')[1:])

#dictionnary to label each mutation
Hap_dict = {}
def load_hap_names(file):
    Hap_list = [x.split("\t") for x in open(file).read().split('\n') if x]
    global Hap_dict
    Hap_dict= {x[0]:x[1] for x in Hap_list}
def gethap(s): return(Hap_dict.get(s,"NA"))

#Haplotype positions in the reference
Hap_posList=[]
def load_hap_positions(file):
    global Hap_posList
    Hap_posList = [int(x) for x in open(file).read().split(',') if x]

working_directory=""
def define_working_directory(path):
    global working_directory
    working_directory = path


    
    

########################################
### Function to load matrix Data #######
def loadData(datafile):
    #load the data
    data=np.load(datafile, allow_pickle=True)
    print(data.shape)

    #prepare the hplotypes and sample number (per unique sequence)
    rowinfos=getRowInfo(data)
    print(rowinfos.shape)

    return(data,rowinfos,data.shape[1])


########################################
### Function to load matrix Data #######
def loadData_filterHap(datafile,Haplotypes_set):
    #load the data
    data=np.load(datafile, allow_pickle=True)
    print(data.shape)

    #prepare the hplotypes and sample number (per unique sequence)
    rowinfos=getRowInfo(data,withhap=True)
    print(rowinfos.shape)

    #plot the sample number hitogram 
    #p=plt.hist(rowinfos['sampleNumber'], bins = 20,log=True)


    #Select onnly the haplotypes we want to project
    filterhap=rowinfos['Haplotype'].isin(Haplotypes_set)
    print("reduce "+str(data.shape[0])+" unique sequences to "+
          str(data[filterhap].shape[0])+" unique sequences ")
    print("reduce "+str(sum(rowinfos['sampleNumber']))+
          " sequences to "+str(sum(rowinfos[filterhap]['sampleNumber']))+" sequences ")
    rowinfos=rowinfos[filterhap]
    return(data[filterhap],rowinfos,data.shape[1])

    
def getRowInfo(mydata,withhap=False):
    rowinfos=pd.DataFrame(index=mydata.index)
    if withhap :
        rowinfos["Haplotype"]=gethapdata(mydata)
    rowinfos["sampleNumber"]=[len(i.split("|")) for i in rowinfos.index]
    return(rowinfos)



#######################################################
### Function to create a PCA (sensitive to memory) ####

def createPCA(datafile):

    #data loading and explosion (for PCA )
    data=np.load(datafile, allow_pickle=True)
    print(data.shape)
    data['sampleNumber']=[i.split("|") for i in data.index]
    data = data.explode('sampleNumber')
    data.index=data['sampleNumber']
    data.drop('sampleNumber', inplace=True, axis=1)
    print(data.shape)

    ipca = IncrementalPCA(n_components=n_components, batch_size=1000)
    ipca.fit(data)

    print(ipca.explained_variance_)
    ipcafile=datafile.replace(".pkl","")+"_ipca.pkl"
    pk.dump(ipca, open(ipcafile,"wb"))
    print(" - file <"+datafile+"> is proccess into <"+ipcafile+">")
    
    
    
    
#Return the reference state of a given position
def getposref(pos):
    return refseq[pos]
#Return the state of a given position for a specific sequence (row)
#Take advantage of the col name format pos _ ref _ alt
def getposfromlist(row,cn):
    for c in cn :
        if row[c]: return c[-1]
    return cn[0][-3]


#Return the state of a given position for a list of sequence (rows)
#Compute the subset of column names with the given position
def getposdata(dat,pos):
    pp="^"+str(pos)+"_"
    collist=list(dat.filter(regex=pp,axis=1).columns)
    if collist==[] : print("Cannot retrieve hap : no information for "+pp)
    return(dat.apply (lambda row: getposfromlist(row,collist), axis=1))

#Return a list of haplotype sequence from a matrix and a list of position
def gethapseqdata(dat):
    a=""
    for p in range(len(Hap_posList)):
        a+=getposdata(dat,Hap_posList[p])
    return(a)


#Return a list of haplotype name from a matrix and a list of position
def gethapdata(dat):
    return([gethap(i) for  i in gethapseqdata(dat)])



def projectData(data,ipca_file):
    
    #load the ipca
    ipca = pk.load(open(ipca_file,'rb'))
    #print("explained variance : \n"+str(ipca.explained_variance_))
  
    batchsize=30000
    batch_projections=[]
    for i in range(round(data.shape[0]/batchsize+0.5)):
        (start,end)=(i*batchsize,(i+1)*batchsize)
        print("project "+str(start)+":"+str(end), end = '')
        batch_projections.append(ipca.transform((data)[start:end]))
        print("->ok")
    full_projection=np.concatenate(batch_projections, axis=0)
  
    
    return(full_projection,ipca.explained_variance_)



def bin2list(element):
    return [x == '1' for x in format(element, str(nb_snp)+'b')]
def list2bin(list_of_bool):
    return(int("".join(map(lambda x : str(int(x)),list_of_bool)), 2))
def transformseq(seq):
    l=[seq[snp[0]-1]==snp[1] for snp in snplist]
    return(list2bin(l))



def constructDataForProjection(datafile,datafiletoproject,Haplotypes_set):
    data=np.load(datafile+".pkl", allow_pickle=True)
    snp_names=list(data.columns.values)
    snplist=[(int(i.split("_")[0]),i.split("_")[2]) for i in snp_names]
    nb_snp=len(snplist)
    print(snplist[0:10])

    #load the ipca
    ipca = pk.load(open(datafile+"_ipca.pkl",'rb'))
    print(ipca.explained_variance_)

    hashTable1={}
    n=0
    with open(datafiletoproject, 'r') as read_obj:
        csv_reader = csv.reader(read_obj, delimiter='\t')
        for row in csv_reader:
            d=transformseq(row[0])
            if d in hashTable1:
                hashTable1[d]=hashTable1[d]+"|"+row[1]
            else:
                hashTable1[d]=row[1]
            n+=1
    print("step 1 DONE : "+str(n)+" sequences correspond to "+str(len(hashTable1))+" unique sequences")


    hashTable2={}
    for k in list(hashTable1.keys()):
        kid=hashTable1[k]
        hashTable2[kid]=bin2list(k)
    print("step 2 DONE : hashtable inverted")
    hashTable1.clear() #for memory


    data=pd.DataFrame.from_dict(hashTable2, orient='index',columns = snp_names)
    print("step 3 DONE : hashtable to matrix")


    p=plt.hist(rowinfos['sampleNumber'], bins = 20,log=True)

    rowinfos=getRowInfo(data)
    filterhap=rowinfos['Haplotype'].isin(Haplotypes_set)
    print("reduce "+str(data.shape[0])+" unique sequences to "+
          str(data[filterhap].shape[0])+" unique sequences")
    print("reduce "+str(sum(rowinfos['sampleNumber']))+
          " sequences to "+str(sum(rowinfos[filterhap]['sampleNumber']))+" sequences ")
    rowinfos=rowinfos[filterhap]



    full_projection=projectOnPCA(data[filterhap],ipca)
    expl_var=ipca.explained_variance_
    snp_nb=nb_snp


