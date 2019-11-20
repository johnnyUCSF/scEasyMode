#!/usr/bin/env python3
"""
Module pymulti
"""

__author__ = "Johnny Yu"
__version__ = "0.1.0"
__license__ = "MIT"


#####################
#####################
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import pickle
from collections import Counter
import timeit
import copy
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import zscore
from scipy.stats import norm
from sklearn.neighbors import DistanceMetric
import datetime
from scipy.spatial.distance import hamming
import gzip
#####################
#####################

#####################
#####################readin from fastq functions
#####################

def split_rawdata(R1,R2,len_10x,len_umi,len_multi,sampname):
    """ this reads in fastq files and definitions of read structure and dumps it as a pickle """
    ###store paired end reads
    reads = []
    ###open paired end files simultaneously
    with gzip.open(R1, "rt") as F1, gzip.open(R2, "rt") as F2:
        for record1,record2 in zip(SeqIO.parse(F1, "fastq"),SeqIO.parse(F2, "fastq")):
            ####extract info from reads
            bc_10x = str(record1.seq[:len_10x])
            umi = str(record1.seq[len_10x:(len_10x+len_umi)])
            r2 = str(record2.seq[:len_multi])
            ####write in
            reads.append([bc_10x,umi,r2])
    ###dump
    pickle.dump(reads, open('pymulti/'+sampname+"_reads.p", "wb" ) )
    
def read_pickle(sampname):
    """ this reads in the pickle data written from split_rawdata """
    ####read in
    readtable = pd.DataFrame(pickle.load(open('pymulti/'+sampname+"_reads.p", "rb" ) ))
    ####format
    readtable.columns = ['cell','umi','multi']
    readtable.set_index('cell')
    ####return
    return(readtable)

#####################
#####################sequencing cleanup functions
#####################

def matchby_hamdist(readtable,bcsmulti,bcs10x):
    """ hamming distance on barcodes. runs in roughly N time. hamming for 10x barcodes takes a long time. """
    #######this matches each item by the hamming distance threshold to the known barcodes for multi and 10x
    ####make translation table
    trantab = make_trantab()
    ####apply to multiseq barcodes
    readtable['multi'] = readtable.apply(lambda row: find_closest(str(row['multi']),bcsmulti,trantab),axis=1)
    print('finished hamming correction for multiseq.',timeit())
    ####apply to 10x barcodes
    readtable['cell'] = readtable.apply(lambda row: find_closest(str(row['cell']),bcs10x,trantab),axis=1)
    print('finished hamming correction for 10x.',timeit())
    ####return
    return readtable
    
def find_closest(searchstr,barcodes,trantab):
    """ finds closest in distance matrix  based on hamming <= 1"""
    if searchstr in barcodes:
        return(searchstr)
    else:
        #######this matches and replaces the item in the row by the multiseq barcode within hamthresh
        dists = fast_hamm(searchstr,barcodes,trantab)
        ###return the barcode with the lowest hamdist if hamdist == 1
        i = dists.index(min(dists))
        if dists[i] <= 1:
            return(barcodes[i])
        else:
            return(searchstr) 
    
def make_trantab():
    """ make translation table for encoding """
    intab = 'AGCTN'
    outtab = '01234'
    trantab = str.maketrans(intab, outtab)
    ###return
    return(trantab)
    
def encode_str(searchstr,trantab):
    """ encoding for the string to numeric for faster """
    encoded = np.array([int(d) for d in searchstr.translate(trantab)])
    ###return
    return(encoded)

def fast_hamm(searchstr,barcodes,trantab):
    """ Make sure that searchstring and searchlist contain equal length of chars before using. 
        Not going to check in here to make faster. 
        Returns hamming distance matrix """
    full_len = len(searchstr)
    encoded = encode_str(searchstr,trantab)
    dists = [hamming(encoded,encode_str(barcode,trantab))*full_len for barcode in barcodes]
    ###return hamming distances
    return(dists)

def pairwise_hamming(pd_series):
    """ used to calculate pairwise hamming across a stack. """
    trantab = make_trantab()
    full_len = len(pd_series.iloc[0])
    ###encode
    encoded = pd_series.map(lambda x: np.array([int(d) for d in x.translate(trantab)]))
    encoded = np.stack(encoded, axis=0)
    ###hamming
    distf = DistanceMetric.get_metric('hamming')
    dist = distf.pairwise(encoded)*full_len
    ###mask diagonals
    np.fill_diagonal(dist, 1000)
    ###find closest
    closest = [np.where(row==np.amin(row))[0] for row in dist]
    ###return
    return(dist,closest)

#####################
#####################check and filter functions
#####################
    
def check_stats(readtable,bcsmulti,bcs10x):
    """ this checks the stats using exact match to multiseq, duplicates, and 10x """
    ####calculate # reads that contain a multiseq barcode
    multirate = 100*(len(readtable[readtable.multi.isin(bcsmulti)]))/len(readtable)
    print("Multiseq rate is: ",multirate)
    ####calculate # reads that are in the 10x bc index
    cellrate = 100*len(readtable[readtable['cell'].isin(bcs10x)])/len(readtable)
    print("10X Cell Barcode rate is: ",cellrate)
    ####calculate duplication rate
    tmp = readtable.drop_duplicates()
    duprate = 100*(len(readtable)-len(tmp))/len(readtable)
    print("Duplication rate is: ",duprate)
    
def filter_readtable(readtable,bcsmulti,bcs10x):
    """ this filters the readtable with exact matches """
    ####subset readtable
    filtd = readtable.drop_duplicates()
    filtd = filtd[filtd.multi.isin(bcsmulti)]
    filtd = filtd[filtd['cell'].isin(bcs10x)]
#     ####check size and shrink to max
#     max_cells = 30000
#     keep = filtd.groupby(by='cell').count()
#     keep.sort_values(by='umi',inplace=True)
#     keep = keep.iloc[-max_cells:].index.tolist()
#     filtd = filtd[filtd.index.isin(keep)]
    ####return
    return(filtd)

def timeit():
    print(datetime.datetime.now())
    
#####################
#####################multiseq correction and cell calling
#####################

def get_sig(row,labels):
    """ get the highest significance value """
    i = row.index(max(row))
    return(labels[i])

def z_ratio(row):
    """ get significance value by zscore difference """
    row.sort()
    top1 = row[-1]
    top2 = row[-2]
    return(top1-top2)

def format_multi_table(filtd):
    """ format multiseq table into a stacked pivot table by cell """
    multi = filtd.groupby(by=["cell","multi"]).count()
    pivot = multi.pivot_table(index='multi', columns='cell', values='umi')
    pivot = pivot.fillna(0)
    return(pivot)

def check_pivot(pivot,sampname):
    """ check before implementing corrections """
    if len(pivot) == 1:
        test = (test.apply(np.log2)).transpose()
        test.to_csv(sampname+'_calls.tsv',sep='\t')
        raise NameError('Everything is one barcode. Calls are printed by log2 counts.')
    else:
        return True

def correct_cutoffs(pivot,thresh_dict,cut):
    """ reduce to zero based on above distributions """
    pivot = pivot.transpose()
    for barcode in list(pivot.columns):
        pivot[barcode] = cut.apply(lambda row: 0.0 if row[barcode]<thresh_dict[barcode] 
                                   else pivot.loc[row.name][barcode],axis=1)
    pivot = pivot.transpose()
    
def correct_simple(filtd,sampname,thresh=False,thresh_dict={},plots):
    """ correct by zscore distributions, plot distributions, and call cells. Saves a file with the calls. """
    ###pre run checks
    pivot = format_multi_table(filtd)
    check_pivot(pivot,sampname)
    if thresh==True: correct_cutoffs(pivot,thresh_dict,cut=correct_simple(filtd,sampname,thresh=False))
    ###raw percentiles
    test = pivot+0.01
    test = test/test.sum()
    if plots==True: sns.clustermap(test,cmap='Blues')
    ###log2
    test = test.apply(np.log2)
    if plots==True: sns.clustermap(test,cmap='Blues')
    ###zscore
    test = test.transpose().apply(zscore)
    ####plot distributions
    for pop in list(test.columns):
        plt.figure()
        test[pop].hist(bins=100,color=np.random.rand(3,))
        plt.title(pop)
        plt.xlim(-5,5)
    ####get significance value and call cells
    new = test.copy()
    new['call'] = test.apply(lambda row: get_sig(list(row),list(test.columns)),axis=1)
    new['sig'] = test.apply(lambda row: z_ratio(list(row)),axis=1)
    ###save
    new.to_csv('pymulti/'+sampname+'_calls.tsv',sep='\t')
    ###return in the event of multiple iterations of corrections
    return(test)

#####################
#####################main function
#####################

def pymulti(R1,R2,bcsmulti,bcs10x,len_10x=16,len_umi=12,len_multi=8,sampname='pymulti_',
            split=True,hamming=False,thresh=False,thresh_dict={},plots=True):
    """ main loop, splits from fastqs and runs through cell calls """
    ###split fastqs and pickle
    os.system('mkdir pymulti')
    if split == True: split_rawdata(R1,R2,len_10x,len_umi,len_multi,sampname)
    ###read in old pickle data
    readtable = read_pickle(sampname)
    #####check duplication multi and 10x rates
    multirate = check_stats(readtable,bcsmulti,bcs10x)
    ####match by hamming distance
    if hamming == True:
        readtable.drop_duplicates(inplace=True)
        matchby_hamdist(readtable,bcsmulti,bcs10x)
    ####filter readtable
    filtd = filter_readtable(readtable,bcsmulti,bcs10x)
    ####implement multiseq correction by within distribution zscores
    correct_simple(filtd,sampname,thresh,thresh_dict,plots)
    