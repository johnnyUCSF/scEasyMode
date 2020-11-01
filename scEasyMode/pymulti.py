#!/usr/bin/env python3

import pandas as pd
import gzip
import sys
import Bio
import pickle
import datetime
import timeit
import os
from collections import Counter
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import zscore, norm
from sklearn.neighbors import DistanceMetric
from scipy.spatial.distance import hamming

sys.setrecursionlimit(10000)

def split_rawdata(R1,R2,len_10x,len_umi,len_multi,sampname,huge,base_start=0):
    """ this reads in fastq files and definitions of read structure and dumps it as a pickle """
    ###store paired end reads
    reads = []
    ###open paired end files simultaneously
    with gzip.open(R1, "rt") as F1, gzip.open(R2, "rt") as F2:
        for record1,record2 in zip(SeqIO.parse(F1, "fastq"),SeqIO.parse(F2, "fastq")):
            ####extract info from reads
            bc_10x = str(record1.seq[:len_10x])
            umi = str(record1.seq[len_10x:(len_10x+len_umi)])
            if len_multi==None:
                r2 = str(record2.seq[base_start:len_multi])
            else:
                r2 = str(record2.seq[base_start:base_start+len_multi])
            ####write in
            reads.append([bc_10x,umi,r2])
    ###if huge file, then use joblib (pickle crashes)
    if huge == True:
        print('file is huge. not saving.')
        return(reads)
    else:
    ###regular pickle dump
        pickle.dump(reads, open('pymulti/'+sampname+"_reads.p", "wb" ) )
        return(reads)

def read_pickle(sampname,reads,huge):
    """ this reads in the pickle data written from split_rawdata """
    if huge == True:
        readtable = pd.DataFrame(reads)
    else:
        ####read in
        readtable = pd.DataFrame(pickle.load(open('pymulti/'+sampname+"_reads.p", "rb" ) ))
    ####format
    readtable.columns = ['cell','umi','multi']
    readtable.set_index('cell')
    ####return
    return(readtable)

'''
sequencing cleanup functions
'''

def matchby_hamdist(readtable,bcsmulti,bcs10x):
    """ hamming distance on barcodes. runs in roughly N time. hamming for 10x barcodes takes a long time. """
    #######this matches each item by the hamming distance threshold to the known barcodes for multi and 10x
    ####make translation table
    trantab = make_trantab()
    ####apply to multiseq barcodes
    if bcsmulti != None:
        readtable['multi'] = readtable.apply(lambda row: find_closest(str(row['multi']),bcsmulti,trantab),axis=1)
        print('finished hamming correction for multiseq.',timeit())
#     ####apply to 10x barcodes
#     readtable['cell'] = readtable.apply(lambda row: find_closest(str(row['cell']),bcs10x,trantab),axis=1)
#     print('finished hamming correction for 10x.',timeit())
#     ####return
#     return readtable

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

'''
check and filter functions
'''

def check_stats(readtable,bcsmulti,bcs10x):
    """ this checks the stats using exact match to multiseq, duplicates, and 10x """
    ####calculate # reads that contain a multiseq barcode
    if bcsmulti != None:
        multirate = 100*(len(readtable[readtable.multi.isin(bcsmulti)]))/len(readtable)
        print("Multiseq rate is: ",multirate)
    ####calculate # reads that are in the 10x bc index
    cellrate = 100*len(readtable[readtable['cell'].isin(bcs10x)])/len(readtable)
    print("10X Cell Barcode rate is: ",cellrate)
    ####calculate duplication rate
    tmp = readtable.drop_duplicates()
    duprate = 100*(len(readtable)-len(tmp))/len(readtable)
    print("Duplication rate is: ",duprate)

def filter_readtable(readtable,bcsmulti,bcs10x,gbc_thresh):
    """ this filters the readtable with exact matches """
    ####subset readtable
    filtd = readtable.drop_duplicates()
    if bcsmulti != None:
        filtd = filtd[filtd.multi.isin(bcsmulti)]
    filtd = filtd[filtd['cell'].isin(bcs10x)]
    if gbc_thresh != None:
        ####
        keep = filtd['multi'].value_counts()>(np.median(filtd['multi'].value_counts())*gbc_thresh)
        keep = keep[keep==True].index.tolist()
        filtd = filtd[filtd.multi.isin(keep)]
        print('Number of GBCs: ',len(keep))
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

'''
multiseq correction and cell calling
'''

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

def correct_simple(filtd,sampname, plots=True,thresh=False,pct_only=False,thresh_dict={}, filter_unmapped_reads=True):
    """ correct by zscore distributions, plot distributions, and call cells. Saves a file with the calls. """
    ###pre run checks
    pivot = format_multi_table(filtd)
    check_pivot(pivot,sampname)
    if thresh==True: correct_cutoffs(pivot,thresh_dict,cut=correct_simple(filtd,sampname,thresh=False))
    ###raw percentiles
    test = pivot+0.01
    test = test/test.sum()
    if plots==True: sns.clustermap(test,cmap='Blues')
    ###option to use pct only if thresholding
    if pct_only != True:
        ###log2
        test = test.apply(np.log2)
        if plots==True: sns.clustermap(test,cmap='Blues')
        ###zscore
        test = test.transpose().apply(zscore)
    else:
        test = test.transpose()
    # plot distributions
    if plots==True:
        for pop in list(test.columns):
            plt.figure()
            test[pop].hist(bins=100,color=np.random.rand(3,))
            plt.title(pop)
            plt.xlim(-5,5)
    # get significance value and call cells
    new = test.copy()
    new['call'] = test.apply(lambda row: get_sig(list(row),list(test.columns)),axis=1)
    new['sig'] = test.apply(lambda row: z_ratio(list(row)),axis=1)
    
    if filter_unmapped_reads:
        old_len = len(new)
        column_names = new.columns
        result_columns = ['cell', 'call', 'sig']
    
        barcodes = [x for x in column_names if x not in result_columns]
        new = new[new[barcodes].sum(axis=1) != 0]
        new_len = len(new)
        print('Removed {} cells that were not mapped to any provided barcode.'.format(old_len - new_len))
    
    # save
    print('saving file ',sampname)
    new.to_csv('pymulti/'+sampname+'_calls.tsv',sep='\t')
    # return in the event of multiple iterations of corrections
    return(test)

def correct_median(filtd,sampname, med_factor, plots=True, filter_unmapped_reads=True):
    # pre run checks
    pivot = format_multi_table(filtd)
    check_pivot(pivot,sampname)
    # subtract the median multiplied by the factor for each barcode (ie remove background)
    testmed = pivot.transpose()-med_factor*(pivot.transpose().apply(np.median))
    ###formatting, set negatives to zero
    testmed = testmed.transpose()
    testmed[testmed<0] = 0
    ###calculate pct per cell
    testmed = testmed/testmed.sum()
    ###formatting
    testmed = testmed.fillna(0)
    if plots == True:
        ###plot
        sns.clustermap(testmed,cmap='Blues',method='weighted',metric='cityblock',row_cluster=False)
    ###formatting
    test = testmed.transpose()
    ####get significance value and call cells
    new = test.copy()
    new['call'] = test.apply(lambda row: get_sig(list(row),list(test.columns)),axis=1)
    new['sig'] = test.apply(lambda row: z_ratio(list(row)),axis=1)
    
    # filter unmapped reads
    
    if filter_unmapped_reads:
        old_len = len(new)
        column_names = new.columns
        result_columns = ['cell', 'call', 'sig']
    
        barcodes = [x for x in column_names if x not in result_columns]
        new = new[new[barcodes].sum(axis=1) != 0]
        new_len = len(new)
        print('Removed {} cells that were not mapped to any provided barcode.'.format(old_len - new_len))
     
    # save
    new.to_csv('pymulti/'+sampname+'_calls.tsv',sep='\t')
    # return in the event of multiple iterations of corrections
    return(test)

'''
main function
'''

def pymulti(R1,R2,bcs10x, len_10x=16, len_umi=12, len_multi=8, base_start=0, med_factor=1.6, gbc_thresh=None, sampname='pymulti_',
            split=True, plots=True, hamming=False, thresh=False, pct_only=False, median_only=False, huge=False, bcsmulti=None,
            reads=None, thresh_dict={}, filter_unmapped_reads=True):
    """ 
    main loop, splits from fastqs and runs through cell calls
        R1 = your Read1 fastq for the multiseq/hashing fraction
        R2 = your Read2 fastq for the multiseq/hashing fraction
        bcsmulti = your whitelist of known multiseq/hashing barcodes
        bcs10x = your whitelist of known cell identifiers
        len10x = the length of the 10x barcodes
        len_umi = the length of the umis on each pair of reads
        len_multi = the length of the multiseq barcode sequence the package will split out
        med_factor = if doing median correction, this is the factor that is multiplied by the median to estimate the amount of background reads per barcode
        sampname = this is the handle with which intermediate files will be saved
        split = True if you have not split it from fastq yet, False if you have split it already and want to rerun the processing. this will save time.
        hamming = True if you want to use match multiseq barcodes within hamming distance of 1 to the multiseq/hashing whitelist
        thresh = True if you want to use a threshold on which to gate reads, different correction method
        pct_only = True if you only want to use the raw percentage of reads to call the barcode
        median_only = True if you want to use the median correction method
        huge = True if your file is huge in which case it will be saved as an h5 file instead of pickling
        thresh_dict = a dictionary of barcodes:threshold which will be used to gate reads on the samples
        
    """
    ###
    if huge == True: print('assuming huge fastqs.')
    # split fastqs and pickle
    if not os.path.isdir('./pymulti'):
        os.mkdir('pymulti')
    
    if split == True: reads = split_rawdata(R1,R2,len_10x,len_umi,len_multi,sampname,huge=huge,base_start=base_start)
    # read in old pickle data
    readtable = read_pickle(sampname,reads=reads,huge=huge)
    #####check duplication multi and 10x rates
    multirate = check_stats(readtable,bcsmulti,bcs10x)
    # match by hamming distance
    if hamming == True:
        readtable.drop_duplicates(inplace=True)
        matchby_hamdist(readtable,bcsmulti,bcs10x)
        #####check duplication multi and 10x rates
        multirate = check_stats(readtable,bcsmulti,bcs10x)
    # filter readtable
    filtd = filter_readtable(readtable,bcsmulti,bcs10x,gbc_thresh)
    ####implement multiseq correction by within distribution zscores
    if median_only == True:
        correct_median(filtd, sampname, med_factor, plots, filter_unmapped_reads=filter_unmapped_reads)
    else:
        correct_simple(filtd, sampname, plots, thresh, pct_only, thresh_dict, filter_unmapped_reads=filter_unmapped_reads)
    return(filtd)
