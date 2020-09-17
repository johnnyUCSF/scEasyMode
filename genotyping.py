import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random
import os
import ast


def format_freemuxvcf(freemuxfile,numsamples):
    ### trim off vcf header for reading
    command = "grep -v '##' "+freemuxfile+' > freemux.tmp'
    os.system(command)
    ### open up vcf from freemux and intersect with vcf from quantseq calling
    freemux_full = pd.read_csv('freemux.tmp',sep='\t',dtype='str')
    ### filter out only those snps with different genotypes between samples
    freemux_full['id'] = freemux_full[freemux_full.columns[1]]+'.'+freemux_full[freemux_full.columns[0]]
    freemux_full.index = freemux_full['id']
    freemux_full = freemux_full.drop(['id'],axis=1)
    freemux_full = freemux_full[freemux_full.columns[-numsamples:]]
    ### subset only those that uniquely distinguish between samples
    for column in freemux_full.columns:
        freemux_full[column] = freemux_full[column].str[:3]
    ###    
    freemux_uniq = freemux_full.apply(lambda row: is_uniq(row),axis=1)
    freemux_uniq = freemux_full[freemux_uniq]
    ### return
    return(freemux_uniq)

def format_demuxvcf(demuxfile,freemux_uniq,numsamples):
    ### trim off vcf header for reading
    command = "grep -v '##' "+demuxfile+' > demux.tmp'
    os.system(command)
    ### open up vcf from demux and intersect with vcf from quantseq calling
    demux_full = pd.read_csv('demux.tmp',sep='\t',dtype='str')
    ### filter out only those snps with different genotypes between samples
    demux_full['id'] = demux_full[demux_full.columns[1]]+'.'+demux_full[demux_full.columns[0]]
    demux_full.index = demux_full['id']
    demux_full = demux_full.drop(['id'],axis=1)
    demux_full = demux_full[demux_full.columns[-numsamples:]]
    ###
    for column in demux_full.columns:
        demux_full[column] = demux_full[column].str[:3]
    ###
    demux_uniq = demux_full[demux_full.index.isin(freemux_uniq.index)]
    ### return
    return(demux_uniq)

def is_uniq(row):
    if len(set(row)) == 1: 
        return(False)
    else:
        return(True)

def compare_snps(demux,freemux):
    matches = []
    for item1 in demux:
        for item2 in freemux:
            if item1 == item2:
                matches.append(1)
            else:
                matches.append(0)
    return(matches)

def get_snptotals(demux_uniq):
    snptotals = []
    for column in demux_uniq.columns.tolist():
        tmptotal = len(demux_uniq[demux_uniq[column]!= '.'])
        snptotals.append(tmptotal)
    return(snptotals)
    
def get_assigndict(sums_matrix):
    assign_dict = {}
    ##
    for row in sums_matrix.index:
        freemux_assign = sums_matrix.columns[sums_matrix.loc[row].tolist().index(sums_matrix.loc[row].max())]
        assign_dict[row] = freemux_assign
    ##
    return(assign_dict)

def cleanup():
    os.system('rm freemux.tmp')
    os.system('rm demux.tmp')
    
def write_dict(dictfile,dictionary):
    f = open(dictfile,"w")
    f.write( str(dictionary) )
    f.close()
    ###
    print('saved as: ',dictfile)
    
def read_dict(dictfile):
    file = open(dictfile,'r')
    contents = file.read()
    dictionary = ast.literal_eval(contents)
    file.close()
    ###
    return(dictionary)
    
def match_demux_freemux(demux_uniq,freemux_uniq):
    ####
    sums = len(demux_uniq.columns.tolist())*len(freemux_uniq.columns.tolist())*[0]
    for snpid in demux_uniq.index.tolist():
        ### compare the two snps between freemux and demux
        snpsum = compare_snps(demux_uniq.loc[snpid],freemux_uniq.loc[snpid])
        for i in range(len(snpsum)):
            sums[i] = sums[i]+snpsum[i]
        ### get total possible number of comparisons
        snptotals = get_snptotals(demux_uniq)
    ### reformat
    sums_matrix = []
    indexer = 0
    for i in range(len(demux_uniq.columns.tolist())):
        tmp_row = []
        for j in range(len(demux_uniq.columns.tolist())):
            tmp_row.append(sums[indexer])
            indexer+=1
        ### divide by total possible number of comparisons
        tmp_row = [x/snptotals[i] for x in tmp_row]
        print(snptotals[i],tmp_row)
        sums_matrix.append(tmp_row)
    ### make into df
    sums_matrix = pd.DataFrame(sums_matrix,index=demux_uniq.columns.tolist(),columns=freemux_uniq.columns.tolist())
    sums_matrix = sums_matrix.transpose()
    ### plot
    sns.heatmap(sums_matrix)
    sns.clustermap(sums_matrix,standard_scale=0)
    ### get assignment dictionary
    assign_dict = get_assigndict(sums_matrix)
    ### cleanup
    cleanup()
    ### return
    return(assign_dict)



