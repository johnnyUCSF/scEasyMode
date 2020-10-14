#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import scanpy as sc
from datetime import datetime
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats.mstats import gmean

def read_species(human=True):
    """ Returns either the human or mouse anndata object
        human = True if you want the human object, false will give you the mouse object
    """
    if human == True:
        return(sc.read_h5ad("DataBySpecies/human.anndata.h5ad"))
    else:
        return(sc.read_h5ad("DataBySpecies/mouse.anndata.h5ad"))

def save(adata,filename):
    adata.write(filename+".h5ad")

def read(filename):
    adata = sc.read_h5ad(filename+".h5ad")
    return(adata)

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

#######################################################
####################################################### Preliminary metadata functions
#######################################################

def overlay_meta(adata,LMOfile,PYMULTIfile):
    """ Automatically overlays the metadata from multiseq to the single cell dataset.
        adata = the anndata object
        LMOfile = the file mapping the cell hashing/multiseq indices to the sample identifiers
        sampname = the sample name that has been used to save files as the prefix

        Returns an anndata object with metadata from the LMOfile overlaid
    """
    ###get barcodes
    adata.obs['barcode'] = adata.obs.index.str[:-2]
    ###get dictionaries of metadata
    LMOdict = get_LMOfile(LMOfile)
    sig,call = get_pymulti(PYMULTIfile)
    ###check indices
    adata = adata[adata.obs.barcode.isin(sig.keys())]
    adata = adata[adata.obs.barcode.isin(call.keys())]
    ###write into .obs
    adata.obs['sig'] = adata.obs.apply(lambda row: sig[row.barcode],axis=1)
    adata.obs['call'] = adata.obs.apply(lambda row: call[row.barcode],axis=1)
    ###write sample name into .obs
    adata.obs['sample'] = adata.obs.apply(lambda row: LMOdict[row.call],axis=1)
    ###return
    return(adata)

def overlay_custom_meta(adata,METAfile,label):
    """ Overlays the metadata specified in METAfile. Note that it must have the barcodes as colname = 'BARCODE' """
    ###read file
    mfile = pd.read_csv(METAfile,sep=',',index_col=0)
    ###convert to dictionary
    mdict = mfile.to_dict()[mfile.columns[0]]
    ###write into .obs
    adata.obs[label] = adata.obs.apply(lambda row: mdict.get(row.barcode),axis=1)
    ###return
    return(adata)

def overlay_gbcs(adata,gbc_df,label1='sig_gbc',label2='call_gbc'):
    ###convert to dictionary
    mdict1 = gbc_df.to_dict()[gbc_df.columns[0]]
    mdict2 = gbc_df.to_dict()[gbc_df.columns[1]]
    ###check indices
    adata = adata[adata.obs.barcode.isin(mdict1.keys())]
    ###write into .obs
    adata.obs[label1] = adata.obs.apply(lambda row: mdict1[row.barcode],axis=1)
    adata.obs[label2] = adata.obs.apply(lambda row: mdict2[row.barcode],axis=1)
    ###return
    return(adata)

def overlay_demux_freemux(adata):
    """ matches freemuxlet demuxlet """
    #### format
    multi = adata.obs[['freemux','demux','barcode']]
    multi = multi.groupby(by=["freemux","demux"]).count()
    pivot = multi.pivot_table(index='freemux', columns='demux', values='barcode')
    pivot = pivot.fillna(0)
    test = pivot.copy()
    test = test.transpose()
    test = test/test.sum()
    #### plot
    sns.clustermap(test)
    #### write to file
    new = test.copy()
    new['call'] = test.apply(lambda row: get_sig(list(row),list(test.columns)),axis=1)
    new['sig'] = test.apply(lambda row: z_ratio(list(row)),axis=1)
    ### save
    sampname = 'freemux.demux'
    print('saving file ',sampname)
    new.to_csv(sampname+'.sig.'+'_calls.tsv',sep='\t')
    #### read into adata object
    new = new.reset_index()
    new.index = new.call
    newdict = new.demux.to_dict()
    adata.obs['cell_match'] = adata.obs.apply(lambda row: newdict[row.freemux],axis=1)
    return(adata)

def format_demuxlet(DEMUXfile):
    """ formats the demuxlet file to two columns. """
    demux = pd.read_csv(DEMUXfile,sep='\t')
    demux.index = demux.BARCODE.str[:-2]
    demux = demux[['SNG.1ST']]
    fname = DEMUXfile+".csv"
    demux.to_csv(fname,sep=',')
    print('Your file has been saved as:',fname)

def format_freemuxlet(FREEMUXfile):
    """ formats the freemux file to two columns. """
    demux = pd.read_csv(FREEMUXfile,sep='\t')
    demux.index = demux.BARCODE.str[:-2]
    demux = demux[['BEST.GUESS']]
    demux['BEST.GUESS'] = demux['BEST.GUESS'].str[:1]
    fname = FREEMUXfile+".bestguess.csv"
    demux.to_csv(fname,sep=',')
    print('Your file has been saved as:',fname)
    demux = pd.read_csv(FREEMUXfile,sep='\t')
    demux.index = demux.BARCODE.str[:-2]
    demux = demux[['DROPLET.TYPE']]
    fname = FREEMUXfile+".droplettype.csv"
    demux.to_csv(fname,sep=',')
    print('Your file has been saved as:',fname)

def get_LMOfile(LMOfile):
    """ read in LMOfile and turn into dictionary.
        using the multiseq barcode sequence as the keys. """
    bcsmulti = pd.read_csv(LMOfile,sep=',',index_col=0,header=None)
    bcsmulti.columns = ['multi']
    return(bcsmulti.to_dict()['multi'])

def get_pymulti(PYMULTIfile):
    """ read in python multi seq file and convert to two dictionaries, one for significance and one for call.
        using the cell barcodes as the keys. """
    bcs = pd.read_csv(PYMULTIfile,sep='\t',index_col=0)
    bcs = bcs[['sig','call']]
    sig = bcs.sig.to_dict()
    call = bcs.call.to_dict()
    return(sig,call)

#######################################################
####################################################### Basic QC
#######################################################

def qcstats(adata,MT_prefix):
    ###stat
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    ###calculate # mitochondrial
    mito_genes = adata.var.index.str.startswith(MT_prefix)
    adata.obs['mt_frac'] = np.sum(
        adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    ###return
    return(adata)

def qcplots(adata):
    #Sample quality plots
    t1 = sc.pl.violin(adata, 'n_genes', groupby='sample', size=2, log=True, cut=0)
    t1 = sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0)
    t2 = sc.pl.violin(adata, 'mt_frac', groupby='sample')
    #Data quality summary plots
    colors = ['mt_frac','sample']
    for color in colors:
        print('plotting '+color+' below.')
        sc.pl.scatter(adata, 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<5000], 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<2500], 'n_counts', 'n_genes', color=color)

def qccheck(adata,suppress_plots,MT_prefix='MT-'):
    adata = qcstats(adata,MT_prefix)
    if suppress_plots == False:
        qcplots(adata)
    return(adata)

def qccounts(adata,threshold=4000):
    #Thresholding decision: counts
    sns.distplot(adata.obs['n_counts'])
    plt.figure()
    sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<threshold],bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>threshold],bins=100)
    plt.figure()

def qcgenes(adata,threshold=1000):
    #Thresholding decision: genes
    sns.distplot(adata.obs['n_genes'], bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<threshold], bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']>threshold], bins=100)
    plt.figure()

def cellcycle(adata,mouse):
    adata.var_names_make_unique()
    ##Score cell cycle and visualize the effect:
    cell_cycle_genes = [x.strip() for x in open('scEasyMode/resources/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    ##add labels for mouse genes
    if mouse == True:
        s_genes = ['MM10___'+s for s in s_genes]
        g2m_genes = ['MM10___'+g2m for g2m in g2m_genes]
    ##calculate
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return(adata)

def qc_all(adata,suppress_plots=False,mouse=False):
    """ Performs standard quality checks on data and visualizes it. Does not do any filtering.
        Checks mitochondrial genes, number of genes, and number of counts
        adata = the anndata object to operate on
        suppress_plots = True if you want to just do the qc without all the plots

        Returns formatted anndata object.
    """
    ###qc
    adata = qccheck(adata,suppress_plots)
    adata.var.index = adata.var.index.str.upper()
    adata = cellcycle(adata,mouse)
    ###plot
    if suppress_plots == False:
        qccounts(adata)
        qcgenes(adata)
    return(adata)

#######################################################
####################################################### Basic Filtering
#######################################################

def annotate_mito(adata,mt_thresh):
    adata.obs['dead'] = adata.obs.apply(lambda row: row.mt_frac>mt_thresh,axis=1)
    return(adata)

def filter_cells(adata,mt_thresh,min_counts=0,max_counts=10000000,min_genes=0):
    ###Filter cells according to identified QC thresholds:
    print('Total number of cells: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, min_counts = min_counts)
    print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, max_counts = max_counts)
    print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
    adata = adata[adata.obs['mt_frac'] < mt_thresh]
    print('Number of cells after MT filter: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, min_genes = min_genes)
    print('Number of cells after gene filter: {:d}'.format(adata.n_obs))
    return(adata)

def filter_genes(adata,min_cells=20):
    print('Total number of genes: {:d}'.format(adata.n_vars))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
    return(adata)

def normalize(adata):
    ###Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()
    ###normalize
    sc.pp.normalize_total(adata,exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    ##Store the full data set in 'raw' as log-normalised data for statistical testing
    adata.raw = adata
    ##return
    return(adata)

def filters(adata,mt_thresh,min_cells,min_counts,max_counts,min_genes,sig_gbc=False,get_hvgs=True,sig_pct=False):
    """ implements many filters at the same time to get you the relevant cells and genes.
        adata = the anndata object to operate on
        mt_thresh = the % mitochondrial above which to filter out cells
        min_cells = the number of cells that a gene has to be expressed in to keep
        min_counts = the number of counts that a cell must have to keep
        max_counts = the number of counts that a cell must be below to keep
        min_genes = the number of genes that a cell must have to keep
        get_hvgs = select HVGs in preparation for downstream clustering
        sig_pct = the percentile of which to remove cells based on the multiseq/hashing

        Returns two formatted anndata objects. The first one is the original without mitochondrial gene filtering, the second one is the one with mitochondrial gene filtering.
    """
    ###annotate dead cells in full data set
    adata = annotate_mito(adata,mt_thresh)
    ###filter by sig from multiseq calls
    if sig_pct != False:
        print('cutting off by sig_pct.')
        adata = adata[adata.obs.sig>sig_pct]
    if sig_gbc != False:
        print('cutting off by sig_gbc.')
        adata = adata[adata.obs.sig_gbc>sig_gbc]
    ###filter mitochondrial high genes out and create clean dataset
    clean = filter_cells(adata,mt_thresh,min_counts=min_counts,max_counts=max_counts,min_genes=min_genes)
    ###filter genes
    adata = filter_genes(adata,min_cells)
    clean = filter_genes(clean,min_cells)
    ###normalize
    adata = normalize(adata)
    clean = normalize(clean)
    ###get hvgs
    if get_hvgs == True:
        adata = define_hvgs(adata)
        clean = define_hvgs(clean)
    ###get shapes
    print(adata.shape)
    print(clean.shape)
    ###return
    return(adata,clean)

def define_hvgs(adata,n_genes=3000):
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_genes)
    print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    sc.pl.highly_variable_genes(adata)
    return(adata)

def filter_by_sig(adata,sig_pct):
    adata = adata[adata.obs.sig>np.percentile(adata.obs.sig,sig_pct)]
    return(adata)

#######################################################
####################################################### Clustering and visualization
#######################################################

def visualize(adata,covariates=['n_counts','n_genes','mt_frac','phase','sample','leiden','dead','sig'],res=0.5,bbknn=False,suppress_plots=False):
    """ calculates visualizations for all cells in pca, umap, diffusion map, and force directed graph 2D space.
        adata = the anndata object on which to act on.
        covariates = the covariates on which to overlay onto the visualizations
        res = the resolution for louvain clustering (higher is more discrete clusters)
        bbknn = True if you want to use bbknn normalization as a batch correction measure
        suppress_plots = True if you just want to calculate the visualizations without all the plots

        Returns an anndata object with all the calculations built in for plotting.
    """
    ###Calculate the visualizations
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    if bbknn == True:
        sc.external.pp.bbknn(adata,batch_key='batch',set_op_mix_ratio=0.1,trim=100)
    else:
        sc.pp.neighbors(adata)
    sc.tl.umap(adata,random_state=10)
    sc.tl.diffmap(adata)
    sc.tl.draw_graph(adata)
    ###cluster
    adata = cluster(adata,res)
    ###plot
    if suppress_plots == False:
        for covariate in covariates:
            sc.pl.pca_scatter(adata, color=covariate,palette=sc.pl.palettes.vega_20)
            sc.pl.umap(adata, color=covariate,palette=sc.pl.palettes.vega_20)
            sc.pl.diffmap(adata, color=covariate, components=['1,2','1,3'],palette=sc.pl.palettes.vega_20)
            sc.pl.draw_graph(adata, color=covariate,palette=sc.pl.palettes.vega_20)
    ###return
    return(adata)

def regress(adata,factors):
    sc.pp.regress_out(adata, factors)
    sc.pp.scale(adata)
    return(adata)

def cluster(adata,res):
    sc.tl.leiden(adata, resolution=res, key_added='leiden', random_state=10)
    return(adata)

def densitymap(adata,sample_key='sample'):
    """ Plots the density of samples across umap space
        adata = anndata object to act on
        key = the key on which to split the dataset and visualize samples individually.

        Returns a plot to the interpreter.
    """
    sc.tl.embedding_density(adata, basis='umap', groupby=sample_key)
    sc.pl.embedding_density(adata, basis='umap', key='umap_density_'+sample_key)

def classify(df,label,newlabel,classifier):
    df[newlabel] = df.apply(lambda row: lambda_classify(row[label],classifier) ,axis=1)
    return(df)

def lambda_classify(search,classifier):
    if search in classifier:
        return('yes')
    else:
        return('no')

def metacell_bylabel(scobject,label,zscore=True):
    clusters = []
    metacells = []
    ###extract metacells
    for cluster in set(scobject.obs[label]):
        ##get cluster
        clust = scobject[scobject.obs[label]==cluster].X
        ##average to metacell
        clust = clust.mean(axis=0)
        ##append to index
        metacells.append(pd.DataFrame(clust)[0].tolist())
        clusters.append(cluster)
    ###format
    merged = pd.DataFrame(metacells,index=clusters,columns=scobject.var.index.tolist())
    ###zscore
    if zscore == True:
        merged = pd.DataFrame(stats.zscore(merged),index=clusters,columns=scobject.var.index.tolist())
    ###return
    return(merged)

#######################################################
####################################################### Demuxlet correction
#######################################################


def correct_demuxlet(adata,cutoff=0.7,label='cell_type'):
    """ This function corrects integrates demuxlet and transcriptome calls by assigning a cell type to each cluster in transcriptome space.
        It uses cutoffs as defined in parameters which refers to proportion of max cell types. """
    cutoff = 0.7
    celldict = {}
    ###
    for clust in set(adata.obs.leiden):
        ##get cluster alone first
        clustdf = adata[adata.obs.leiden==clust].obs
        ##calculate proportions
        props = clustdf.groupby(label).count()['barcode']/(len(clustdf))
        ##see if majority based on cutoff
        majorcell = props[props>cutoff].index.tolist()
        ##change all cell_type values to the most abundant value
        cells = clustdf.index.tolist()
        if len(majorcell) == 1:
            ##add to dictionary
            celldict[clust] = majorcell[0]
        else:
            celldict[clust] = 'delete'
        ##print a log
        print('leiden cluster ',clust,' has a majority cell type of ',majorcell)
    ###replace cell type calls
    adata.obs[label] = adata.obs.apply(lambda row: replace_celltypes(celldict,row,label),axis=1)
    ###remove cells flagged as non confident calls (probable doublets)
    adata = adata[adata.obs[label]!='delete']
    ###return object
    return(adata)

def replace_celltypes(celldict,row,label):
    """ Replaces cell type with max celltype in anndata object, but preserves original call if no max cell type is found. """
    if row['leiden'] in celldict:
        return(celldict[row['leiden']])
    else:
        return(row[label])

#######################################################
####################################################### MixNSeq tools
#######################################################

def calc_proportions(adata,label='cell_type',treatment='treatment',vehicle='DMSO',deseq_norm=False):
    ##calculate proportion by group
    ####get lengths of each treatment dataset
    sizedict={}
    for drug in set(adata.obs[treatment]):
        onedf = adata[adata.obs[treatment]==drug].obs
        sizedict[drug] = len(onedf)
    ####deseq normalize
    if deseq_norm==True:
        results = pd.DataFrame()
        ##count cells by type
        for celltype in set(adata.obs[label]):
            ##get cluster alone first
            onedf = adata[adata.obs[label]==celltype].obs
            ##count
            props = pd.DataFrame(onedf.groupby([treatment]).count()['barcode'])
            ##append
            results[celltype] = props['barcode']
        ##create pseudo reference
        pseudo = gmean(results.iloc[:,:],axis=0)
        ##calculate ratios
        tmp = results/pseudo
        ratios = tmp.median(axis=1)
        ##divide
        results = results.transpose()/ratios.tolist()
        ###divide everything by the vehicle
        results = results.div(results[vehicle],axis=0)
        return(results)
    ####regular normalization
    else:
        ####calculate relative proportions by cell line to get fitness relative to vehicle
        ##save data to results
        results = pd.DataFrame()
        for celltype in set(adata.obs[label]):
            ##get cluster alone first
            onedf = adata[adata.obs[label]==celltype].obs
            ##calculate proportions
            props = pd.DataFrame(onedf.groupby([treatment]).count()['barcode'])
            props = pd.DataFrame(props.apply(lambda row: row['barcode']/sizedict[row.name],axis=1))
            ##write into results dataframe
            results[celltype] = props[0]
        ###divide everything by the vehicle
        results = results.transpose()
        results = results.div(results[vehicle],axis=0)
        ###scale everything to max
        results = results/results.max()
        ###
        return results

def plot_fitness(fitness):
    ###plot fitness scores
    palette = sns.color_palette("Purples_d",n_colors=int(len(fitness)*1.3))
    ###
    for treatment in set(fitness.columns):
        ###format
        subset = pd.DataFrame(fitness[treatment])
        subset.sort_values(by=treatment,inplace=True)
        subset['cell type'] = subset.index
        ###plot
        sns.set(style="white")
        plt.figure(figsize=(10,6))
        sns.barplot(x='cell type',y=treatment, data=subset,palette=palette)
        ##formatting
        plt.title('Fitness by Cell Line after '+treatment+' treatment')
        plt.xticks(rotation=90)
        plt.ylabel('Fitness')
        plt.xlabel('Cell Type')
        plt.tight_layout()
        plt.savefig('figures/'+treatment+'.pdf',dpi=200)

def get_PT_path(indata,timelabel,origin_group,dcN,min_search=True,fit_genes=True):
    """ Does the heavy lifting on pseudotime calculation after you define the axis you want to search on. Also fits linear model to predict genes along PT trajectory. """
    #####this is a particular DC. You mustknow the one you want prior.
    countsN = 50
    ###formatting
    adata = indata.copy()
    # adata.var.drop(labels=['n_cells'],inplace=True,axis=1)
    sc.pp.filter_genes(adata, min_counts=countsN)
    adata.X[adata.X<0] = 0
    print(adata.shape,'finding root cell ',datetime.now())
    #Find the T0 cell with the lowest/highest DC1 value to act as root for the diffusion pseudotime and compute DPT
    stem_mask = np.isin(adata.obs[timelabel], origin_group)
    if min_search == True:
        max_stem_id = np.argmin(adata.obsm['X_diffmap'][stem_mask,dcN])
    else:
        max_stem_id = np.argmax(adata.obsm['X_diffmap'][stem_mask,dcN])
    root_id = np.arange(len(stem_mask))[stem_mask][max_stem_id]
    adata.uns['iroot'] = root_id
    #Compute dpt
    print('tracing pseudotime ',datetime.now())
    sc.tl.dpt(adata)
    ####fit model
    x = pd.DataFrame(adata.X,index=adata.obs_names,columns=adata.var_names)
    y = adata.obs.dpt_pseudotime
    ####
    if fit_genes == True:
        print('fitting genes to trajectory ',datetime.now())
        from sklearn import linear_model
        reg = linear_model.RidgeCV(alphas=np.array([1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.e+0, 1.e+1,1.e+2, 1.e+3, 1.e+4, 1.e+5, 1.e+6]),
                                   cv=None, fit_intercept=True, normalize=True,
                                   scoring=None, store_cv_values=False)
        reg.fit(x,y)
        ###print metrics
        print(reg.score(x,y),reg.intercept_)
        ###save values to new gene list
        dp1genes = pd.DataFrame(reg.coef_,index=adata.var_names,columns=['glm_coef'])
        dp1genes.sort_values(by='glm_coef',ascending=False,inplace=True)
        ###see only high expressors
        dp1genes_hi = dp1genes[dp1genes.index.isin(adata.var.index.tolist())]
        ###
        genes = dp1genes_hi.head(n=5).index.tolist()
        ###visualize
        for gene in genes:
            sc.pl.diffmap(adata, components='1,2', color=gene)
    else:
        dp1genes = None
    #Visualize pseudotime over variables
    sc.pl.diffmap(adata, components='1,2', color='dpt_pseudotime')
    sc.pl.diffmap(adata, components='1,2', color=timelabel)
    sc.pl.diffmap(adata, components='1,2', color='phase')
    ###returns the pseudotime calculated df
    return(adata,dp1genes)

def vis_pseudostable(adata):
    """ Histogram across PT state to allow definition of pseudotime """
    # Density Plot and Histogram of all arrival delays
    plt.figure(figsize=(15, 15))
    sns.distplot(adata.obs.dpt_pseudotime, hist=True, kde=True,
                 bins=int(len(adata)/20), color = 'darkblue',
                 hist_kws={'edgecolor':'black'},
                 kde_kws={'linewidth': 4})
    plt.xlim(-0.1,1.1)
    plt.yscale('log',nonposy='clip')
    plt.ylim(0.09,10**2)


def vis_pseudostable_covariates(adata,metadata):
    """ Histogram across PT state to allow definition of pseudotime overlaid with the covariates you want to see. """
    plt.figure(figsize=(15, 15))
    for meta_i in range(len(metadata)):
        meta = metadata[meta_i]
        ###plot individual
        plt.figure()
        subpops = list(set(adata.obs[meta]))
        subpops.sort()
        for subpop in subpops:
            tmp =adata[adata.obs[meta]==subpop]
            ###plot
            sns.distplot(tmp.obs.dpt_pseudotime, hist = False, kde = True,
                     kde_kws = {'linewidth': 3},
                     label = subpop)
            plt.title(meta)
        plt.xlim(-0.1,1.1)
        plt.yscale('log',nonposy='clip')
        plt.ylim(0.09,10**2)

def categorize_pseudostable(adata,PT_label,states):
    """ Once pseudotimes have been calculated, takes a user defined list of cutoffs between 0 and 1 to categorize cells into pseudostates. """
    adata.obs['pseudostate'] = adata.obs.apply(lambda row: lambda_categorize(row[PT_label],states),axis=1)
    return(adata)

def lambda_categorize(PT,states):
    i = 1
    for cutoff in states:
        if PT <= cutoff:
             return(str(i))
        i+=1
