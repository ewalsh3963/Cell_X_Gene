#/usr/bin/python3
# Data Science 
# Evan Walsh (evan.walsh@capsida.com)

import os
import sys
import glob
import scanpy as sc
import numpy as np
import pandas as pd
import json
import matplotlib.gridspec as gridspec
import scipy
from matplotlib import pyplot as plt
import pdb
import kneed
import re
import OS_Tools



class loadAdata:
    def __init__(self, out):
        self.out = out
    
    def loadAdata(self):
        adataFile=os.path.join(self.out, 'AdataFiles','1_adata.h5ad')
        if not os.path.exists(adataFile):
            sys.exit("ERROR: " + str(adataFile) + " object does not exist...")
        else: 
            adata = sc.read_h5ad(adataFile)
        return adata, adataFile

class lowD: 
    def __init__(self, adata, nHVGs, theta, res, nComps, nNeigh):
        self.adata = adata
        self.nHVGs = nHVGs
        self.theta = theta
        self.res = res
        self.nComps = nComps
        self.nNeigh = nNeigh
    
    def recipe(self):
        self.adata.layers["raw"] = self.adata.X.copy()
        self.adata.layers["sqrt_norm"] = np.sqrt(sc.pp.normalize_total(self.adata, inplace=False)["X"])
        sc.experimental.pp.recipe_pearson_residuals(self.adata, n_top_genes = self.nHVGs, inplace = True, theta = self.theta)
        
        cmd = """ sc.experimental.pp.recipe_pearson_residuals(self.adata, n_top_genes = self.nHVGs, inplace = True, theta = self.theta) """
        
        return cmd
        # return self.adata, cmd

    def normalization(self):
        ## save raw and square root transformed matrices 
        self.adata.layers["raw"] = self.adata.X.copy()
        self.adata.layers["sqrt_norm"] = np.sqrt(sc.pp.normalize_total(self.adata, inplace=False)["X"])
        sc.experimental.pp.normalize_pearson_residuals(self.adata, inplace=True, theta=self.theta)
        cmd = """ sc.experimental.pp.normalize_pearson_residuals(self.adata, inplace=True, theta=self.theta) """
        return cmd
        # return self.adata, cmd

    def feature_selection(self, adata):
        sc.experimental.pp.highly_variable_genes(adata, flavor="pearson_residuals", n_top_genes=self.nHVGs, layer='raw', inplace=True)
        cmd = """ sc.experimental.pp.highly_variable_genes(adata, flavor="pearson_residuals", n_top_genes=self.nHVGs, layer='raw', inplace=True) """        
        return cmd

    def PCA(self, adata):
        sc.pp.pca(adata, n_comps=self.nComps)
        return adata 

    def UMAP(self, adata, n_pcs):
        sc.pp.neighbors(adata, n_neighbors=self.nNeigh, n_pcs=n_pcs)
        sc.tl.umap(adata)
        return adata

    def clustering(self, adata):
        for i in self.res:
            keyName = 'leiden_' + str(i)
            sc.tl.leiden(adata, resolution=i, key_added = keyName)
        return adata

class plotQC:
    def __init__(self, adata, nComps):
        self.adata = adata
        self.nComps = nComps
    
    def MeanVar(self, ax_object):
        ax_object.set_xlabel("Avg expression level")
        ax_object.set_ylabel("Residual Variance")
        
        # get discrete colormap
        cmap = plt.get_cmap('viridis', 2)
        scatter = ax_object.scatter(np.log(self.adata.var['means']), self.adata.var['residual_variances'], s=6, edgecolors='none', c=self.adata.var.highly_variable, cmap=cmap)

        return scatter

    def iqrVar(self, ax_object):
        ax_object.set_xlabel("Residual IQR")
        ax_object.set_ylabel("Residual Variance")
        
        cmap = plt.get_cmap('viridis', 2)
        scatter = ax_object.scatter(self.adata.var['scaled_iqr'], self.adata.var['residual_variances'], s=6, edgecolors='none', c=self.adata.var.highly_variable, cmap=cmap)
        
        return scatter

    def PCAcor(self, ax_object):
        ## measure the correlation between PC loadings and sequencing depth 
        pcIndex = ['PC_' + str(f) for f in np.arange(self.nComps)]
        pcLoadings = pd.DataFrame((self.adata.obsm['X_pca']), index=self.adata.obs.index, columns=pcIndex)

        PCs = pcLoadings.columns
        pcLoadings = pd.merge(pcLoadings, self.adata.obs['tscp_count'], left_index=True, right_index=True)

        ## pearson correlation
        pearValues = []
        spearValues = []
        for PC in PCs:
            pearValues.append(abs(scipy.stats.pearsonr(pcLoadings[PC], pcLoadings['tscp_count'])[0]))
            spearValues.append(abs(scipy.stats.spearmanr(pcLoadings[PC], pcLoadings['tscp_count'])[0]))

        pearValues = pd.DataFrame(pearValues, index=pcIndex, columns=['PearsonR'])
        spearValues = pd.DataFrame(spearValues, index=pcIndex, columns=['SpearmanR'])
        corValues = pd.merge(pearValues, spearValues, left_index=True, right_index=True)
        corValues = corValues[0:31] ## subset to first 30 PCs
        #corValues['PC'] = corValues.index

        PClist = [re.sub("^PC_", "", x) for x in corValues.index]
        PCdat = pd.DataFrame([float(i) for i in PClist], index=corValues.index, columns=['PC'])
        corValues = pd.merge(corValues, PCdat, left_index=True, right_index=True)

        # Width of a bar 
        width = 0.3

        # Plotting
        ax_object.bar(corValues['PC'], corValues['PearsonR'] , width, label='Pearson')
        ax_object.bar(corValues['PC'] + width, corValues['SpearmanR'], width, label='Spearman')
        
        ax_object.set_xlabel('Principal Component')
        ax_object.set_ylabel('| R |')
        return ax_object

    def PCAratio(self, ax_object):
        pcDat = pd.DataFrame(np.log2(self.adata.uns['pca']['variance_ratio']))
        # pcIndex = ['PC_' + str(f) for f in np.arange(self.nComps)]
        # pcDat['pcIndex'] = pcIndex
        
        
        pcDat.columns = ['VarianceRatio']
        pcDat = pcDat[0:31]
        pcDat['PC'] = np.arange(31)

        ax_object.scatter(pcDat['PC'], pcDat['VarianceRatio'], s=10)
        ax_object.set_xlabel('Prinicpal Component')
        ax_object.set_ylabel('Variance ratio (log2)')

        return ax_object

    def umap(self, key, title, ax_object):
        sc.pl.umap(self.adata, color=key, title=title, ax=ax_object, show=False)

def main(args, run_logs):
    ## Get arguments from argparse
    outroot, study, recipe, nHVGs, theta = args.o, args.s, args.recipe, args.nHVGs, args.theta
    pdfDir = os.path.join(outroot, "lowD"); OS_Tools.ensure_directory(pdfDir, critical=False) 
    ## initialize Log tools
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    ##################################################################################################
    ## load the adata object
    adataClass = loadAdata(outroot)
    adata, adataFile = adataClass.loadAdata()
    my_logs.append_logs(breaker)
    my_logs.append_logs("Starting the low-dimensional transformation analysis module")
    my_logs.append_logs("Loading adata object from file: " + adataFile + '\n')
    
    ##################################################################################################
    ## run the low dimensional transformation pipeline
    res, nComps, nNeigh = args.res, args.nComps, args.nNeigh
    res = res.split()
    res = [eval(i) for i in res]
    lowDClass = lowD(adata, nHVGs, theta, res, nComps, nNeigh)
    
    if recipe:
        my_logs.append_logs("Running full pipeline for HVG selection and normalization by analytic Pearson residuals:")
        cmd = lowDClass.recipe()
        my_logs.append_logs('\tRecipe command:' + cmd)
    else:
        ## 1) pearson residual normalization
        my_logs.append_logs("Performing pearson residual normalization...:")
        cmd = lowDClass.normalization()
        my_logs.append_logs("\tNormalization command:" + cmd)
        ## 2) feature selection
        my_logs.append_logs("\nPerforming high-variance feature selection...:")
        cmd = lowDClass.feature_selection(adata)
        my_logs.append_logs("\tFeature selection command:" + cmd)

        ## 3) PCA
        my_logs.append_logs("\nPerforming Principal Compnents Analysis...:")
        my_logs.append_logs("\tNumber of principal compnents calculated: " + str(nComps))
        lowDClass.PCA(adata)
    
    pcDat = pd.DataFrame(np.log2(adata.uns['pca']['variance_ratio']))
    pcDat.columns = ['VarianceRatio']
    pcDat = pcDat[0:31]
    pcDat['PC'] = np.arange(31)
    kneedle = kneed.KneeLocator(pcDat['PC'], pcDat['VarianceRatio'], S=1.0, curve="convex", direction="decreasing")
    pc_lim = kneedle.knee

    pc_lim_dat = {"PC_threshold": str(pc_lim)}
    outfile = os.path.join(pdfDir, 'pc_lim.json')
    with open(outfile, 'w') as f:
        json.dump(pc_lim_dat, f, ensure_ascii=False)
    
    ## 4) compute neighborhood graph for UMAP
    my_logs.append_logs("\nPerforming UMAP manifold embedding and cell clustering...:\n")
    if args.pc_lim:
        lowDClass.UMAP(adata, n_pcs=pc_lim) ## UMAP module includes calculating the neighborhood and UMAP graph
    else:
        lowDClass.UMAP(adata, n_pcs=30) ## UMAP module includes calculating the neighborhood and UMAP graph


    ## 6) run cell clustering 
    lowDClass.clustering(adata)
    clusterDF = pd.DataFrame(adata.obs.loc[:, adata.obs.columns.str.startswith('leiden')].apply(np.unique, axis=0).apply(len), columns=['nClusters']).reset_index(inplace=False).rename(columns = {'index':'resolution'})
    my_logs.append_logs('Number of clusters generated:')
    my_logs.append_logs(clusterDF)
    my_logs.append_logs('\n')

    ## calculate IQR and variance of scaled values for plotting functions
    iqr = np.apply_along_axis(scipy.stats.iqr, 0, adata.X)
    var = np.apply_along_axis(np.var, 0, adata.X)
    adata.var = pd.merge(adata.var, pd.DataFrame(iqr, index=adata.var.index, columns=['scaled_iqr']), left_index=True, right_index=True)
    adata.var = pd.merge(adata.var, pd.DataFrame(var, index=adata.var.index, columns=['scaled_variation']), left_index=True, right_index=True)

    ##################################################################################################
    ## save the updated adata object
    if args.overwrite:
        my_logs.append_logs('Overwriting the adata object...')
        my_logs.append_logs('Saving adata object as : ' + adataFile)
        my_logs.append_logs(breaker)
        adata.write(adataFile)
    else:
        ## list adata files
        adataDir = os.path.join(outroot, 'AdataFiles')
        adataFiles = sorted(glob.glob(adataDir + '/[0-9]*.h5ad'))
        numFiles = len(adataFiles)
        fileID = numFiles + 1
        ## save the adata file
        outfile = os.path.join(adataDir, str(fileID) + '_adata.h5ad')
        my_logs.append_logs('Saving adata object as : ' + outfile)
        my_logs.append_logs(breaker)
        adata.write(outfile)

if __name__ == "__main__":
    main()


    
