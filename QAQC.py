# !/bin/python3
# Evan Walsh (evan.walsh@capsida.com)

import os
import sys
import re
import scanpy as sc
import csv
import numpy as np
from threading import Thread
import pandas as pd
from pybiomart import Server
import scvi 
from matplotlib.backends.backend_pdf import PdfPages
import glob
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib_venn import venn2
from scipy.stats import fisher_exact
from matplotlib_venn import venn2_unweighted
from scipy.signal import find_peaks, peak_prominences
import seaborn as sns
import pdb
import scrublet as scr
import kneed
import math
import OS_Tools

"""
## QAQC Steps:
##  1) load DGE into an adata object
##  2) run scrublet (for doublet detection)
##  3) filter doublets and cells that do not make sequencing depth requirement (quantiative outliers)
"""

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

class plotQC:
    def __init__(self, adata):
        self.adata = adata
    
    def violin(self, key, ax_object):
        sc.pl.violin(self.adata, keys=key, groupby='sample', show = False, rotation=90, ax=ax_object)
        return ax_object

    def scatter(self, feature_x, feature_y, label_x, label_y, ax_object, thresholds = None):
        
        ## set axis labels
        ax_object.set_xlabel(label_x)
        ax_object.set_ylabel(label_y)        
        ax_object.tick_params(axis='x', labelrotation=90)
        
        if thresholds is not None:
            self.adata.obs['Flag'] = 0
            self.adata.obs.loc[(self.adata.obs.tscp_count > thresholds[3]) | (self.adata.obs.gene_count > thresholds[1]), 'Flag'] = 1
            self.adata.obs.loc[(self.adata.obs.tscp_count < thresholds[2]) | (self.adata.obs.gene_count < thresholds[0]), 'Flag'] = 1
            
            
            ax_object.axhline(y=thresholds[0], xmin=0.0, xmax=1.0, color='grey', linestyle='--')
            ax_object.axvline(x=thresholds[2], ymin=0.0, ymax=1.0, color='grey', linestyle='--')
            ax_object.axhline(y=thresholds[1], xmin=0.0, xmax=1.0, color='grey', linestyle='--')
            ax_object.axvline(x=thresholds[3], ymin=0.0, ymax=1.0, color='grey', linestyle='--')

            # get discrete colormap
            cmap = plt.get_cmap('viridis', 2)
            scatter = ax_object.scatter(self.adata.obs[feature_x], self.adata.obs[feature_y], s=10, edgecolors='none', c=self.adata.obs["Flag"], cmap=cmap)

        else:
            scatter = ax_object.scatter(self.adata.obs[feature_x], self.adata.obs[feature_y], s=10, edgecolors='none')
        
        return scatter
    
    def venn_diagram(tscpDat, adata):
        tg_barcodes = set(tscpDat['barcode'].unique().tolist())
        wt_barcodes = set(adata.obs.index.to_list())
        # wt_barcodes = set([re.sub('__.*','',x) for x in wt_barcodes])

        pre_val1 = len(tg_barcodes - wt_barcodes)
        pre_val2 = len(wt_barcodes - tg_barcodes)
        pre_val3 = len(set(tg_barcodes) & set(wt_barcodes))

        return [pre_val1, pre_val2, pre_val3]

class scrublet:
    def __init__(self, adata, args):
        self.adata = adata
        self.args = args
    
    def run_scrublet(self):
        exp_doub_rate = self.args.exp_doub_rate
        sim_ratio = self.args.sim_ratio
        knn_metric = self.args.knn_metric
        
        sc.external.pp.scrublet(adata = self.adata, sim_doublet_ratio = sim_ratio, expected_doublet_rate = exp_doub_rate, knn_dist_metric = knn_metric, mean_center = True, threshold  = None)
        
        cmd = """ sc.external.pp.scrublet(adata = self.adata, sim_doublet_ratio = sim_ratio, expected_doublet_rate = exp_doub_rate, knn_dist_metric = knn_metric, mean_center = True, threshold  = None) """

        ## print histogram of doublet scores 
        singletonScores = self.adata.obs[self.adata.obs['predicted_doublet'] == False].doublet_score.to_numpy()
        obsScores = self.adata.obs[self.adata.obs['predicted_doublet'] == True].doublet_score.to_numpy()
        simScores = self.adata.uns['scrublet']['doublet_scores_sim']

        ## initialize plot 
        fig, ax = plt.subplots()

        bins = np.linspace(start = 0, stop = 1, num=100)
        ax.hist(singletonScores, bins = bins, stacked=True, density=True, histtype='barstacked')
        ax.hist(obsScores, bins = bins, stacked=True, density=True, histtype='barstacked')
        ax.hist(simScores, bins = bins, stacked=True, density=True, histtype='barstacked')

        ## create legend and labels
        labels = ["singleton","obsDoublet", "simDoublet"]
        ax.legend(labels, loc='upper left')
        ax.set_xlabel('Doublet Score')
        ax.set_ylabel('Density')
        ax.set_title('Scrublet Score Distribution')

        return self.adata, fig, cmd

def main(adata, args, id, organism, my_logs):
    # Get arguments from argparse
    outroot, study = args.o, args.s
    pdfDir = os.path.join(outroot, "QAQC"); OS_Tools.ensure_directory(pdfDir, critical=False) 

    ## initialize Log tools and directories
    # my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120
    my_logs.append_logs(breaker)
    my_logs.append_logs("Running QAQC script")

    ##################################################################################################
    ## load the adata object
    # adataClass = loadAdata(outroot)
    # adata, adataFile = adataClass.loadAdata()
    my_logs.append_logs(breaker)
    my_logs.append_logs("Starting the QAQC analysis module")
    # my_logs.append_logs("Loading adata object from file: " + adataFile + '\n')
    
    ##################################################################################################
    ## Calculate QC metrics including % UMIs mapping to MT genes
    adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
    unfilt_dat = adata.obs.groupby('donor_id').size().reset_index()
    unfilt_dat.columns = ['donor_id','total']

    var_names = adata.var_names.to_list() 
    if organism != 'Homo sapiens':
        human_mt_ensembl = [
            "ENSG00000198695",
            "ENSG00000198712",
            "ENSG00000198727",
            "ENSG00000198763",
            "ENSG00000198786",
            "ENSG00000198804",
            "ENSG00000198840",
            "ENSG00000198886",
            "ENSG00000198888",
            "ENSG00000198899",
            "ENSG00000198938",
            "ENSG00000212907",
            "ENSG00000228253"
        ]

        ## define attributes to recover based on species of interest
        ortho_species = organism
        ortho_species_biomart = (ortho_species[0] + re.sub('.* ','',ortho_species)).lower()
        att1 = ortho_species_biomart + '_homolog_ensembl_gene'   
        att2 = ortho_species_biomart + '_homolog_associated_gene_name'
        
        ## define the mart/serve to use for querying
        server = Server(host='http://www.ensembl.org') #server.list_marts()
        mart = server['ENSEMBL_MART_ENSEMBL'] # mart.list_datasets()
        dataset = mart['hsapiens_gene_ensembl'] 
        dat = dataset.query(attributes=['ensembl_gene_id', att1, att2], filters={'link_ensembl_gene_id': human_mt_ensembl})
        dat = dat.dropna()
        mt_genes = dat[dat.columns[2]].to_list()
        adata.var['mt'] = [x in mt_genes for x in var_names]
    else:
        mt_genes = [x for x in adata.var['feature_name'].to_list() if x.startswith('MT-')]
        adata.var['mt'] = adata.var['feature_name'].isin(mt_genes)
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    ##################################################################################################
    ## Remove cell barcodes at the low-end of the transcript threshold (> min_UMIs, min_genes)
    min_genes, max_genes, min_transcripts, max_transcripts = args.min_genes, args.max_genes, args.min_transcripts, args.max_transcripts
    min_genes = math.floor(np.quantile(adata.obs.n_genes_by_counts, float(min_genes)))
    min_transcripts = math.floor(np.quantile(adata.obs.raw_sum, float(min_transcripts)))
    max_genes = math.floor(np.quantile(adata.obs.nnz, float(max_genes)))
    max_transcripts = math.floor(np.quantile(adata.obs.raw_sum, float(max_transcripts)))

    ## Calculate the total num of unfiltered cells per sample
    filter_dat = adata.obs.groupby('donor_id').size().reset_index()
    filter_dat.columns = ['donor_id','Unfiltered Num of Cells']

    ## Remove cells based on minimum num of genes threshold
    sc.pp.filter_cells(adata, min_genes=min_genes);
    filter_tmp_1 = adata.obs.groupby('donor_id').size().reset_index()
    filter_tmp_1.columns = ['donor_id','Filtered Num of Cells']
    filter_tmp_1['Statistic'] = "Min num of genes"
    filter_tmp_1['Value'] = min_genes
    filter_tmp_1 = pd.merge(filter_dat, filter_tmp_1)
    filter_tmp_1['Num of cells removed'] = filter_tmp_1['Unfiltered Num of Cells'] - filter_tmp_1['Filtered Num of Cells']
    filter_tmp_1 = filter_tmp_1.drop(["Unfiltered Num of Cells", "Filtered Num of Cells"], axis = 1)
    filter_tmp_1 = filter_tmp_1[['donor_id','Statistic','Value','Num of cells removed']]

    ## Reset Total
    filter_dat = adata.obs.groupby('donor_id').size().reset_index()
    filter_dat.columns = ['donor_id','Unfiltered Num of Cells']

    ## Remove cells based on minimum num of UMIs threshold
    sc.pp.filter_cells(adata, min_counts=min_transcripts);
    filter_tmp_2 = adata.obs.groupby('donor_id').size().reset_index()
    filter_tmp_2.columns = ['donor_id','Filtered Num of Cells']
    filter_tmp_2['Statistic'] = "Min num of transcripts"
    filter_tmp_2['Value'] = min_transcripts
    filter_tmp_2 = pd.merge(filter_dat, filter_tmp_2)
    filter_tmp_2['Num of cells removed'] = filter_tmp_2['Unfiltered Num of Cells'] - filter_tmp_2['Filtered Num of Cells']
    filter_tmp_2 = filter_tmp_2.drop(["Unfiltered Num of Cells", "Filtered Num of Cells"], axis = 1)
    filter_tmp_2 = filter_tmp_2[['donor_id','Statistic','Value','Num of cells removed']]

    ## Reset Total
    filter_dat = adata.obs.groupby('donor_id').size().reset_index()
    filter_dat.columns = ['donor_id','Unfiltered Num of Cells']

    ##################################################################################################
    ## add files used to construct adata object to log file 
    breaker = "=" * 120
    my_logs.append_logs(breaker)
    my_logs.append_logs("Starting the Quality Analysis and Quantily control module...:\n")
    my_logs.append_logs("Creating an adata object with scanpy for sc-analysis...:\n")
    my_logs.append_logs('Num of unfiltered cells:' + str(adata.shape[0]))
    my_logs.append_logs('Num of unfiltered genes:' + str(adata.shape[1]))

    ## remove low abundant genes
    sc.pp.filter_genes(adata, min_cells=3)

    ## run scrublet for doublet detection
    if args.scrublet and not args.solo:
        scrub = scr.Scrublet(adata.X)
        scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

        if args.scrt:
            scrub.call_doublets(threshold=float(args.scrt))
        
        umap_dat = pd.DataFrame(scrub._embeddings['UMAP'])
        umap_dat.columns = ['UMAP_1','UMAP_2']
        scores_dat = pd.DataFrame({'doublet_score': scrub.doublet_scores_obs_, 'predicted_doublet': scrub.predicted_doublets_})
        dat = pd.concat([umap_dat, scores_dat], axis = 1)
        dat = pd.concat([adata.obs.reset_index()[['raw_sum','nnz']], dat], axis = 1)

        ## Save doublet data
        outfile=os.path.join(pdfDir, 'doubletData.csv')
        dat.to_csv(outfile, sep=',', header = 1, index = False)
        
        simDat = pd.DataFrame({"simulated_doublet":scrub.doublet_scores_sim_})
        outfile=os.path.join(pdfDir, 'sim_doubletDat.csv')
        simDat.to_csv(outfile, sep=',', header = 1, index = False)

        ax2 = sns.kdeplot(data=simDat["simulated_doublet"], bw_adjust = 1)
        x = ax2.lines[0].get_data()[0]
        y = ax2.lines[0].get_data()[1]

        # compute second derivative
        smooth_d2 = np.gradient(np.gradient(y))
        infls = np.where(np.diff(np.sign(smooth_d2)))[0]
        new_x = x[np.where((x >= x[infls[1]]))] 
        new_y = y[np.where((x >= x[infls[1]]))] 
        kneedle = kneed.KneeLocator(new_x, new_y, S=1.0, curve="convex", direction="decreasing")

        plt.close()

        tab = pd.DataFrame([scrub.expected_doublet_rate, scrub.overall_doublet_rate_, scrub.threshold_, kneedle.knee, scrub.detectable_doublet_fraction_ ])
        tab.columns = ['value']
        tab['item'] = ['Expected Doublet Rate', "Observed doublet rate", "Score Threshold", "Density Threshold", "Detectable Doublet Fraction"]
        tab = tab[['item','value']]
        outfile=os.path.join(pdfDir, 'summary_doubletData.csv')
        tab.to_csv(outfile, sep=',', header = 1, index = False)

        if args.scrt:
            scrub.call_doublets(threshold=float(args.scrt))
        else:
            scrub.call_doublets(threshold=float(kneedle.knee))

        scrub.plot_histogram();
        plt.savefig(os.path.join(pdfDir,'doublet_score_histogram.png'))
        plt.close()
        adata.obs['predicted_doublet'] =  scrub.predicted_doublets_
    
    elif args.solo:
        nCells = adata.obs.shape[0]
        if nCells <= 2000:
            nEpochs = 400
        elif nCells <= 5000:
            nEpochs = 200
        elif nCells <= 10000:
            nEpochs = 100
        elif nCells <= 10000:
            nEpochs = 75
        elif nCells > 50000:
            nEpochs = 25
        scvi.model.SCVI.setup_anndata(adata)
        vae = scvi.model.SCVI(adata)
        vae.train(max_epochs = nEpochs, early_stopping=True, early_stopping_patience=20, early_stopping_monitor="elbo_validation")
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train(max_epochs = nEpochs, early_stopping=True, early_stopping_patience=20, early_stopping_monitor="elbo_validation")

        scores_dat = solo.predict()
        scores_dat['prediction'] = solo.predict(soft = False)
        scores_dat.index = scores_dat.index.map(lambda x: x[:-2])

        n_samp = int(adata.obs.shape[0] * 0.4)
        barcodes = adata.obs.sample(n_samp).index.to_list()
        adata_tmp = adata[barcodes].copy()
        sc.pp.normalize_total(adata_tmp, target_sum = 1e4)
        sc.pp.log1p(adata_tmp)
        sc.tl.pca(adata_tmp)
        sc.pp.neighbors(adata_tmp)
        sc.tl.umap(adata_tmp)
        # sc.tl.leiden(adata_tmp, resolution = 0.5)

        umap_dat = pd.DataFrame(adata_tmp.obsm['X_umap'])
        umap_dat.index = adata_tmp.obs.index
        umap_dat.columns = ['UMAP_1','UMAP_2']
        umap_dat = pd.merge(umap_dat, scores_dat, left_index = True, right_index = True)
        dat = pd.merge(umap_dat, adata_tmp.obs[['raw_sum','nnz']], left_index = True, right_index = True)
        outfile=os.path.join(pdfDir, 'doubletData.csv')
        dat.to_csv(outfile, sep=',', header = 1, index = False)

        overall_doublet_rate = round(100* (sum(scores_dat['prediction'] == 'doublet') / scores_dat.shape[0]), 2)
        tab = pd.DataFrame([None, overall_doublet_rate, None, None, None])
        tab.columns = ['value']

        tab['item'] = ['Expected Doublet Rate', "Observed doublet rate", "Score Threshold", "Density Threshold", "Detectable Doublet Fraction"]
        tab = tab[['item','value']]
        outfile=os.path.join(pdfDir, 'summary_doubletData.csv')
        tab.to_csv(outfile, sep=',', header = 1, index = False)
        adata.obs['predicted_doublet'] =  np.where(scores_dat.prediction == 'doublet', True, False)

    else:
        adata.obs['predicted_doublet'] = False
    
    adata.obs['predicted_doublet_inv'] = adata.obs['predicted_doublet'].apply(lambda x: not x)
    filter_tmp_3 = adata.obs.groupby('donor_id')['predicted_doublet_inv'].sum().reset_index() 
    filter_tmp_3.columns = ['donor_id','Filtered Num of Cells']
    filter_tmp_3['Statistic'] = "Doublet Threshold"
    filter_tmp_3['Value'] = min_genes
    try:
        filter_tmp_3['Value'] = scrub.threshold_
    except NameError:
        filter_tmp_3['Value'] = None

    filter_tmp_3 = pd.merge(filter_dat, filter_tmp_3)
    filter_tmp_3['Num of cells removed'] = filter_tmp_3['Unfiltered Num of Cells'] - filter_tmp_3['Filtered Num of Cells']
    filter_tmp_3 = filter_tmp_3.drop(["Unfiltered Num of Cells", "Filtered Num of Cells"], axis = 1)
    filter_tmp_3 = filter_tmp_3[['donor_id','Statistic','Value','Num of cells removed']]

    ## Reset Total
    filter_dat = adata.obs.groupby('donor_id').size().reset_index()
    filter_dat.columns = ['donor_id','Unfiltered Num of Cells']

    ##################################################################################################
    ## remove the cells that do not meet gene and transcript max cut offs

    ## Remove cells with high num of genes
    sc.pp.filter_cells(adata, max_genes=max_genes);
    filter_tmp_4 = adata.obs.groupby('donor_id').size().reset_index()
    filter_tmp_4.columns = ['donor_id','Filtered Num of Cells']
    filter_tmp_4['Statistic'] = "Max num of genes"
    filter_tmp_4['Value'] = max_genes
    filter_tmp_4 = pd.merge(filter_dat, filter_tmp_4)
    filter_tmp_4['Num of cells removed'] = filter_tmp_4['Unfiltered Num of Cells'] - filter_tmp_4['Filtered Num of Cells']
    filter_tmp_4 = filter_tmp_4.drop(["Unfiltered Num of Cells", "Filtered Num of Cells"], axis = 1)
    filter_tmp_4 = filter_tmp_4[['donor_id','Statistic','Value','Num of cells removed']]

    ## Reset Total
    filter_dat = adata.obs.groupby('donor_id').size().reset_index()
    filter_dat.columns = ['donor_id','Unfiltered Num of Cells']

    ## Remove cells with high num of UMIs
    sc.pp.filter_cells(adata, max_counts=max_transcripts);
    filter_tmp_5 = adata.obs.groupby('donor_id').size().reset_index()
    filter_tmp_5.columns = ['donor_id','Filtered Num of Cells']
    filter_tmp_5['Statistic'] = "Max num of transcripts"
    filter_tmp_5['Value'] = max_transcripts
    filter_tmp_5 = pd.merge(filter_dat, filter_tmp_5)
    filter_tmp_5['Num of cells removed'] = filter_tmp_5['Unfiltered Num of Cells'] - filter_tmp_5['Filtered Num of Cells']
    filter_tmp_5 = filter_tmp_5.drop(["Unfiltered Num of Cells", "Filtered Num of Cells"], axis = 1)
    filter_tmp_5 = filter_tmp_5[['donor_id','Statistic','Value','Num of cells removed']]
    
    ## Remove Doublets
    barcodes_to_keep = adata.obs[adata.obs['predicted_doublet'] == False].index
    adata = adata[barcodes_to_keep]

    filterDat = pd.concat([filter_tmp_1, filter_tmp_2 , filter_tmp_3, filter_tmp_4, filter_tmp_5])
    filterDat = pd.merge(filterDat, unfilt_dat)
    filterDat['Percent of cells removed'] = 100 * (filterDat['Num of cells removed'] / filterDat['total'])
    filterDat = filterDat.drop('total', axis = 1)
    outfile = os.path.join(pdfDir, 'summary_filterData.csv')
    filterDat.to_csv(outfile, sep=',', header = True, index = False)

    outfile = os.path.join(pdfDir, 'unfilteredData.csv')
    unfilt_dat.to_csv(outfile, sep=',', header = True, index = False)

    ## summarize thresholds in the log file
    my_logs.append_logs("Adata object generation...:")
    my_logs.append_logs("Filtering thresholds...:")
    my_logs.append_logs('\tNum of genes (min, max) ' + str(min_genes) + ' ' + str(max_genes))
    my_logs.append_logs('\tNum of tscps (min, max) ' + str(min_transcripts) + ' ' + str(max_transcripts) + '\n')
        
    try: 
        adata.obs = adata.obs.drop(['barcode', 'Sublib'], axis = 1)
    except:
        pass

    ##################################################################################################
    ## save QC data and adata file
    sc.pp.filter_genes(adata, min_cells=1, inplace = True)

    ## QAQC 
    return adata
    ## save the adata file
    # adataDir = os.path.join(outroot, 'AdataFiles')
    # adataFile = os.path.join(adataDir, '1_adata.h5ad')
    # if args.overwrite:
    #     my_logs.append_logs('Overwriting the adata object...')
    #     my_logs.append_logs('Saving adata object as : ' + adataFile)
    #     my_logs.append_logs(breaker)
    #     adata.write(adataFile)
    # else:
    #     ## list adata files
    #     adataDir = os.path.join(adataDir, 'AdataFiles')
    #     adataFiles = sorted(glob.glob(adataDir + '/[0-9]*.h5ad'))
    #     numFiles = len(adataFiles)
    #     fileID = numFiles + 1
    #     ## save the adata file
    #     outfile = os.path.join(adataDir, str(fileID) + '_adata.h5ad')
    #     my_logs.append_logs('Saving adata object as : ' + outfile)
    #     my_logs.append_logs(breaker)
    #     adata.write(outfile)
    
if __name__ == "__main__":
    main()

