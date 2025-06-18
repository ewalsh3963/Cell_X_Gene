# !/bin/python3
# Evan Walsh (evan.walsh@capsida.com)

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import math
import seaborn as sns
import re
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gseapy
import glob
import subprocess
from pyensembl import EnsemblRelease
from collections import Counter
from pybiomart import Server
from pybiomart import Dataset
import pdb
import OS_Tools

"""
## GSEA Steps:
##  1) load marker gene file 
##  2) run differential expression
##  3) run GSEA
##  4) annotate cell types
"""

class loadAdata:
    def __init__(self, out):
        self.out = out
    def loadAdata(self):
        ## change the working directory to the adataFiles dir 
        adataDir = os.path.join(self.out, 'AdataFiles')
        os.chdir(adataDir)
        
        ## list adata files
        adataFiles = sorted(glob.glob('[0-9]*.h5ad'))
        adataFile = os.path.join(adataDir, adataFiles[len(adataFiles)-1])
        
        # adataFile=os.path.join(self.out, self.title, 'AdataFiles','2_adata.h5ad')
        if not os.path.exists(adataFile):
            sys.exit("ERROR: " + str(adataFile) + " object does not exist...")
        else: 
            adata = sc.read_h5ad(adataFile)
        return adata, adataFile

class loadMarkerGeneDB:
    def __init__(self, outdir, db_species, tissues):
        self.outdir = outdir
        self.db_species = db_species
        self.tissues = tissues
        self.dbDict = {'All': 
                "/home/ewalsh/Cell_X_Gene/CellMarkers/Cell_marker_All.xlsx",
                'Human': 
                "/home/ewalsh/Cell_X_Gene/CellMarkers/Cell_marker_Human.xlsx",
                'Mouse':
                "/home/ewalsh/Cell_X_Gene/CellMarkers/Cell_marker_Mouse.xlsx",
                'SingleCell':
                "/home/ewalsh/Cell_X_Gene/CellMarkers/Cell_marker_Seq.xlsx"
                }
        self.PangloDB_url = "/home/ewalsh/Cell_X_Gene/CellMarkers/PanglaoDB_markers_27_Mar_2020.tsv"
    
    def getDatabase(self):
        return self.dbDict[self.db_species]
    
    def loadDatabase(self, file):
            try:
                dat = pd.read_excel(file, header = 0, index_col = None)
            except:
                raise Exception('Marker gene database file not found...exiting')
            return dat
    def PangloDB(self, file):
        try:
            dat = pd.read_csv(file, sep='\t', header = 0, index_col = None)
        except:
            raise Exception('Marker gene database file not found...exiting')
        return dat

class DifferentialExpression:
    def __init__(self, adata, geneDict, numDEGs):
        self.adata = adata
        self.numDEGs = numDEGs
        self.geneDict = geneDict

    def runDE(self, groupVar, layer):
        sc.tl.rank_genes_groups(self.adata, groupby=groupVar, corr_method = 'benjamini-hochberg', method = 'wilcoxon', layer = layer,  pts = True)

    def GSEA(self, group):
        ## get the genes upregulated in each cluster and cross-reference them against marker gene database
        try:
            glist = sc.get.rank_genes_groups_df(self.adata, group=group, pval_cutoff = 0.05, log2fc_min = 0)['names'].squeeze().str.strip().tolist()
        except: 
            return
        try: 
            enr_res = gseapy.enrichr(gene_list=glist[:self.numDEGs], 
                                     gene_sets=self.geneDict, 
                                     background = self.adata.X.shape[1] - len(glist), 
                                     cutoff=1, 
                                     organism='human')
        except ValueError:
            dat = None
            return dat
        ## Return the gene set enrichment output
        if enr_res.results.shape[0] == 0:
            dat = None
        else:
            dat = enr_res.results
            dat['score'] = -np.log(dat["Adjusted P-value"]) * dat['Odds Ratio']
            dat = dat.sort_values('score', ascending=False)
            dat['group'] = group
        return dat

class plotAnnotatation:
    def __init__(self, adata, termDat):
        self.adata = adata
        self.termDat = termDat

    def gseaEnrichment(self, cluster, groupVar):
        self.termDat['log_pval_adj'] = -1*np.log10(self.termDat['Adjusted P-value'])
        dat = self.termDat.loc[(self.termDat['group'] == cluster) & (self.termDat['resolution'] == groupVar)]
        dat = dat.sort_values('log_pval_adj', ascending=False)
        dat = dat[dat['Adjusted P-value'] < 0.05]

        plt.figure(layout='constrained')
        plt.bar(dat.Term, dat.log_pval_adj)
        plt.xlabel("Cell type")
        plt.ylabel("-log10(adj. p-value)")
        plt.title('Resolution: ' + groupVar + '\nCell cluster: ' + cluster)
        plt.xticks(rotation=45, ha='right')
        plt.axhline(y=-1*np.log10(0.05), color='grey', linestyle='--')

        return plt
    
    def gseaEnrichment_score(self, cluster, groupVar):
        # self.termDat['log_pval_adj'] = -1*np.log10(self.termDat['Adjusted P-value'])
        dat = self.termDat.loc[(self.termDat['group'] == cluster) & (self.termDat['resolution'] == groupVar)]
        dat = dat.sort_values('score', ascending=False)

        plt.figure(layout='constrained')
        plt.bar(dat.Term, dat.score)
        plt.xlabel("Cell type")
        plt.ylabel("-log10(adj. p-value) * Odds Ratio")
        plt.title('Resolution: ' + groupVar + '\nCell cluster: ' + cluster)
        plt.xticks(rotation=90, ha='right')
        return plt

    def umap(self, key, title, ax_object, fontsize= 14, fontoutline = 4, legend='right margin', layer=None, add_outline = True):
        sc.pl.umap(self.adata, color=key, title=title, add_outline=add_outline, legend_loc=legend, legend_fontsize=fontsize, legend_fontoutline=fontoutline, frameon=False, ax=ax_object, show=False, layer=layer)

class findOrthos:
    def __init__(self, database, db_species, ortho_species):
        self.database = database
        self.db_species = db_species
        self.ortho_species = ortho_species
    
    def symbol_to_ensembl(self):
        speciesDict = {'Mouse': 'mus_musculus', 'Human': 'homo_sapiens'}
        if self.db_species == 'Mouse':
            data = EnsemblRelease(release = 108, species ='mus_musculus')
        else:
            data = EnsemblRelease(release = 108, species ='homo_sapiens')

        ## check if the EnsemblRelease data has been downloaded 
        dbList = subprocess.run(["pyensembl", "list"], shell = False)
        
        if re.search('108', dbList) and re.search(speciesDict[self.db_species], dbList):
            print('Ensembl database found...')
        else: 
            print("Ensembl database NOT found...downloading and indexing now")
            data.download()
            data.index()
        
        ## initialize ensembl ID column
        self.database['ensembl_id'] = None
        for i in np.arange(self.database.shape[0]):
            gene_symbol = self.database.Symbol[i]
            if isinstance(gene_symbol, str):
                try: 
                    self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = data.gene_ids_of_gene_name(gene_symbol)[0]                    
                except ValueError:
                    self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = None
            else:
                self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = None
        
        return self.database

    def symbol_to_ensembl_panglo(self):
        speciesDict = {'Mouse': 'mus_musculus', 'Human': 'homo_sapiens'}
        if self.db_species == 'Mouse':
            data = EnsemblRelease(release = 108, species ='mus_musculus')
        else:
            data = EnsemblRelease(release = 108, species ='homo_sapiens')

        ## check if the EnsemblRelease data has been downloaded 
        dbList = subprocess.run(["pyensembl", "list"], shell=False)
        
        if re.search('108', dbList) and re.search(speciesDict[self.db_species], dbList):
            print('Ensembl database found...')
        else: 
            print("Ensembl database NOT found...downloading and indexing now")
            data.download()
            data.index()
        
        ## filter marker gene database 
        if self.db_species == "Mouse":
            self.database = self.database[self.database['species'].str.contains("Mm")]
            self.database['official gene symbol'] = self.database['official gene symbol'].apply(lambda x: x[0] + x[1:].lower())
            self.database = self.database.rename(columns = {'official gene symbol':'Symbol'})
            self.database = self.database.reset_index().drop('index', axis=1)
        else:
            self.database = self.database.rename(columns = {'official gene symbol':'Symbol'})
            self.database = self.database.reset_index().drop('index', axis=1)

        ## initialize ensembl ID column
        self.database['ensembl_id'] = None
        for i in np.arange(self.database.shape[0]):
            gene_symbol = self.database.Symbol[i]
            if isinstance(gene_symbol, str):
                try: 
                    self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = data.gene_ids_of_gene_name(gene_symbol)[0]                    
                except ValueError:
                    self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = None
            else:
                self.database.loc[self.database['Symbol']==gene_symbol, 'ensembl_id'] = None        
        return self.database

    def queryBioMart(self, database):
        ## define attributes to recover based on species of interest
        att1 = self.ortho_species + '_homolog_ensembl_gene'   
        att2 = self.ortho_species + '_homolog_associated_gene_name'
        
        ## define the mart/serve to use for querying
        server = Server(host='http://www.ensembl.org', use_cache=False) #server.list_marts()
        mart = server['ENSEMBL_MART_ENSEMBL'] # mart.list_datasets()
        
        if self.db_species == 'Mouse':
            dataset = mart['mmusculus_gene_ensembl']
        else:
            dataset = mart['hsapiens_gene_ensembl']

        queryList = database['ensembl_id'].tolist()
        queryList = [i for i in queryList if i is not None] ## remove values that are None
        
        ## split the query list into chunks of size querySize
        querySize = 100; queryDF = pd.DataFrame()
        queryList = [queryList[i:i + querySize] for i in range(0, len(queryList), querySize)]
        queryFinal = []
        for i in range(0, len(queryList), querySize):
            queryFinal.append(queryList[i:i+querySize])        
        queryFinal = queryFinal[0]
        for item in queryFinal:
            ## run biomart query
            tmp = dataset.query(attributes=['ensembl_gene_id', att1, att2], filters={'link_ensembl_gene_id': item})
            try: 
                queryDF = pd.concat([queryDF, tmp])
            except ValueError:
                pass
        
        queryDF.columns = ['ref_ensembl','ortho_ensembl','ortho_symbol']
        return queryDF
        
def main(args, run_logs):
    ## Get arguments from argparse
    outroot, study, markerFile, species, db_species, tissues, find_ortho, numDEGs, PangloDB = args.o, args.s, args.markerFile, args.species, args.db_species, args.tissues, args.find_ortho, args.numDEGs, args.PangloDB
    pdfDir = os.path.join(outroot, "GSEA"); OS_Tools.ensure_directory(pdfDir, critical=False) 
    ## initialize Log tools
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    ## load the adata object with Scanpy
    workingDir = os.getcwd()
    adataClass = loadAdata(outroot)
    adata, adataFile = adataClass.loadAdata()
    os.chdir(workingDir)

    my_logs.append_logs(breaker)
    my_logs.append_logs("Starting the GSEA analysis module")
    my_logs.append_logs("Loading adata object from file: " + adataFile + '\n')

    ##################################################################################################
    ## load the marker gene database file 
    
    if not PangloDB: ## use normal marker gene database opposed to Panglo database
        if not markerFile: ## if a marker gene file is not inputted download the file from S3 
            markerGene = loadMarkerGeneDB(outroot, db_species, tissues)
            markerFile = markerGene.getDatabase()
            my_logs.append_logs("Marker gene database...:" + '\n' + markerFile)
            db = markerGene.loadDatabase(markerFile)
        else: 
            my_logs.append_logs("Marker gene database...:" + '\n' + markerFile)
            markerGene = loadMarkerGeneDB(outroot, db_species, tissues)
            db = markerGene.loadDatabase(markerFile)
        
        ## filter the marker gene database for correct cell and tissue types 
        db = db[db['cancer_type'] == "Normal"]
        db = db[db['tissue_type'].isin(tissues)]
        db = db.reset_index()
        db = db.drop_duplicates()
        ##################################################################################################
        ## convert gene symbols/IDs to the species of sample 
        if find_ortho:
            ## convert human gene symbols to macaque gene symbols
            ortho_species = species
            ortho_species_biomart = (ortho_species[0] + re.sub('.*_','',ortho_species)).lower()
            orthoClass = findOrthos(database=db, db_species=db_species, ortho_species=ortho_species_biomart)
            ## get the ensembl IDs from gene symbols for reference species
            db = orthoClass.symbol_to_ensembl()
            db = db[['cell_name','Symbol','ensembl_id']].dropna().drop_duplicates()
            ## get the orthologus gene symbols and ensembl IDs
            orthoDB = orthoClass.queryBioMart(db) 
            orthoDB = orthoDB.dropna().drop_duplicates()
            ## merge the two data
            db = db.merge(orthoDB, left_on='ensembl_id', right_on='ref_ensembl', how='inner')
            ## convert marker gene dataframe to a dictionary 
            db = db[db['ortho_symbol'].notna()]
            geneDict = db.groupby('cell_name')['ortho_symbol'].apply(list).to_dict()
        else:
            ## convert marker gene dataframe to a dictionary
            db = db[db['Symbol'].notna()]
            geneDict = db.groupby('cell_name')['Symbol'].apply(list).to_dict()
    else:
        markerGene = loadMarkerGeneDB(outroot, db_species, tissues)
        db = markerGene.PangloDB()
        if find_ortho:
            ortho_species = species
            ortho_species_biomart = (ortho_species[0] + re.sub('.* ','',ortho_species)).lower()
            orthoClass = findOrthos(database=db, db_species=db_species, ortho_species=ortho_species_biomart)
            ## get the ensembl IDs from gene symbols for reference species
            db = orthoClass.symbol_to_ensembl_panglo()
            db = db[['cell type','Symbol','ensembl_id']].dropna().drop_duplicates().rename(columns = {'cell type':'cell_name'})
            ## get the orthologus gene symbols and ensembl IDs
            orthoDB = orthoClass.queryBioMart(db) 
            orthoDB = orthoDB.dropna().drop_duplicates()
            ## merge the two data
            db = db.merge(orthoDB, left_on='ensembl_id', right_on='ref_ensembl', how='inner')
            ## convert marker gene dataframe to a dictionary 
            db = db[db['ortho_symbol'].notna()]
            geneDict = db.groupby('cell_name')['ortho_symbol'].apply(list).to_dict()
        else:
            ## convert marker gene dataframe to a dictionary
            db = db[db['Symbol'].notna()]
            geneDict = db.groupby('cell_name')['Symbol'].apply(list).to_dict()
    
    ## save the the marker gene database to a csv file
    outfile = os.path.join(pdfDir, 'GSEA_markers.csv')
    db.to_csv(outfile, sep=',', header=True, index=False)

    ##################################################################################################
    ## initialize the differential expression / GSEA class
    DE = DifferentialExpression(adata, geneDict = geneDict, numDEGs=numDEGs)
    groupVars = [x for x in adata.obs.columns.to_list() if 'leiden' in x] ## get the cell meta data associated with cell clustering (leiden algorithm)

    ## run differential expression and store in adata object 
    my_logs.append_logs("\nRunning differential expression...:")
    my_logs.append_logs(["\t>grouping by " + "\'" + x + "\'" for x in groupVars])
    
    ##################################################################################################
    ## run differential expression and gene set enrichment analysis for all clusterings 
    termDat = pd.DataFrame(); # pred = dict()
    for groupVar in groupVars:
        try:
            DE.runDE(groupVar = groupVar, layer='sqrt_norm')
        except ZeroDivisionError:
            groupVars = [x for x in groupVars if x != groupVar]
            continue
        clusters = np.unique(adata.obs[groupVar]); 
        for clust in clusters:
            dat = DE.GSEA(group = clust)
            if dat is not None:
                dat['resolution'] = groupVar
                termDat = pd.concat([termDat, dat])
    
    ## convert the p-value for plotting in outputs
    termDat['log10-adj-pvalue'] = -1*np.log10(termDat['Adjusted P-value'])

    ##################################################################################################
    ## iterate through the cell type annotations and pick the most appropriate cell type based on GSEA scores 
    k=0; termTmp = termDat[termDat['resolution'] == sorted(termDat['resolution'].unique())[0]]
    ## Get the highest scoring cell type for cluster set with the lowest resolution (fewest number of clusters)
    scoreDat = termTmp.iloc[termTmp.reset_index().groupby(['group'])['log10-adj-pvalue'].idxmax()][['resolution','group','Term', 'log10-adj-pvalue']]
    while k < scoreDat.shape[0]: ## iterate through every cluster annotation
        res_one = scoreDat['resolution'].iloc[k] ## current clustering resolution
        cluster = scoreDat['group'].iloc[k] ## current cluster from that resolution
        CT = scoreDat['Term'].iloc[k] ## current cell type for the cluster 
        res_comp_list = groupVars[groupVars.index(res_one) + 1:] ## comparison resolutions
        for res_two in res_comp_list:
            ## Compare the score of annotation cluster ${cluster} from resolution $(res_one) to all the clusters from $(res_two) that map to the same cells
            t1 = pd.DataFrame(adata.obs[adata.obs[res_one] == cluster][res_two].value_counts())
            t1.columns = [res_two]
            t1['pct'] = t1[res_two] / sum(t1[res_two])
            ## if less than 10% of cells from a cluster at $res_two are associated with this current cluster disregard for comparison
            comp_clusters = t1[t1['pct'] > 0.1].index.to_list()
            for cluster_two in comp_clusters:
                ## get the term Dat for this cluster at this resolution 
                termDat_tmp =  termDat.loc[(termDat['resolution'] == res_two) & (termDat['group'] == (cluster_two))]
                ## sort data frame based on p-value
                termDat_tmp = termDat_tmp.sort_values("Adjusted P-value").reset_index().drop('index', axis = 1)
                ## check if the top term (cell type) is the same as the previous cell type
                try: 
                    newCT = termDat_tmp['Term'].iloc[0]
                    # if the new cell type is different than the parent cell type 
                    if newCT != CT:
                        ## compare the p-value of newCT and CT for this cluster at this resolution
                        try:
                            newCT_pval = termDat_tmp[termDat_tmp['Term'] == newCT].reset_index()['log10-adj-pvalue'].iloc[0]
                            CT_pval = termDat_tmp[termDat_tmp['Term'] == newCT].reset_index()['log10-adj-pvalue'].iloc[0]
                        except:
                            pdb.set_trace()
                            
                        ## if the new cell type term is a more than 1.5X more significant update the cell type
                        if newCT_pval > 1.5 * CT_pval:
                            newDat = termDat_tmp.iloc[0:1][['resolution','group','Term', 'log10-adj-pvalue']]
                            scoreDat = pd.concat([scoreDat, newDat])
                            break
                except IndexError:
                    continue
        k = k + 1

    ## saving the gsea results to text file 
    outfile = os.path.join(pdfDir, 'gseaRes.csv')
    my_logs.append_logs('Saving GSEA results as...: ' + outfile)
    termDat = termDat[['resolution', 'group','Gene_set','Term','Overlap','Genes','P-value','Adjusted P-value','log10-adj-pvalue','Odds Ratio', 'score']]
    termDat.to_csv(outfile, sep=',', header=True, index=False)
    outfile = os.path.join(pdfDir, 'gseaMaxRes.csv')
    scoreDat.to_csv(outfile, sep=',', header=True, index=False)
    
    ##################################################################################################
    ## filter the gsea resuts for the best (lowest p-value)
    my_logs.append_logs('\nCell type by cluster:')
    my_logs.append_logs(scoreDat)
    adata.obs['celltype'] = "NA"
    adata.obs['log10-adj-pvalue'] = 0
    for i in np.arange(scoreDat.shape[0]):
        var = scoreDat['resolution'].iloc[i]
        clust = scoreDat['group'].iloc[i]
        adata.obs['celltype'] = np.where(adata.obs[var] == clust, scoreDat['Term'].iloc[i], adata.obs['celltype'])
        adata.obs['log10-adj-pvalue'] = np.where(adata.obs[var] == clust, scoreDat['log10-adj-pvalue'].iloc[i], adata.obs['log10-adj-pvalue'])

    ##################################################################################################
    ## Z-score annotation: utilizes the score_genes function from scanpy and annotates each cell individually
    ## based on relative expression of marker genes. Has been found to be useful for retinal datasets.
    adata.layers["pearson"] = adata.X.copy()
    adata.X = adata.layers['sqrt_norm'] ## move the raw counts back into matrix slot
    celltypes = list(geneDict.keys())
    for ct in celltypes:
        ct_markers = geneDict[ct]
        try:
            try:
                ## The score is the average expression of a set of genes subtracted with the average expression of a reference set of genes. 
                ## The reference set is randomly sampled from the gene_pool for each binned expression value.
                sc.tl.score_genes(adata, ct_markers, score_name='{}_scores'.format(re.sub(r'\s+','_',ct)), use_raw = False)
                score_mean = np.mean(adata.obs['{}_scores'.format(re.sub(r'\s+','_',ct))])
                score_sd = np.std(adata.obs['{}_scores'.format(re.sub(r'\s+','_',ct))])
                ## calculate a z-score from the score_genes function with (x - mean(x) / sd(x))
                adata.obs['{}_z_scores'.format(re.sub(r'\s+','_',ct))] = (adata.obs['{}_scores'.format(re.sub(r'\s+','_',ct))] - score_mean) / score_sd
            except KeyError:
                continue           
        except ValueError:
            continue
    
    ## get the highest scored cell type for each cell (barcode)
    barcode_list = adata.obs.index.to_list()
    dat = adata.obs[[x for x in adata.obs.columns if 'z_scores' in x]].reset_index().rename(columns = {'index':'barcode'}).melt('barcode')
    dat = dat.iloc[dat.groupby('barcode')['value'].agg(pd.Series.idxmax)]

    ## assign barcodes and merge with adata.obs
    dat['celltype_z_score'] = dat['variable'].apply(lambda x: re.sub("_z_scores","",x))
    adata.obs = pd.merge(adata.obs.reset_index(), dat[['barcode','celltype_z_score']], left_on = 'index', right_on = 'barcode').set_index('index')

    groupVars = [x for x in adata.obs.columns if "leiden" in x]
    for groupVar in groupVars:
        t1 = adata.obs.groupby([groupVar,'celltype_z_score']).size().reset_index()
        t1 = t1.rename(columns = {0:'size'})
        t2 = pd.merge(t1 ,t1.groupby(groupVar)['size'].sum().reset_index().rename(columns =  {'size':'total'}))
        t2['fraction'] = t2['size'] / t2['total']
        t2 = t2.iloc[t2.groupby([groupVar])['fraction'].agg(pd.Series.idxmax)][[groupVar,'celltype_z_score']].rename(columns = {'celltype_z_score':'{}_celltype'.format(groupVar)})
        t2[groupVar] = t2[groupVar].apply(str)
        adata.obs = pd.merge(adata.obs.reset_index(), t2, on = groupVar).set_index('index')
        adata.obs = adata.obs.loc[barcode_list]

    ##################################################################################################
    ## save the updated adata object
    adata.X = adata.layers['raw'] ## move the raw counts back into matrix slot
    if args.overwrite:
        my_logs.append_logs('\nOverwriting the adata object...')
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
        my_logs.append_logs('\nSaving adata object as : ' + outfile)
        my_logs.append_logs(breaker)
        adata.write(outfile)
        
if __name__ == "__main__":
    main()