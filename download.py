import os
import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scanpy as sc
import pdb
import QAQC
import scipy
import kneed
import gseapy
import subprocess
import cellxgene_census
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from pyensembl import EnsemblRelease
from matplotlib.pyplot import rc_context
import glob
import OS_Tools



def main(args, run_logs):
    """ Downalod dataset(s) from Cell X Gene"""
    ## Load the adata files and integrate them
    outroot, study, dataset_ids, organisms, ref_organism, obs_value_filter = args.o, args.s, args.id, args.org, args.ref_org, args.obs_value_filter

    ## initialize Log tools and directories
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120
    my_logs.append_logs(breaker)
    my_logs.append_logs("Initializing Cell X Gene download script")

    #####################################
    ## Downalod sc-RNAseq data with cellxgene_census.open_soma()
    # metaData = pd.read_csv(metadata_file, sep=',', header = 0)
    census = cellxgene_census.open_soma(census_version="2024-07-01")
    adata_list = []
    for id, organism in zip(dataset_ids, organisms):

        obs_string = [f"dataset_id == '{id}'"]
        for key, value in obs_value_filter.items():
            obs_string.extend([f"{key}=='{value}'"])

        my_logs.append_logs("Downloading from cellxgene_census with following obs_value_filters:")
        obs_string = ' '.join(obs_string)
        my_logs.append_logs(obs_string)
        adata = cellxgene_census.get_anndata(census=census, organism=organism, obs_value_filter=obs_string)
        
        ####################################################################################
        ## Run QAQC
        adata = QAQC.main(adata, args, id, organism, my_logs)

        if organism != ref_organism and ref_organism is not None:
            t1 = pd.DataFrame(adata.var['feature_name'])
            outfile = os.path.join(outroot, "{}_vars.csv".format(study))
            t1.to_csv(outfile, header = None, index = False)
            
            ## Get human gene symbols and read it in using biomart
            subprocess.call(["/usr/bin/Rscript", 
                            "/home/ewalsh/Cell_X_Gene/biomart_convert.R", 
                            outroot,
                            study,
                            organism])

            outfile = os.path.join(outroot, "{}_homo_sapiens_vars.csv".format(study))
            new_vars = pd.read_csv(outfile, sep=',', header = 0)
            new_vars.columns = ['species_symbol','species_ensembl_id','hsapien_symbol','hsapiens_ensembl_id']
            new_vars = new_vars[~new_vars['species_symbol'].duplicated(keep=False)] ## Remove gene if it shows up > 1
            new_vars = new_vars[~new_vars['hsapien_symbol'].duplicated(keep=False)] ## Remove ref orthology if it shows up > 1
            new_vars = new_vars.dropna()
            
            # Filter the AnnData object for genes with a reference organism ortholog
            org_symbols = adata.var['feature_name'].to_list()
            new_symbols = new_vars['species_symbol'].to_list()
            symbols_to_keep = list(set(org_symbols) & set(new_symbols))
            adata = adata[:, adata.var['feature_name'].isin(symbols_to_keep)]

            ## Change gene names to human symbol
            adata.var = pd.merge(adata.var, new_vars, left_on = 'feature_name', right_on = 'species_symbol', how = 'left').set_index('hsapien_symbol').drop(["species_symbol", "species_ensembl_id", "hsapiens_ensembl_id"], axis = 1)
            adata.var = adata.var.drop('feature_name', axis = 1)
        else:
            adata.var = adata.var.set_index('feature_name') 
        adata_list.append(adata)

    adata = adata_list[0].concatenate(adata_list[1::])
    
    ## save the adata object to its output directory 
    adataDir = os.path.join(outroot, 'AdataFiles')
    OS_Tools.ensure_directory(adataDir, critical = False)
    os.chdir(adataDir)
        
    ## list adata files
    adataFiles = sorted(glob.glob('[0-9]*.h5ad'))
    numFiles = len(adataFiles)
    fileID = numFiles + 1 

    outfile = os.path.join(adataDir, str(fileID) + '_adata.h5ad')
    my_logs.append_logs('Saving adata object as : ' + outfile)
    my_logs.append_logs(breaker)
    adata.write(outfile)

if __name__ == "__main__":
    main()
