import threading, time
from queue import Queue
import os
import sys
import re
from matplotlib import pyplot as plt
import argparse
import OS_Tools
import glob
import yaml
import scanpy as sc
import pdb
from quarto_python import quarto_profile, quarto_direct, preview_profile

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

def main(args, run_logs):
    ## Get arguments from argparse
    outroot, study, quarto_dir = args.o, args.s, args.o

    ## initialize Log tools
    my_logs = OS_Tools.Logs(run_logs)
    breaker = "=" * 120

    ############################################
    ## 1) Create a sample report for every sample in the /scanpy/$study directory 
    scanpyDir = os.path.join("~/scratch/scanpy", study)
    current_working_directory = os.getcwd()
    adataDir = os.path.join(scanpyDir, 'AdataFiles')
    os.chdir(adataDir)
    
    ## list adata files
    adataFiles = sorted(glob.glob('[0-9]*.h5ad'))
    adataFile = os.path.join(adataDir, adataFiles[len(adataFiles)-1])
    if not os.path.exists(adataFile):
        sys.exit("ERROR: " + str(adataFile) + " object does not exist...")

    os.chdir(current_working_directory)
    
    ## yml file parameters:
    output_dir = f"outdir/{study}/"
    yaml_dict = {'project':{'output-dir':output_dir}, 
                    'adata_file':adataFile, 
                    'study': study,
                    "rootDir": scanpyDir}
    yaml_write_file = os.path.join(quarto_dir, f'_quarto-{study}.yml')
    with open(yaml_write_file, 'w') as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)
        
    ######################################################################
    ## Render a report for this sample 
    os.chdir(quarto_dir)
    quarto_profile(f"{study}")

    ######################################################################
    try: 
        index_files = OS_Tools.find_files(parent_dir = os.path.join(quarto_dir, 'outdir', study), extension = 'index.html', check=None, wd=os.getcwd(), recursive=True)
        if type(index_files) != list:
            index_files = [index_files]
        for index_file in index_files:
            OS_Tools.delete_file_if_exists(index_file)
    except FileNotFoundError:
        pass
    try:
        about_files = OS_Tools.find_files(parent_dir = os.path.join(quarto_dir, 'outdir', study), extension = 'about.html', check=None, wd=os.getcwd(), recursive=True)
        if type(about_files) != list:
            about_files = [about_files]
        for about_file in about_files:
            OS_Tools.delete_file_if_exists(about_file)
    except FileNotFoundError:
        pass

if __name__ == "__main__":
    main()


