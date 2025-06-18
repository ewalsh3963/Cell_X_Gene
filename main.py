# !/bin/python3
# QBDS
# Evan Walsh (evan.walsh@capsida.com)

import argparse
import logging
import os
from datetime import datetime
import re
import sys
import time
import glob
import lowD
import GSEA
import Report
import download
import pandas as pd
import json
import pdb
import sys
import OS_Tools

class CommandLine:
    """ Handle the command line, usage and help requests """

    def __init__(self, in_opts=None):
        """ CommandLine constructor: Implements a parser to interpret the command line argv string using argparse. """
        self.parser = argparse.ArgumentParser(description="%(prog)s is the pipeline for downloading sc-RNAseq data with Cell X Gene. \n"
        "Modules in this pipeline include download, \n",
        epilog='CAPSIDA BIOTHERAPEUTICS sc-RNAseq analysis pipeline', 
        add_help=True, 
        prefix_chars='-', 
        usage='python3 %(prog)s [-h] [-v] [all, download, lowD, GSEA, Report] ...')
     
        ## version 
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0 - Capsida Biotherapeutics - Evan Walsh (evan.walsh@capsida.com)')

        subparsers = self.parser.add_subparsers(title="programs", 
                                                dest="program",
                                                metavar="[all, download, lowD, GSEA, Report]")
        ## Full Pipeline
        self.all_parse = subparsers.add_parser("all",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'all\' runs all modules to process sc-RNAseq data from Cell X Gene:\n"
                                                )
        self.all_args()

        ## Download
        self.download_parse = subparsers.add_parser("download",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'download\' download sc-RNAseq anndata from Cell X Gene and run QAQC")
        self.download_args()

        ## lowD
        self.lowD_parse = subparsers.add_parser("lowD",
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                description="\'lowD\' is designed to embed cells in a intepretable low-dimensional manifold.\n"
                                                    "This transformation intends to capture high-dimensional variation in gene expression and project it into a interpretable\n"
                                                    "low-dimensional embedding. This module includes generation of plots to ensure that biological (rather than technical)\n"
                                                    "variation is captured in embeddings.")
        self.lowD_args()

        ## GSEA 
        self.GSEA_parse = subparsers.add_parser("GSEA",
                                                    formatter_class=argparse.RawTextHelpFormatter,
                                                    description="\'Gene Set Enrichment Analysis (GSEA)\' is meant to use DE analysis results to determine cell types.\n"
                                                    "Cell type is determined based on the overlap between DE genes and marker database genes")
        self.GSEA_args()

        ## Report 
        self.Report_parse = subparsers.add_parser("Report",
                                                    formatter_class=argparse.RawTextHelpFormatter,
                                                    description="\'Report\' generates a Quarto markdown report of the data with tables and visualizations")
        self.Report_args()

        if in_opts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(in_opts)
                
    def all_args(self):
        self.all_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.all_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='~/scratch')
        self.all_parse.add_argument('-overwrite', action='store_false', required = False, help='Overwrite the adata object after each module or save new [default = False]')

        ## download parameters
        self.all_parse.add_argument("-obs_value_filter", type=json.loads, default={}, help="python expression with selection conditions to fetch cells meeting a criteria")
        
        ## filtering thresholds 
        self.all_parse.add_argument('-min_genes', metavar='--min_genes', action='store', type=str, help='minimum number of genes detected in a cell for downstream analysis [default: quantile(nGenes, 0.05)]', required=False, default="0.05")
        self.all_parse.add_argument('-max_genes', metavar='--max_genes', action='store', type=str, help='maximum number of genes detected in a cell for downstream analysis [default: quantile(nGenes, 0.95)]', required=False, default="0.95")
        self.all_parse.add_argument('-min_transcripts', metavar='--min_transcripts', action='store', type=str, help='minimum number of transcripts detected in a cell for downstream analysis [default: quantile(nTscps, 0.05)]', required=False, default="0.05")
        self.all_parse.add_argument('-max_transcripts', metavar='--max_transcripts', action='store', type=str, help='maximum number of transcripts detected in a cell for downstream analysis [default: quantile(nGenes, 0.95)]', required=False, default="0.95")

        ## doublet detection 
        self.all_parse.add_argument('-solo', action='store_true', help='whether on not to run SOLO for doublet detection and removal [default: True]', required=False)
        self.all_parse.add_argument('-scrublet', action='store_false', help='whether on not to run scrublet for doublet detection and removal [default: True]', required=False)
        self.all_parse.add_argument('-scrt', metavar='--scrublet_threshold', action='store', help='doublet score threshold', required=False, default= None)
        self.all_parse.add_argument('-no_scrublet', action='store_false', help='whether on not to run scrublet for doublet detection and removal [default: True]', required=False)
        self.all_parse.add_argument('-exp_doub_rate', metavar='--expected_doublet_rate', action='store', type=int, help='estimated doublet rate [default: 0.05]', required=False, default=0.05)
        self.all_parse.add_argument('-sim_ratio', metavar='--sim_doublet_ratio', action='store', type=int, help='Number of doublets to simulate relative to the number of observed transcriptomes [default = 2]', required=False, default=2)
        self.all_parse.add_argument('-knn_metric', metavar='--knn_dist_metric', action='store', type=str, help='distance metric used when finding nearest neighbors', required=False, default='euclidean')

        ## low dimensional transformation parameters 
        self.all_parse.add_argument('-nHVGs', metavar='--num_highvar_genes', action='store', type=int, help='number of HVGs to select for DimRed [default: 2000]', required=False, default=2000)
        self.all_parse.add_argument('-recipe', action='store_true', help='run the full scanpy DimRed recipe', required=False)
        self.all_parse.add_argument('-no_recipe', action='store_false', dest='recipe', help='run Scanpy DimRed recipe stepwise', required=False)
        self.all_parse.add_argument('-theta', metavar='--theta', action='store', type=int, help='value of theta parameter in NB regression [default: 100]', required=False, default=100)
        self.all_parse.add_argument('-nComps', metavar='--num_of_comps', action='store', type=int, help='number of principal components to compute [default: 50]', required=False, default=50)
        self.all_parse.add_argument('-nNeigh', metavar='--num_of_neighbors', action='store', type=int, help=' size of local neighborhood used for manifold approximation [default = 15]', required=False, default=15)
        self.all_parse.add_argument('-res', metavar='--resolution', action='store', type=str,  nargs="*", help='value of the resolution parameter for clustering  [default: [0.1, 0.25, 0.5, 0.75]]', required=False, default='0.1 0.25 0.5 0.75')
        self.all_parse.add_argument('-pc_lim', action='store_false', help='number of principal components to use in low D transform (if True will calculate limit with kneeldle.kneed if false 30)', required=False, default=50)

        ## GSEA parameters
                
        ## PangloDB
        self.all_parse.add_argument('-PangloDB', action='store_true', help='use the PangloDB cell marker database', required=False)

        ## gsea parameters
        self.all_parse.add_argument('-markerFile', metavar='--markerFile', action='store', type=str, help='Marker gene see [see **** for format]', required=False, default=None)
        self.all_parse.add_argument('-tissues', metavar='--tissues', action='store', type=str, help='Tissues to filter marker genes for', nargs="*", required=False, default=None)
        self.all_parse.add_argument('-db_species', metavar='--db_species', action='store', type=str, help='Species to use for marker gene database [default: Human]', required=False, choices=['All','Human', 'Mouse', 'SingleCell'], default='Human')
        self.all_parse.add_argument('-find_ortho', action='store_true', help='find orthologous marker gene names', required=False)
        self.all_parse.add_argument('-no_find_ortho', action='store_false', dest='recipe', help='find orthologous marker gene names', required=False)
        self.all_parse.add_argument('-numDEGs', metavar='--num_of_DEGs', action='store', type=int, help='Top # of DEGs to use in Gene Set Enrichment analysis [default: all]', required=False, default=None)


    def download_args(self):
        self.download_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.download_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='~/scratch')
        self.download_parse.add_argument('-id', metavar='--dataset_id', action='store', help='ID for dataset(s) to be downloaded', nargs = "*", required=True, type=str)
        self.download_parse.add_argument('-org', metavar='--organism', action='store', help='organism [speices] of datasets', nargs = "*", required=True, type=str)
        self.download_parse.add_argument('-ref_org', metavar='--reference_organism', action='store', help='if multiple datasets from different species gene names will changed to reference organism gene symbols', required=False, type=str)
        self.download_parse.add_argument('-overwrite', action='store_false', required = False, help='Overwrite the adata object after each module or save new [default = False]')

        self.download_parse.add_argument("-obs_value_filter", type=json.loads, default={}, help="python expression with selection conditions to fetch cells meeting a criteria")

        ## filtering thresholds 
        self.download_parse.add_argument('-min_genes', metavar='--min_genes', action='store', type=str, help='minimum number of genes detected in a cell for downstream analysis [default: quantile(nGenes, 0.05)]', required=False, default="0.05")
        self.download_parse.add_argument('-max_genes', metavar='--max_genes', action='store', type=str, help='maximum number of genes detected in a cell for downstream analysis [default: quantile(nGenes, 0.95)]', required=False, default="0.95")
        self.download_parse.add_argument('-min_transcripts', metavar='--min_transcripts', action='store', type=str, help='minimum number of transcripts detected in a cell for downstream analysis [default: quantile(nTscps, 0.05)]', required=False, default="0.05")
        self.download_parse.add_argument('-max_transcripts', metavar='--max_transcripts', action='store', type=str, help='maximum number of transcripts detected in a cell for downstream analysis [default: quantile(nGenes, 0.95)]', required=False, default="0.95")

        ## doublet detection 
        self.download_parse.add_argument('-solo', action='store_true', help='whether on not to run SOLO for doublet detection and removal [default: True]', required=False)
        self.download_parse.add_argument('-scrublet', action='store_false', help='whether on not to run scrublet for doublet detection and removal [default: True]', required=False)
        self.download_parse.add_argument('-scrt', metavar='--scrublet_threshold', action='store', help='doublet score threshold', required=False, default= None)
        self.download_parse.add_argument('-no_scrublet', action='store_false', help='whether on not to run scrublet for doublet detection and removal [default: True]', required=False)
        self.download_parse.add_argument('-exp_doub_rate', metavar='--expected_doublet_rate', action='store', type=int, help='estimated doublet rate [default: 0.05]', required=False, default=0.05)
        self.download_parse.add_argument('-sim_ratio', metavar='--sim_doublet_ratio', action='store', type=int, help='Number of doublets to simulate relative to the number of observed transcriptomes [default = 2]', required=False, default=2)
        self.download_parse.add_argument('-knn_metric', metavar='--knn_dist_metric', action='store', type=str, help='distance metric used when finding nearest neighbors', required=False, default='euclidean')
    
    def lowD_args(self):
        self.lowD_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.lowD_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='~/scratch')

        ## low dimensional transformation parameters 
        self.lowD_parse.add_argument('-nHVGs', metavar='--num_highvar_genes', action='store', type=int, help='number of HVGs to select for DimRed [default: 2000]', required=False, default=2000)
        self.lowD_parse.add_argument('-recipe', action='store_true', help='run the full scanpy DimRed recipe', required=False)
        self.lowD_parse.add_argument('-no_recipe', action='store_false', dest='recipe', help='run Scanpy DimRed recipe stepwise', required=False)
        self.lowD_parse.add_argument('-theta', metavar='--theta', action='store', type=int, help='value of theta parameter in NB regression [default: 100]', required=False, default=100)
        self.lowD_parse.add_argument('-nComps', metavar='--num_of_comps', action='store', type=int, help='number of principal components to compute [default: 50]', required=False, default=50)
        self.lowD_parse.add_argument('-nNeigh', metavar='--num_of_neighbors', action='store', type=int, help=' size of local neighborhood used for manifold approximation [default = 15]', required=False, default=15)
        self.lowD_parse.add_argument('-res', metavar='--resolution', action='store', type=str,  nargs="*", help='value of the resolution parameter for clustering  [default: [0.1, 0.25, 0.5, 0.75]]', required=False, default='0.1 0.25 0.5 0.75')
        self.lowD_parse.add_argument('-pc_lim', action='store_false', help='number of principal components to use in low D transform (if True will calculate limit with kneeldle.kneed if false 30)', required=False, default=50)

        ## Adata save params
        self.lowD_parse.add_argument('-overwrite', action='store_false', required = False, help='Overwrite the adata object after each module or save new [default = False]')

    def GSEA_args(self):
        self.GSEA_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.GSEA_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='~/scratch')
        self.GSEA_parse.add_argument('-species', metavar='--species', action='store', required = True, help='Species of sample in format: Genus species')

        ## PangloDB
        self.GSEA_parse.add_argument('-PangloDB', action='store_true', help='use the PangloDB cell marker database', required=False)

        ## gsea parameters
        self.GSEA_parse.add_argument('-markerFile', metavar='--markerFile', action='store', type=str, help='Marker gene see [see **** for format]', required=False, default=None)
        self.GSEA_parse.add_argument('-tissues', metavar='--tissues', action='store', type=str, help='Tissues to filter marker genes for', nargs="*", required=False, default=None)
        self.GSEA_parse.add_argument('-db_species', metavar='--db_species', action='store', type=str, help='Species to use for marker gene database [default: Human]', required=False, choices=['All','Human', 'Mouse', 'SingleCell'], default='Human')
        self.GSEA_parse.add_argument('-find_ortho', action='store_true', help='find orthologous marker gene names', required=False)
        self.GSEA_parse.add_argument('-no_find_ortho', action='store_false', dest='recipe', help='find orthologous marker gene names', required=False)
        self.GSEA_parse.add_argument('-numDEGs', metavar='--num_of_DEGs', action='store', type=int, help='Top # of DEGs to use in Gene Set Enrichment analysis [default: all]', required=False, default=None)

        ## Adata save params
        self.GSEA_parse.add_argument('-overwrite', action='store_false', required = False, help='Overwrite the adata object after each module or save new [default = False]')    

    def Report_args(self):
        self.Report_parse.add_argument('-s', metavar='--study', action='store', type=str, required=True, help='Study title')
        self.Report_parse.add_argument('-o', metavar='--output_dir', action='store', help='output file path', required=False, type=str, default='~/scratch')

class LogInitializer:
    """ Create Logs depending on the program chosen. """
    def __init__(self, program, out):
        self.program = program
        self.out = out

    def logger_check(self):
        """ Create and check logs. """
        run_logs = os.path.join(self.out, "{}_scanpy.log".format(self.program))  # create logs file name
        OS_Tools.check_file(run_logs, critical=False, rewrite=True)  # check logs file - create new log
        return run_logs

    def param_logger(self, short, run_logs, title, args):
        """ Add input parameters to logs. """
        breaker = '=' * 120
        my_logs = OS_Tools.Logs(run_logs)
        my_logs.append_logs("{} {} ran with the following parameters:\n{}\n{}\n{}\n\n"
                            .format("RunSpipe.py", self.program, breaker, "\n".join(f'{k} = {v}' for k, v in vars(args).items()),breaker))

        my_logs.append_logs("{}:\n{} | {}\n".format(short, title, datetime.now()))

def main(command=None):
    """ Run specified program(s). """
    if command is None:
        my_command = CommandLine()  # read options from the command line
    else:
        my_command = CommandLine(command)  # interpret the list passed from the caller of main

    ## check if outroot exists 
    my_command.args.o = os.path.abspath(my_command.args.o)
    OS_Tools.ensure_directory(my_command.args.o, critical=False) 
    outroot = my_command.args.o

    # Handle top-level options
    if my_command.args.program == "all":
        """ Run all steps of the sc-RNASeq analysis pipeline  """

        init_log = LogInitializer("all", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
        init_log.param_logger("Running all sc-RNAseq analysis steps [download]", run_logs, my_command.args.s, my_command.args)

        ###################################################################################################
        ## Run Cell X Gene download
        my_logs.append_logs(breaker)
        my_logs.append_logs("Downloading fastq files from GEO using fastq_dump")
        outdir = os.path.join(outroot, 'expdata', my_command.args.s)
        
        ######################
        ## Run download and QAQC programs
        my_command.args.root = outroot
        my_command.args.o = outdir
        download.main(my_command.args, run_logs)
        
        ######################
        ## Run lowD program 
        lowD.main(my_command.args, run_logs)

        ######################
        ## Run GSEA program
        GSEA.main(my_command.args, run_logs)

        ######################
        ## Run report program
        Report.main(my_command.args, run_logs)
       
    elif my_command.args.program == "download":
        """ Download sc-RNAseq dataset(s) from Cell X Gene and run QAQC """
        outdir = os.path.join(outroot, 'scanpy', my_command.args.s)
        init_log = LogInitializer("download", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Download fastq files from GEO using fastq_dump", run_logs, my_command.args.s, my_command.args)
        
        ## Clean the adata objects that may already be present 
        adataDir = os.path.join(outdir, 'AdataFiles'); OS_Tools.ensure_directory(adataDir, critical = False)
        OS_Tools.clear_directory(adataDir, extensions = 'h5ad')

        ######################
        ## Run download program 
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        download.main(my_command.args, run_logs)

    elif my_command.args.program == "lowD":
        """ Low dimenionsal transformation of single-cell RNAseq data"""
        outdir = os.path.join(outroot, 'scanpy', my_command.args.s)
        init_log = LogInitializer("lowD", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Low dimenionsal transformation of single-cell RNAseq data", run_logs, my_command.args.s, my_command.args)
        
        ## Run 
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        lowD.main(my_command.args, run_logs)

    elif my_command.args.program == "GSEA":
        """ Cell Type annotation via Gene Set Enrichment analysis."""
        outdir = os.path.join(outroot, 'scanpy', my_command.args.s)
        init_log = LogInitializer("GSEA", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Annotate cell types with gene set enrichment analysis", run_logs, my_command.args.s, my_command.args)
        
        ## Run
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        GSEA.main(my_command.args, run_logs)

    elif my_command.args.program == "Report":
        """ Generate a Quarto markdown report of the analysis """
        outdir = os.path.join(outroot, 'CXG_Markdown')
        init_log = LogInitializer("Report", outdir)
        run_logs = init_log.logger_check()
        ## add command and paramters to 
        my_logs = OS_Tools.Logs(run_logs)
        breaker = "=" * 120
        my_logs.append_logs(breaker)
        my_logs.append_logs(" ".join(sys.argv))
                
        init_log.param_logger("Generate quarto markdown analysis summary report", run_logs, my_command.args.s, my_command.args)
        
        ## Run
        my_command.args.root = my_command.args.o
        my_command.args.o = outdir
        Report.main(my_command.args, run_logs)




if __name__ == "__main__":
    main()


# !import code; code.interact(local=vars())

###########################################################################
## TEST 
# python3 /home/ewalsh/Cell_X_Gene/main.py download -s MULTI_TEST \
# -id 856c1b98-5727-49da-bf0f-151bdb8cb056 842dd04e-5907-4da0-84af-90279bbf0902 \
# -org "Homo sapiens" "Mus musculus" \
# -ref_org "Homo sapiens" 

# python3 /home/ewalsh/Cell_X_Gene/main.py download -s TEST \
# -id 856c1b98-5727-49da-bf0f-151bdb8cb056 \
# -org "Homo sapiens" 

# python3 /home/ewalsh/Cell_X_Gene/main.py lowD -s TEST
 
# python3 /home/ewalsh/Cell_X_Gene/main.py GSEA -s TEST \
# -tissues Retina \
# -species "Homo sapiens" \
# -db_species "Human"

# python3 /home/ewalsh/Cell_X_Gene/main.py Report -s TEST
