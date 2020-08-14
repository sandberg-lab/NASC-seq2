#!/usr/bin/env python

import os, argparse, sys, subprocess, time, glob, ntpath, shutil
from joblib import Parallel, delayed
import pandas as pd
import yaml

## # TODO 14.8.2020:

### Yaml parameter naming
### File handling when files already exist?
### Cleanup of temporary files
### Reporting structure to the user and logging
### List of required packages and versions
### Flag / partial run handling
### Optional VCF handling

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-y','--yamlfile',required=True)
	o = parser.parse_args()

## yaml file handling and extraction
with open(o.yamlfile, 'r') as stream:
    try:
        yamldata=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

## Assign yaml options to variables (using zUMIs yaml file for ease of use)
experimentdir = yamldata['out_dir']
scriptpath = yamldata['NASC-seq']['script_dir']
numCPU = yamldata['num_threads']
mem_limit = yamldata['mem_limit']
verbose = yamldata['NASC-seq']['verbose']
gtf = yamldata['reference']['GTF_file']
fasta = yamldata['NASC-seq']['fasta']
vcf = yamldata['NASC-seq']['vcf']
isoform = yamldata['NASC-seq']['isoform']
stitcher_exec = yamldata['NASC-seq']['stitcher_exec']
NASCflag = yamldata['NASC-seq']['NASC-stage']
python_exec = 'python3'
## File and folder handling functions
def safe_mkdir(f):
	if not os.path.exists(f):
		os.mkdir(f)

def run_cmd(cmd,commandlogfile,verbose):
    commandlogfile.write('%s\n' % " ".join(cmd))
    if verbose=='FALSE':
        subprocess.call(" ".join(cmd),shell=True)

def sort_bam(infile,outfile,numCPU,mem_limit,commandlogfile,verbose):
    sortThreads = int(int(numCPU)/10) ## No need for 100-200 CPUs for this process...
    sortMem = ''.join([int(mem_limit)/sortThreads),'G']
    run_cmd(['samtools sort',infile,'-o',outfile,'-m',sortMem,'-@',sortThreads],commandlogfile,verbose=verbose)
    run_cmd(['samtools index',outfile,'-@',sortThreads],commandlogfile,verbose=verbose)

## Prepare directories and files
safe_mkdir(os.path.join(experimentdir,'NASC-seq'))
safe_mkdir(os.path.join(experimentdir,'NASC-seq','logfiles'))
commandlogfile=open(os.path.join(experimentdir,'NASC-seq','logfiles','commandlog.txt'), 'a')

## Processing
if NASCflag=='stitcher' or NASCflag=='all':
    infile = os.path.join(experimentdir,yamldata['project']+'.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam')
    run_cmd(['samtools index',infile,'-@',int(numCPU)],commandlogfile,verbose=verbose)
    outfile = os.path.join(experimentdir,yamldata['project']+'.stitched.bam')
    logfile = os.path.join(experimentdir,'NASC-seq','logfiles','stitcherlog.txt')
    run_cmd(['nohup',python_exec,stitcher_exec,'--input',infile,'--o',outfile,'--g',gtf,'--iso',isoform,'--t',int(numCPU),'>',logfile,'2>&1'],commandlogfile,verbose=verbose)
    outfileSorted = os.path.join(experimentdir,yamldata['project']+'.stitched.sorted.bam')
    sort_bam(outfile,outfileSorted,numCPU,mem_limit,commandlogfile,verbose)

if NASCflag=='tag' or NASCflag=='all':
    infile = os.path.join(experimentdir,yamldata['project']+'.stitched.sorted.bam')
    outfile = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.bam')
    mutfile = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.mutationfile.hdf5')
    logfile = os.path.join(experimentdir,'NASC-seq','logfiles','taglog.txt')
    run_cmd(['nohup',python_exec,os.path.join(scriptpath,'tag_molecules.py'),'-i',infile,'-o',outfile,'-g',gtf,'-f',fasta,'-mut',mutfile,'-t',int(numCPU),'>', logfile, '2>&1'],commandlogfile,verbose=verbose)
    outfileSorted = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.sorted.bam')
    sort_bam(outfile,outfileSorted,numCPU,mem_limit,commandlogfile,verbose)

if NASCflag=='extract' or NASCflag=='all':
    infile = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.sorted.bam')
    outfile = os.path.join(experimentdir,yamldata['project']+'.moleculeInformation.h5')
    logfile = os.path.join(experimentdir,'NASC-seq','logfiles','extractlog.txt')
    run_cmd(['nohup',python_exec,os.path.join(scriptpath,'extract_tags.py'),'-i',infile,'-o',outfile,'-g',gtf,'-t',int(numCPU),'>',logfile , '2>&1'],commandlogfile,verbose=verbose)

if NASCflag=='estim_pc' or NASCflag=='all':
    infile = os.path.join(experimentdir,yamldata['project']+'.moleculeInformation.h5')
    logfile = os.path.join(experimentdir,'NASC-seq','logfiles','estimlog.txt')
    run_cmd(['nohup',python_exec,os.path.join(scriptpath,'estim_pc.py'),'-h5',infile,'-t',int(numCPU),'>',logfile,'2>&1'],commandlogfile,verbose=verbose)

if NASCflag=='hyptest' or NASCflag=='all':
    infile = os.path.join(experimentdir,yamldata['project']+'.moleculeInformation.h5')
    logfile = os.path.join(experimentdir,'NASC-seq','logfiles','hyptestlog.txt')
    run_cmd(['nohup',python_exec,os.path.join(scriptpath,'do_htest.py'),'-h5',infile,'-t',int(numCPU),'>',logfile,'2>&1'],commandlogfile,verbose=verbose)
