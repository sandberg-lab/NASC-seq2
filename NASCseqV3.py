#!/usr/bin/env python
## GJH

import os, argparse, sys, subprocess, time, glob, ntpath, shutil
from joblib import Parallel, delayed
import pandas as pd
import yaml

## 3.9.20 TODO :

### File handling when files already exist?
### Cleanup of temporary files
### Reporting structure to the user and logging
### Error handling confirmation
### Move h5 file to NASC-seq folder for cleaner root?
###		Same for bam files produced in NASC-seq process?

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
python_exec = 'python3' ## Should probably be changed to yaml option...
R_exec = yamldata['Rscript_exec']
## File and folder handling functions
def safe_mkdir(f):
	if not os.path.exists(f):
		os.mkdir(f)

def run_cmd(cmd,commandlogfile,verbose):
    commandlogfile.write('%s\n' % " ".join(cmd))
    if str(verbose)=='False':
        subprocess.call(" ".join(cmd),shell=True)

def sort_bam(infile,outfile,numCPU,mem_limit,commandlogfile,verbose):
    sortThreads = int(int(numCPU)/10) ## No need for 100-200 CPUs for this process...
    sortMem = ''.join([str(int(int(mem_limit)/sortThreads)),'G'])
    run_cmd(['samtools sort',infile,'-o',outfile,'-m',sortMem,'-@',str(sortThreads)],commandlogfile,verbose=verbose)
    run_cmd(['samtools index',outfile,'-@',str(sortThreads)],commandlogfile,verbose=verbose)

def check_logfile(logfile,patterns = ['Error','error']):
	with open(logfile) as file:
    contents = file.read()
	x = 0
	for(pattern in patterns):
    	if pattern in contents:
        	x+=1
	if x>0:
		print('Errors have been detected... please check your logfile here: %s' % logfile)

## Prepare directories and files
safe_mkdir(os.path.join(experimentdir,'NASC-seq'))
safe_mkdir(os.path.join(experimentdir,'NASC-seq','logfiles'))
commandlogfile=open(os.path.join(experimentdir,'NASC-seq','logfiles','commandlog.txt'), 'a')

## Processing
if NASCflag=='stitcher' or NASCflag=='all':
	infile = os.path.join(experimentdir,yamldata['project']+'.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam')
	run_cmd(['samtools index',infile,'-@',str(numCPU)],commandlogfile,verbose=verbose)
	outfile = os.path.join(experimentdir,yamldata['project']+'.stitched.bam')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','stitcherlog.txt')
	indexfile = yamldata['barcodes']['barcode_file']
	run_cmd(['nohup',python_exec,stitcher_exec,'--input',infile,'--o',outfile,'--g',gtf,'--iso',isoform,'--t',str(numCPU),'--cells',indexfile,'>',logfile,'2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	outfileSorted = os.path.join(experimentdir,yamldata['project']+'.stitched.sorted.bam')
	sort_bam(outfile,outfileSorted,str(numCPU),mem_limit,commandlogfile,verbose)
	print('Finished stitching reads')

if NASCflag=='tag' or NASCflag=='all':
	infile = os.path.join(experimentdir,yamldata['project']+'.stitched.sorted.bam')
	outfile = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.bam')
	mutfile = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.mutationfile.hdf5')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','taglog.txt')
	run_cmd(['nohup',python_exec,os.path.join(scriptpath,'tag_molecules.py'),'-i',infile,'-o',outfile,'-g',gtf,'-f',fasta,'-mut',mutfile,'-t',str(numCPU),'>', logfile, '2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	outfileSorted = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.sorted.bam')
	sort_bam(outfile,outfileSorted,numCPU,mem_limit,commandlogfile,verbose)
	print('Finished tagging conversions')

if NASCflag=='qc' or NASCflag=='all':
	infile = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.sorted.bam')
	outfileSpike = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.conversionRates.diySpike.rds')
	outfile = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.conversionRates.rds')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','convratelog.txt')
	logfileSpike = os.path.join(experimentdir,'NASC-seq','logfiles','convratespikelog.txt')
	run_cmd(['nohup',R_exec,os.path.join(scriptpath,'conversionQC.R'),infile,outfile,str(numCPU),yamldata['NASC-seq']['spikeID'],'TRUE','>', logfile, '2>&1'],commandlogfile,verbose=verbose)
	run_cmd(['nohup',R_exec,os.path.join(scriptpath,'conversionQC.R'),infile,outfileSpike,str(numCPU),yamldata['NASC-seq']['spikeID'],'FALSE','>', logfileSpike, '2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	check_logfile(logfileSpike)
	print('Finished conversion rate QC')

if NASCflag=='extract' or NASCflag=='all':
	infile = os.path.join(experimentdir,yamldata['project']+'.stitched.tagged.sorted.bam')
	outfile = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.moleculeInformation.h5')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','extractlog.txt')
	run_cmd(['nohup',python_exec,os.path.join(scriptpath,'extract_tags.py'),'-i',infile,'-o',outfile,'-g',gtf,'-t',str(numCPU),'>',logfile , '2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	print('Finished extracting information from bam file and creating h5 file')

if NASCflag=='estim_pc' or NASCflag=='all':
	infile = os.path.join(experimentdir,'NASC-seq',yamldata['project']+'.moleculeInformation.h5')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','estimlog.txt')
	run_cmd(['nohup',python_exec,os.path.join(scriptpath,'estim_pc.py'),'-h5',infile,'-t',str(numCPU),'>',logfile,'2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	print('Finished estimating pc and pe for each cell')

if NASCflag=='hyptest' or NASCflag=='all':
	infile = os.path.join(experimentdir,,'NASC-seq'yamldata['project']+'.moleculeInformation.h5')
	logfile = os.path.join(experimentdir,'NASC-seq','logfiles','hyptestlog.txt')
	run_cmd(['nohup',python_exec,os.path.join(scriptpath,'do_htest.py'),'-h5',infile,'-t',str(numCPU),'>',logfile,'2>&1'],commandlogfile,verbose=verbose)
	check_logfile(logfile)
	print('Finished hypothesis testing')
	print('All done!')
