This folder contains scripts for our cobursting analysis.

To run the plotting scripts in this folder, first install the required python packages:  
`python3 -m pip install matplotlib seaborn pandas anndata scipy numpy tqdm joblib`  

Then unsplit the input tables (corresponds to the files created by read.py, except in csv format due to github file size limitations):  
`cat ../data_tables/fibroblast_cobursting/mousefibroblasts_newrna_altallele_umicounts.csv.bz2.partaa ../data_tables/fibroblast_cobursting/mousefibroblasts_newrna_altallele_umicounts.csv.bz2.partab > mousefibroblasts_newrna_altallele_umicounts.csv.bz2`  
`cat ../data_tables/fibroblast_cobursting/mousefibroblasts_newrna_refallele_umicounts.csv.bz2.partaa ../data_tables/fibroblast_cobursting/mousefibroblasts_newrna_refallele_umicounts.csv.bz2.partab > mousefibroblasts_newrna_refallele_umicounts.csv.bz2`  

Then the two scripts that calculate correlations and plot them (takes ~20 minutes each):  
`python3 make_figures.py`  
`python3 make_figures2.py`  
Ignore the warnings regarding unexpressed genes (ConstantInputWarning). The script have been tested and run on both Linux Ubuntu python3.6 and on MacOSX python 3.8.
