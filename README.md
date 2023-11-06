# NASC-seq2
This repository contains scripts and tabular data related to a forthcoming manuscript that describes the RNA library preparation method NASC-Seq2 and presents data analysis using such data to study transcriptional burst kinetics and gene cobursting.

The tabular data, meaning count tables and estimated burst kinetics from the study, can be found in folder `data_tables` of this repository.

For script, there are three catogories, corresponding to the sections of the manuscript:

In the folder `data_processing`, you can find scripts and instructions on how to go from NASC-Seq2 sequencing output to molecular counts.

In the folder `kinetics_lookup`, you can find script for generating initial burst kinetic estimates from molecular counts. These estimates can be refined using the scripts in the folder `inference_ML`. To see the script we used to plot the resulting data, look in the folder `kinetics_plotting`, these script will also use some of the output of the scripts in `misc_kinetics`.

In the folder `cobursting`, you can find script that were used for our cobursting analysis.
