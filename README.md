## "Single-cell new RNA sequencing reveals principles of transcription at the resolution of individual bursts".
This repository contains scripts and tabular data related to 4sU-based single-cell RNA-seq for improved inference of transcriptional bursting kinetics, and transcriptome-wide analysis of co-bursting.

We provide the following tabular data: New RNA counts tables (at allelic resolition) and transcriptome-wide inferred kinetic parameters. Both these are found in the folder `data_tables` of this repository.

For script, there are three catogories, corresponding to the sections of the manuscript:

In the folder `data_processing`, you can find scripts and instructions on how to go from NASC-Seq2 sequencing output to molecular counts.

In the folder `kinetics_lookup`, you can find script for generating burst kinetic estimates from molecular counts, and code we used for plots in the manuscript, simulations and format conversions.

In the folder `cobursting`, you can find script that were used for our cobursting analysis.

Some of the folders contain additional readme files for further information.
