This demo starts from fastq files from a NASC-Seq2 experiment, aligns the reads, cluster them by molecule and processes them into an anndata file with separated counts for new and old RNA. The demo is written for Linux.

First install stitcher and zUMIs, for example by running (takes ~2 minutes):  
`git clone https://github.com/AntonJMLarsson/stitcher.py`  
`git clone https://github.com/sdparekh/zUMIs.git`  
Also install STAR from https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz or some other part of https://github.com/alexdobin/STAR that matches the version in zUMIs well enough.  

Download the mouse genome and index the genome (takes about ~20 minutes):  
`wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`  
`gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`  
`mkdir mm39_starindex`  
`STAR --runThreadN 50 --runMode genomeGenerate --genomeDir mm39_starindex --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa`  

Decompress the relevant files (takes ~10 seconds):  
`bunzip2 ../../data_tables/NASCseqPF3.final_annot.gtf.bz2 NASCseqPF3.final_annot.intervals.json.bz2 NASCseqPF3.final_annot.refskip.json.bz2`  

And now run zUMIs:  
`zUMIs/zUMIs.sh -c -y demo_zUMIs_NASCseq2.yaml`  

Install the dependencies needed for NASCseqV2.py:  
`python3 -m pip install joblib pandas numpy h5py scipy pysam pygtrie seaborn pyfaidx portion pyvcf pyyaml`  
`Rscript -e 'install.packages(c("data.table", "dtplyr", "dplyr", "parallel", "Rsamtools", "tictoc"))'`

And then the NASC-Seq2-specific steps:  
`python3 NASCseqV2.py -y demo_zUMIs_NASCseq2.yaml`
