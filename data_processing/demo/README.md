To run this pipeline from the demo fastq files on Linux, first install stitcher and zUMIs, for example by running (takes ~2 minutes):  
`git clone https://github.com/AntonJMLarsson/stitcher.py`  
`git clone https://github.com/sdparekh/zUMIs.git`  
Also install STAR 2.7.3a from https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz

Download the mouse genome and index the genome (takes about ~20 minutes):  
`wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`  
`gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`  
`mkdir mm39_starindex`  
`STAR --runThreadN 50 --runMode genomeGenerate --genomeDir mm39_starindex --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa`  

Decompress the relevant files (takes ~10 seconds):  
`bunzip2 ../../data_tables/NASCseqPF3.final_annot.gtf.bz2 NASCseqPF3.final_annot.intervals.json.bz2 NASCseqPF3.final_annot.refskip.json.bz2`  

And now run zUMIs:  
`zUMIs/zUMIs.sh -c -y demo_zUMIs_NASCseq2.yaml`  
