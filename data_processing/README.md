# NASC-seq2 scripts

## NASCseqV2.py

This function allows you to easily run some (or all) of the steps of the NASC-seq2 analysis pipeline.
To do this, the following information should be added to the zUMIs yaml file:

```yaml

NASC-seq:
  NASC-stage: 'all' ## Possible values ('stitcher','tag','extract','estim_pc','hyptest','all') ## Use 'all' to run all steps directly.
  script_dir: '~/programs/NASC-seqV3' ## Directory where to find scripts for NASC-seq data analysis
  verbose: False ## Should commands actually be ran. Used for testing...
  vcf: 'vcfFilePath.vcf' ## Optional, add in vcf file to be used for allele assignment
  fasta:  '~/resources/genomes/human/hg38.primary_assembly.sorted.fa' ## The fasta file corresponding to the genome used for alignment of the data.
  isoform: '~/programs/stitcher.py/hg38_unique_intervals_for_isoforms.json' ## .json file with isoform information.
  junction: '~/programs/stitcher.py/hg38_unique_refskip_for_isoforms.json' ## .json file with junction information.
  stitcher_exec: '~/programs/stitcher.py/stitcher.py' ## Location of stitcher.py
  spikeID: 'diySpike' ## Pattern identifying spike reads that are to be treated as such.
  
```
In addition to the produced datafiles, logfiles of both executed commands and warning messages/error messages are automatically saved to the experiment directory.
Please note that some of the zUMIs settings (such as num_threads and num_memory) are used in NASCseqV3.py. If you want to use less resources, you will have to change these settings there. 

### Dependencies

#### Python:
* joblib
* pandas
* numpy
* yaml
* h5py
* scipy
* multiprocessing
* pysam
* pygtrie
* itertools
* vcf
* gc
* seaborn
* pyfaidx
* portion
* collections

#### R:
* data.table
* dtplyr
* dplyr
* parallel
* Rsamtools
* tictoc

### Usage

NASCseqV2.py [-y yamlfile]

**arguments:**
```
-h, --help show this help message and exit
-y yamlfile, yamlfile with both zUMIs and NASC-seq options
```


## tag_molecules.py

### Usage
tag_molecules.py [-h] [-i input] [-o output] [-g GTF] [-f FASTA]
                        [-v VCF] [-mut MUTATIONS] [-t threads]
                        [--contig contig]

**arguments:**
```
  -h, --help show this help message and exit
  -i input, --input input Input .bam file
  -o output, --output output Output .bam file
  -g GTF, --gtf GTF     gtf file with gene information
  -f FASTA, --fasta FASTA Fasta file for your reference genome
  -v VCF, --vcf VCF     vcf file with genotype information
  -mut MUTATIONS, --mutations MUTATIONS Output hdf5 file with mutation information
  -t threads, --threads threads Number of threads
  --contig contig       Restrict stitching to contig
```
## extract_tags.py

### Usage
extract_tags.py [-h] [-i input] [-o output] [-g gtf] [-t threads]
                       [--contig contig]
                      
**arguments:**
```
  -h, --help            show this help message and exit
  -i input, --input input Input .bam file
  -o output, --output output Output .h5 file
  -g gtf, --gtf gtf     gtf file with gene information
  -t threads, --threads threads Number of threads
  --contig contig       Restrict extraction to contig
```
## estim_pc.py

### Usage
estim_pc.py [-h] [-h5 input] [-t threads]

**arguments:**
```
  -h, --help            show this help message and exit
  -h5 input, --hdf5 input .h5 file to process
  -t threads, --threads threads Number of threads
```
To run the separate parts of the pipeline, you can either modify the NASC-stage flag in the yaml file, or call on the individual scripts directly by using the information below.
