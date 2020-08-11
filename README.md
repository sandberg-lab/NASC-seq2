# NASC-seq2 scripts

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
