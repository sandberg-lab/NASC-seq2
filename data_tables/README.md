### Data tables folder.

This folder contain new RNA counts matrices and inferred transcriptional bursting kinetics for primary fibroblasts and K562 cells. Note that the primary fibroblast data were derived from new RNA profiling in close to 10,000 cells, whereas K562 data were from around 600 cells. 

### new RNA count matrices



### Inferred transcriptome-wide bursting kinetics


### Inferred co-bursting


### General information on decompressing

Files named [something].bz2.part[suffix] can be decompressed in mac/linux bash using "cat [something].bz2.part* |bunzip2 -c > [something]", for example "cat fibroblast_cobursting/mousefibroblasts_newrna_altallele_umicounts.csv.bz2.part* |bunzip2 -c > fibroblast_cobursting/mousefibroblasts_newrna_altallele_umicounts.csv. On Windows the command should be "copy [something].bz2.partaa+[something].bz2.partab [something].bz2" followed by using a decompression tool.
