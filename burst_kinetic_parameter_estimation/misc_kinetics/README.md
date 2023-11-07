txburstML_pkl_to_csv.py converts the output from https://github.com/sandberg-lab/txburst/txburstML.py into the same csv format as the output of transient_ml_v2.py. For example, to convert mystreadystatekineticsML.pkl, run (requires pandas):  
`python3 txburstML_pkl_to_csv.py mystreadystatekineticsML.pkl mystreadystatekineticsML.csv`

anndata_to_csv.py converts the output of data_processing/write_anndata_object.py into count csv format. For example, to extract the layer new from mydata.h5ad, install the python packages anndata and pandas and then run:  
`python3 ../misc_kinetics/anndata_to_csv.py mydata.h5ad new > mydata_newrna.csv`  

simulate_data.py creates a count csv file based on a given grid distribution of burst kinetic parameters. For example, to create a 10x10 count matrix for 10 different k<sub>on</sub> values from 0.1 to 0.9 and 5 different k<sub>syn</sub> values from 50 to 100, the command would be (takes ~30 minutes, requires numpy):  
`python3 simulate_data.py 10 10 --open_rate 0.1 0.9 10 --transcribe_rate 50 100 5 --table_out simulation_matrix_out.csv` 
Similarly you can set k<sub>off</sub> values with --close_rate etc.
