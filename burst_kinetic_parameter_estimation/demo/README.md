This demo shows how to use the script in this repository to go from a molecular count table, for single-cell NASC-Seq2 data or similar nascent RNA labeling data, to burst kinetic parameters of the telegraph model, that is, k<sub>on</sub>, k<sub>off</sub> and k<sub>syn</sub> where k<sub>on</sub> is the rate of transiting from OFF to ON state, k<sub>off</sub> is the rate from ON to OFF and k<sub>syn</sub> is the rate of transcription within the ON state. The demo runs off csv files, check at the bottom for using anndata files.

To run the code, you will need python3 (available at https://www.python.org/downloads/) and install python packages (see https://pip.pypa.io/en/stable/installation/ if pip is missing):  
`python3 -m pip install pandas scipy numpy joblib`  
The code has been tested on both Linux Ubuntu 18.04.6 with python 3.6.9 and Mac OS X 11.7.8 with python 3.8.0, however the scripts that require Arb (parameter_table_calc_from_prob_3proxies.py, transient_ml_v2.py) have only been run on Linux.

To calculate degradation rates for the demo data and a fictional labeling time of 1.23 hours, run (takes ~1 second):   
`python3 ../calculate_degradation_rate.py --bootstrap 50 demo_newrna_reconstructedmolecules.csv demo_oldrna_reconstructedmolecules.csv 1.23 demo_degradation_output.csv`  
In addition to the output table, it will also tell the mean degradation rate in stdout.  

Either use the existing hyp1f1_arb_vector.so file (already compiled for linux ubuntu) or recompile it using:  
`sudo apt install libflint-arb-dev`  # or one of the other installation options at https://arblib.org/setup.html  
`make`  

The main kinetic-parameters-to-summary-statistics table we used for the study can be found under the name table3prox_2h_7deg0.065.trimmed.tsv. To use a different degradation rate and/or time, for example a degradation rate of 0.258/h and a time of 1.23 hours, run (takes 2-3 minutes):  
`python3 ../parameter_table_calc_from_prob_3proxies.py --open_rate 0.002 50 8 --transcribe_rate 1 200 4 --degrade_rate 0.258 0.258 1 --close_rate 0.25 500 5 --time 1.23 1.23 1 --proc 260 --prec 10000 --table_out demo_m_table.tsv`  
Then remove extreme values that get in the way of linear interpolation (takes ~1 second):  
`python3 ../trim_lookup_table.py demo_m_table.tsv`  
This will create `demo_m_table.trimmed.tsv`, this file goes after -m in bootstrap_nonzero_three_estimate_lookup_one_csv_v5.py. Replace 8 with 73 (3rd argument for --open_rate), 4 with 38 (3rd argument for --transcribe_rate) and 5 with 55 (3rd argument for --close_rate) for more accurate downstream results (takes 20-30 minutes), and adjust --proc based on your CPUs. 

Now you can run the initial estimation of the burst kinetic parameters (takes ~20 seconds):  
`python3 ../bootstrap_nonzero_three_estimate_lookup_one_csv_v5.py demo_newrna_reconstructedmolecules.csv demo_kinetics_from_lookup.tsv -m demo_m_table.trimmed.tsv --lookup_mode interp nearest noborders --proc 4 1 --reuse_lookups --bootstraps 50`

And then refine the numbers (especially lowly expressed genes) by running (takes ~25 minutes):  
`python3 ../transient_ml_v2.py -i demo_newrna_reconstructedmolecules.csv --guesses demo_kinetics_from_lookup.tsv -d 0.258 --prec 10000 --time 1.23 --threads 200 -o demo_kinetics_from_maximum_likelihood.tsv`

Done. Now an example of how to reuse one of the plotting scripts. First install some python packages:  
`python3 -m pip install matplotlib seaborn statsmodels`  
To look at the parameter distribution, after filtering on the bootstrap distributions, you can run (takes ~3 seconds):  
`python3 -Wignore ../kinetics_plotting/plot_parameter_distributions_with_bootstrap_filters_v3.py demo_kinetics_from_maximum_likelihood.csv --degradation demo_degradation_output.csv --expression_level_basis demo_kinetics_from_lookup.tsv -r 1 2 50 --conf 50 --single_plot --draw_CI -o demo_fictional_parameter_distribution_plot.pdf --sep ,`  
Ignore the gmean values in the output.

If you start with an anndata file from data_processing/write_anndata_object.py called say mydata.h5ad, convert them first with:  
`python3 -m pip install anndata pandas`  
`python3 ../misc_kinetics/anndata_to_csv.py mydata.h5ad new > mydata_newrna.csv`  
`python3 ../misc_kinetics/anndata_to_csv.py mydata.h5ad old > mydata_oldrna.csv`  
and use them in place of demo_newrna_reconstructedmolecules.csv and demo_oldrna_reconstructedmolecules.csv.
