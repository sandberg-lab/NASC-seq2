To calculate degradation rates for the demo data and a fictional labeling time of 1.23 hours, run (tales ~1 second): 
`python3 ../calculate_degradation_rate.py --bootstrap 50 demo_newrna_reconstructedmolecules.csv demo_oldrna_reconstructedmolecules.csv 1.23 demo_degradation_output.csv`
In addition to the output table, it will also tell the mean degradation rate in stdout.

Either use the existing .so file (already compiled for linux ubuntu) or compile it using:
`sudo apt install libflint-arb-dev`  # or one of the other installation options at https://`arblib.org/`setup.html`
`make`

The main kinetic-parameters-to-summary-statistics table we used for the study can be found under the name table3prox_2h_7deg0.065.trimmed.tsv. To use a different degradation rate and/or time, for example a degradation rate of 0.258/h and a time of 1.23 hours, run (takes 2-3 minutes):
`python3 ../parameter_table_calc_from_prob_3proxies.py --open_rate 0.002 50 8 --transcribe_rate 1 200 4 --degrade_rate 0.258 0.258 1 --close_rate 100 100 1 --time 1.23 1.23 1 --proc 260 --prec 10000 --table_out demo_m_table.tsv`
Then remove extreme values that get in the way of linear interpolation (takes ~1 second):
`python3 ../trim_lookup_table.py demo_m_table.tsv`
This will create `demo_m_table.trimmed.tsv`, this file goes after -m in bootstrap_nonzero_three_estimate_lookup_one_csv_v5.py Replace 8 with 73 and 4 with 38 for more accurate downstream results (takes 20-30 minutes), and adjust --proc based on your CPUs. 

Now you can run the initial estimation of the burst kinetic parameters (takes ~20 seconds):
`python3 ../bootstrap_nonzero_three_estimate_lookup_one_csv_v5.py demo_newrna_reconstructedmolecules.csv demo_kinetics_from_lookup.tsv -m demo_m_table.trimmed.tsv --lookup_mode interp nearest noborders --proc 4 1 --reuse_lookups --bootstraps 50`

And then refine the numbers (especially lowly expressed genes) by running (takes ~25 minutes):
`python3 ../transient_ml_v2.py -i demo_newrna_reconstructedmolecules.csv --guesses demo_kinetics_from_lookup.tsv -d 0.258 --prec 10000 --time 1.23 --threads 200 -o demo_kinetics_from_maximum_likelihood.tsv`

Done. Now an example of how to reuse one of the plotting scripts. To look at the parameter distribution, after filtering on the bootstrap distributions, you can run (takes ~3 seconds):
`python3 -Wignore ../kinetics_plotting/plot_parameter_distributions_with_bootstrap_filters_v3.py demo_kinetics_from_maximum_likelihood.csv --degradation demo_degradation_output.csv --expression_level_basis demo_kinetics_from_lookup.tsv -r 1 2 50 --conf 50 --single_plot --draw_CI -o demo_fictional_parameter_distribution_plot.pdf --sep ,`
Due to the filtering and the few genes in the demo file, not all parameters can be plotted. Ignore the 'posx and posy should be finite values' warning and the gmean values.