BBN_Graphs.py has all the functions necessary to generate diagnostics and element plots. If the .py is opened, each function has notes on how to use it and what exactly it does. 

The only functions that actually need to be run are the process_npz where all the mass files are processed and created into new mass files with additional information. It takes one input of a folder location where all files that end with ".npz" will be processed. It also takes an input name_indicator. Each mass file will be saved into a new file with the name name_indicator-mass file name. This name_indicator should be used in later functions to help identity the processed mass files 

The next function to be ran is the BBN_Abundances which generates the BBN abundacnes for all the mass files and generates the diagnostic plots. This function takes inputs mass file location, processed mass file location, the same name_indicator used in the previous function, and a new name indicator BBN_abun_indicator. The BBN Abundances for each mass file will be saved with the name BBN_abun_indicator-mass file name. 

The last function to be run will be element_plots which generates some essential plots of BBN Abundances versus lifetimes. This function takes inputs BBN_abun_location, mass file location, and the same BBN_abun_indicator used in the previosu function. 

Out of all the functions in this .py, these are the only ones that would need to be run to generate plots/diagnostics we need. 

Current plots as of 8/17 can be found in BBN_Graphs.ipynb
