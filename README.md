# MSCR
Code used to generate the results in "Incorporating Memory into Continuous-Time Spatially Explicit Capture-Recapture Models" by Clara Panchaud, Ruth King, David Borchers, Hannah Worthington, Ian Durbach and Paul van Dam-Bates.

The folder "Functions" contains all of the functions to define both the standard continous-time SCR model and the extended memory MSCR version. 

The folder "Simulations" contains two subfolders, "Simulations_MSCR" and "Simulations_OU". The former contains code to simulate data from the MSCR model and fit it with both MSCR and SCR, while the latter does the same with data simulated from an OU process. Both also contain a file to calculate and plot the AC PDFs of individuals. 

The folder "American Martens" contains data from an American martens study, the analysis of this data set by both models and a file that plots the AC PDFs.
