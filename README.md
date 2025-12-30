# Monte-Carlo-Molecular-Modelling
This repository contains the codes developed for the Monte Carlo part of the course on Molecular Modelling of the Master in Complex Systems and Biophysics at UB.

## simulacion_codes
This folder contains all the codes needed to simulate a 2D Ising model on a square lattice of size L.  
First of all the initial configuartion of spins must be created with the **configuration.f90** code:
- Inputs: L, file_name (to generate), random_seed (used: 123)
- Outputs: File containing the initial configuration of spins
  
After that, the **MC_simulation.f90** code is used to simulate the system: 
- Inputs: L, T, NMCS (# of Monte Carlo steps), NMEAS (every how many steps E and M are saved, and random_seed), file_name (containing the initial spin configuration)
- Outputs: File containining the $E$ and $M$ per particle time series

More detailed information about how to compile and execute them can be found in the codes as a comment.

## analysis_and_plots_codes
This folder contains all the codes used to analyse and plot the generated data. **appendix_T_2_1.py** was used to plot the time series in lin-log scale with $T=2.1$.
**binning_code.ipynb** was used for the rest of the plots. It contains the following functions to help the creation of the binning plots, the time series plots and the corresponding fit to obtain the statistical errors and the autocorrelation time. 
- binning
- detect_plateau
- model_fixed_plateau
- get_fit_params
- plot_binning_with_fit
- make_combined_plot

## results2DIssing.ods
Open office document containing the values of $\langle E\rangle$ and $\langle |M|\rangle$ per particle for each $L$ and $T$.
