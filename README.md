# DiffuseBunching

This repository contains computational code for computing patterns of diffuse bunching around kinks and notches. 

The directory `data` contains two CSV files with simulated data, one for bunching around a pure kink (`sample_kink.csv`), and the other around a notch (`sample_notch.csv`). 
Both simulations use the same bracket threshold ($k = 300$) and marginal tax rates ($t_0 = 0.1$ below the threshold and $t_1 = 0.2$ above the threshold). 

To run the estimation, run the script `matlab/main.m` in Matlab. 
It should estimate parameters for the elasticity ($e$), lumpiness parameter ($\mu$) and, in the case of the notch, the size of the notch ($dT$). 

The script will produce figures displaying both the underlying data (histogram) and the simulated density under the best-fit parameters. 
