# Code to generate results in the PBLA paper

Here we provide the necessary pieces of code to, alongside the files in the data folder, replicate the results found in the PBLA paper. 
The files here are:

* `tdc_ed.R` - code to generate Eichner Dietz method results for the Tristan Da Cunha data, and generate plots comparing the three methods
* `tdc_pbla_mcmc.c` - MCMC code for PBLA method with Tristan da Cunha data
* `tdc_pbla_lh.h` - PBLA likelihood for Tristan Da Cunha MCMC (required to run `tdc_pbla_mcmc.c`)

Code to generate deterministic Ebola data analysis, based on the method of 
> C.L. Althaus. Estimating the reproduction number of Ebola virus (EBOV) during
> the 2014 outbreak in West Africa. PLoS Currents Outbreaks, 6, 2014,

was adapted code made available at their Github repository, available at https://github.com/calthaus/Ebola
