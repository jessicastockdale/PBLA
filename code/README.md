# Code to generate results in the PBLA paper

Here we provide the necessary pieces of code to, alongside the files in the data folder, replicate the results found in the PBLA paper. 
The files here are:

* `tdc_ed.R` - code to generate Eichner Dietz method results for the Tristan Da Cunha data, and create plots comparing the three methods
* `tdc_pbla_mcmc.c` - MCMC code for PBLA method with Tristan da Cunha data
* `tdc_pbla_lh.h` - PBLA likelihood for Tristan Da Cunha MCMC (required to run `tdc_pbla_mcmc.c`)
* `ebola_pbla_and_althaus.R` - code to generate PBLA and Althaus analysis for Ebola virus data, and create plots comparing the analyses
* `fam_ed_MAPs.R` - code to generate Eichner Dietz method MAPs for (simulated) Foot and Mouth outbreak data
* `fam_PBLA.R` - code to generate PBLA method MAPs and profile log likelihoods for (simulated) Foot and Mouth outbreak data
* `lh_foropt.cpp` - PBLA log likelihood, for Foot and Mouth analysis
* `lh_foropt_post.cpp` - PBLA posterior, for Foot and Mouth analysis

Code to generate deterministic Ebola data analysis, based on the method of 
> C.L. Althaus. Estimating the reproduction number of Ebola virus (EBOV) during
> the 2014 outbreak in West Africa. PLoS Currents Outbreaks, 6, 2014,

was adapted from code made available at their Github repository, https://github.com/calthaus/Ebola


