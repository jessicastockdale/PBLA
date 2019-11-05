# Code to generate results in the PBLA paper

Here we provide the necessary pieces of code to, alongside the files in the data folder, replicate the results found in the PBLA paper. 
The files here are:

* `tdc_ed.R` - code to generate Eichner Dietz method results for the Tristan Da Cunha data, and create plots comparing the three methods
* `tdc_pbla_mcmc.c` - MCMC code for PBLA method with Tristan da Cunha data
* `tdc_pbla_lh.h` - PBLA likelihood for Tristan Da Cunha MCMC (required to run `tdc_pbla_mcmc.c`)
* `ebola_pbla_and_althaus.R` - code to generate PBLA and Althaus analysis for Ebola virus data, and create plots comparing the analyses
* `bij_ed_MAPs.R` - code to generate Eichner Dietz method MAPs for a heterogenous mixing model such as in the FMD analysis
* `bij_PBLA.R` - code to generate PBLA method MAPs and profile log likelihoods for a heterogenous mixing model such as in the FMD analysis
* `bij_results.R` - code to generate heterogenous mixing model plots, comparing PBLA/ED/DA-MCMC
* `lh_foropt.cpp` - PBLA log likelihood, for heterogenous mixing model analysis
* `lh_foropt_post.cpp` - PBLA posterior, for heterogenous mixing model analysis

The last 5 files are for analysis using a heterogeneous mixing model such as that used in the Foot and Mouth outbreak example. The code is set up for the same model as the FMD analysis: infection rate beta_ij from infector i to infectee j as a function of the Euclidean distance between them and the numbers of cows and sheep on each farm. However, the code may be adapted for any model desired (remember that each parameter introduced will greatly increase the time taken to find the MLE etc.).

Code to generate deterministic Ebola data analysis, based on the method of 
> C.L. Althaus. Estimating the reproduction number of Ebola virus (EBOV) during
> the 2014 outbreak in West Africa. PLoS Currents Outbreaks, 6, 2014,

was adapted from code made available at their Github repository, https://github.com/calthaus/Ebola


