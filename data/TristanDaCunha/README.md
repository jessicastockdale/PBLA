# Tristan Da Cunha respiratory disease data

This folder contains two data files:
* `tdc_jitteredtimes.txt` - removal times of the 254 individuals in the dataset. The first 40 individuals are those who were 
    infected, and non-infected individuals are given a removal time of -1000. Times have been jittered to avoid numerical error
    in the PBLA algorithm.
* `tdc_agegroups.txt` - age grouping of the 254 individuals in the dataset, ordered as in the removal times file. 1 = infants, 
    2 = children and 3 = adults.
