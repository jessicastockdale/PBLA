//////////////////////////////////////////////////////////////////////////////
// This is MCMC code for the PBLA likelihood with the Tristan da Cunha data //
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
// GSL Libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

////////////////////////////////////////////////////////////////////////

// structure defined for each individual.
struct individual {
  int label;
  double r;
  int t; // 1 = infant, 2 = child, 3 = adult
};

// include the functions needed to compute log likelihood etc
#include "tdc_pbla_lh.h"

// Define the variable need to create random number generator
gsl_rng *r_seed;

////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

////////////////////////////////////////////////////////////////////////
    // get things set up
	int N = 254; // size of the population
	int iter = 100000 ; // number of MCMC iterations
	int burn = 10000; // burn in period
	int seed = time(NULL); // seed

	// define starting values for model parameters
	double beta1 = 0.00451, beta2 = 0.00181, beta3 = 0.00131, gamma = 0.371;

    // priors
	double mu_beta1 = 0.0000001, nu_beta1 = 0.00001;
	double mu_beta2 = 0.0000001, nu_beta2 = 0.00001;
	double mu_beta3 = 0.0000001, nu_beta3 = 0.00001;
	double mu_gamma = 0.0001, nu_gamma = 0.001;

    // tuning parameters
	double sigma_beta1 = 0.005;
	double sigma_beta2 = 0.005;
	double sigma_beta3 = 0.005;
	double sigma_gamma = 0.2;

	// Environment for the random number generator */
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r_seed = gsl_rng_alloc(T);
	gsl_rng_set(r_seed,seed);

	// define some job variables.
	int i, j, k;
	time_t t1, t2;

////////////////////////////////////////////////////////////////////////

    // Finally, read in the data
	FILE *F10 = fopen("../data/TristanDaCunha/tdc_jitteredtimes.txt", "r");
	gsl_vector *obs_data = gsl_vector_alloc(N);
	gsl_vector_fscanf(F10, obs_data);

	FILE *F11 = fopen("../data/TristanDaCunha/tdc_agegroups.txt", "r");
	gsl_vector *type_data = gsl_vector_alloc(N);
	gsl_vector_fscanf(F11,type_data);

	printf("The observed data have been read ... \n");

////////////////////////////////////////////////////////////////////////

    // calculate number of infectives nI
	int nI=0;
    for (i=0; i<N; i++){
        if (gsl_vector_get(obs_data, i) != -1000.0) {
            nI++;
        }
        else{break;}
    }
    printf("nI = %d \n", nI);

    // create a structure which will handle the data;
	struct individual data[N];
	for (i=0; i<N; i++) {

        //removal time
		data[i].r = gsl_vector_get(obs_data, i);
		if (data[i].r == -1000.0) data[i].r = GSL_POSINF;
		printf(" r=%lf \n",data[i].r );

        // type - infant=1, child=2 or adult=3
		data[i].t = gsl_vector_get(type_data, i);

        }

    //output files
    FILE *F1 = fopen("beta1.txt", "w");
    FILE *F2 = fopen("beta2.txt", "w");
    FILE *F3 = fopen("beta3.txt", "w");
    FILE *F4 = fopen("gamma.txt", "w");

////////////////////////////////////////////////////////////////////////


    // Check the lh looks sensible
    printf("lh = %lf \n", log_lh(data, beta1, beta2, beta3, gamma, N, nI));
    system("pause");

    // required values for MCMC
    double loglh_can, loglh_cur, log_denom, log_nomin, beta1_can, beta2_can, beta3_can, gamma_can;
    loglh_cur = log_lh(data, beta1, beta2, beta3, gamma, N, nI);
    double u;

////////////////////////////////////////////////////////////////////////

    /// MCMC loop starts here
	printf(" - Loop starts here ... \n");
	int it, rep;
	time(&t1);

	for (it = 0; it < (iter + burn); it++) {
	  printf("    - Percentage completed:"); printf(" %.0lf \r", 100 * ((double)it / ((double)iter + (double)burn) ));

    // update beta 1

    //propose a candidate value
		beta1_can = beta1 + gsl_ran_gaussian(r_seed, sigma_beta1);

		if (beta1_can > 0.0){
			// compute likelihood ratio;
			loglh_can = log_lh(data, beta1_can, beta2, beta3, gamma, N, nI);

            log_denom = loglh_cur + log(gamma_prior(beta1, mu_beta1, nu_beta1));
			log_nomin = loglh_can + log(gamma_prior(beta1_can, mu_beta1, nu_beta1));

			u = gsl_ran_flat(r_seed, 0.0, 1.0);
            // accept/reject:
			if (log(u) < log_nomin - log_denom) {
				loglh_cur = loglh_can;
				beta1 = beta1_can;
			}
		}


    // update beta 2

    //propose a candidate value
		beta2_can = beta2 + gsl_ran_gaussian(r_seed, sigma_beta2);

		if (beta2_can > 0.0) {
			// compute likelihood ratio;
			loglh_can = log_lh(data, beta1, beta2_can, beta3, gamma, N, nI);

            log_denom = loglh_cur + log(gamma_prior(beta2, mu_beta2, nu_beta2));
			log_nomin = loglh_can + log(gamma_prior(beta2_can, mu_beta2, nu_beta2));

			u = gsl_ran_flat(r_seed, 0.0, 1.0);
            // accept/reject:
			if (log(u) < log_nomin - log_denom) {
				loglh_cur = loglh_can;
				beta2 = beta2_can;
			}
		}


    //   update beta 3

    //propose a candidate value
		beta3_can = beta3 + gsl_ran_gaussian(r_seed, sigma_beta3);

		if (beta3_can > 0.0) {
			// compute likelihood ratio;
			loglh_can = log_lh(data, beta1, beta2, beta3_can, gamma, N, nI);

            log_denom = loglh_cur + log(gamma_prior(beta3, mu_beta3, nu_beta3));
			log_nomin = loglh_can + log(gamma_prior(beta3_can, mu_beta3, nu_beta3));

			u = gsl_ran_flat(r_seed, 0.0, 1.0);
            // accept/reject:
			if (log(u) < log_nomin - log_denom) {
				loglh_cur = loglh_can;
				beta3 = beta3_can;
			}
		}

    // update gamma

    //propose a candidate value
		gamma_can = gamma + gsl_ran_gaussian(r_seed, sigma_gamma);

		if (gamma_can > 0.0) {
			// compute likelihood ratio;
			loglh_can = log_lh(data, beta1, beta2, beta3, gamma_can, N, nI);

            log_denom = loglh_cur + log(gamma_prior(gamma, mu_gamma, nu_gamma));
			log_nomin = loglh_can + log(gamma_prior(gamma_can, mu_gamma, nu_gamma));

			u = gsl_ran_flat(r_seed, 0.0, 1.0);
            // accept/reject:
			if (log(u) < log_nomin - log_denom) {
				loglh_cur = loglh_can;
				gamma = gamma_can;
			}
		}


		/// write output to files.
		if (it>0){
		       fprintf(F1, "%lf \n", beta1);
		       fprintf(F2, "%lf \n", beta2);
		       fprintf(F3, "%lf \n", beta3);
               fprintf(F4, "%lf \n", gamma);
		}


	}
	time(&t2);


	printf("\n* Calculations complete ... Time needed: %.3lf min\n", difftime(t2,t1)/60.0);


	return 0;
}
