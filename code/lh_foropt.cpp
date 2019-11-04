// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include <RcppGSL.h>
//using namespace Rcpp;

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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cdf.h>

/// Pochhammer symbol function
double poch(int x, int k){
    double prod, j;
    if (k==0){
        prod = 1;
    }
    else if (k==1){
        prod = x;
    }
    else{
        prod = x;
        for (j = 1; j < k; j++){
              prod = prod*(x-j);
        }
    }
    return prod;
}

// structure defined for each individual.
struct individual {
  double r;
};

/// Likelihood function for PBLA Foot and Mouth, spatial beta ///
// b[i,j] = b0*K(i,j)*(ep(nci^ze)+nsi^ze)*(xi(ncj^ze)+nsj^ze)
// K(i,j) = spatial kernel
// nci = #cows on farm i
// nsi = #sheep on farm i


// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
 double log_lh(const RcppGSL::Vector &params, const RcppGSL::Vector &R2, RcppGSL::Matrix &pmat, int N, int n, int m, const RcppGSL::Vector &n_c, const RcppGSL::Vector &n_s) {

    // define required variables
    double lh=0.0;
    double denom, chij, psij, value, sum, exp1, exp2;
    int j, i, k, l, p;

    /// params to be estimated
    double beta0 = gsl_vector_get(params,0);
    double g = gsl_vector_get(params,1);
    double v = gsl_vector_get(params,2);
    double ep = gsl_vector_get(params,3);
    double xi = gsl_vector_get(params,4);
    double ze = gsl_vector_get(params,5);


    /// beta matrix
    gsl_matrix *b = gsl_matrix_alloc(N, N);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            gsl_matrix_set(b,i,j,v/(gsl_matrix_get(pmat,i,j) + pow(v,2))*(ep*(pow(gsl_vector_get(n_c,i),ze))+
                    pow(gsl_vector_get(n_s,i),ze))*(xi*pow((gsl_vector_get(n_c,j)),ze)+
                    pow(gsl_vector_get(n_s,j),ze)));
        }
    }
    /// beta matrix - multiply through by beta0
    gsl_matrix_scale(b, beta0);

    /// B
    gsl_vector *B = gsl_vector_calloc(n);
    for (i=0; i<n; i++) {
            for (j=n; j<N; j++) {
                gsl_vector_set(B, i, gsl_vector_get(B,i)+gsl_matrix_get(b,i,j) );
            }
    }

    /// delta
    gsl_vector *d = gsl_vector_alloc(n);
    for (i=0; i<n; i++) {
        gsl_vector_set(d, i, g + gsl_vector_get(B,i) );
    }

  /// begin lh calculation
    if (gsl_matrix_min (b)<0 || g<0){lh = -10000000;}
    else{
    // don't need sum term for initial infective- just don't do loop for first infective
    for (j = 1; j < n; j++){
            chij = 0;
            psij=1;
            // calculate chi(j) term and psi(j) together
            for (i = 0; i < n; i++){
                if (i != j){
                    sum=0.0;
                    denom=0.0;
                    if (R2->data[j] < R2->data[i]){
                        for (l = 0; l < m; l++){
                          exp1=0.0;
                          for (p = 0; p < (l+1); p++){
                                if (m>1){
                                      exp1 = exp1 + gsl_sf_choose(l, p)*(pow(R2->data[i]-R2->data[j],l-p))*poch(m-1+p,p)/(pow(gsl_vector_get(d,j)+gsl_vector_get(d,i),p));
                                 }else{exp1 = 1;}
                          }
                          denom = denom + exp1*((pow(gsl_vector_get(d,i),m))*pow(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j),l-m) - pow(gsl_vector_get(d,i),l))/gsl_sf_fact(l);
                          sum = sum + exp1*pow(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j),l)/gsl_sf_fact(l);
                            }
                        sum = sum*pow(gsl_vector_get(d,i)/(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j)),m)*exp(-gsl_vector_get(d,i)*(R2->data[i]-R2->data[j]))*pow(gsl_vector_get(d,j)/(gsl_vector_get(d,j) + gsl_vector_get(d,i)),m);
                        denom  = 1 + denom*exp(-gsl_vector_get(d,i)*(R2->data[i]-R2->data[j]))*pow(gsl_vector_get(d,j)/(gsl_vector_get(d,j) + gsl_vector_get(d,i)),m);
                        psij = psij*(denom);
                    }
                    else{
                        for (l = 0; l < m; l++){
                          exp2=0.0;
                          for (p = 0; p < m; p++){
                                if (m>1){
                                      exp2 = exp2 + gsl_sf_choose(m-1, p)*(pow(R2->data[j]-R2->data[i],m-1-p))*poch(l+p,p)/(pow(gsl_vector_get(d,j)+gsl_vector_get(d,i),p));
                                }else{exp2 = 1;}
                          }
                          denom = denom + exp2*(pow(gsl_vector_get(d,i)/(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j)),m-l)-1)*pow(gsl_vector_get(d,i),l)/pow(gsl_vector_get(d,i)+gsl_vector_get(d,j),l+1);
                          sum = sum + exp2*pow(gsl_vector_get(d,i)/(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j)),m)*pow(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j),l)/pow(gsl_vector_get(d,i)+gsl_vector_get(d,j),l+1);
                          }
                        denom  = 1 - (1-pow(gsl_vector_get(d,i)/(gsl_vector_get(d,i)+gsl_matrix_get(b,i,j)),m))*gsl_cdf_gamma_P(R2->data[j]-R2->data[i],(double)m,1/gsl_vector_get(d,j)) + denom*pow(gsl_vector_get(d,j),m)*
                                exp(-gsl_vector_get(d,j)*(R2->data[j]-R2->data[i]))/gsl_sf_fact(m-1);
                        sum = sum*exp(-gsl_vector_get(d,j)*(R2->data[j]-R2->data[i]))*pow(gsl_vector_get(d,j),m)/gsl_sf_fact(m-1);
                        psij = psij*(denom);
                    }
                    chij = chij + (sum/denom)*gsl_matrix_get(b,i,j);
                }

            }

            lh = lh + log(chij) + log(psij);
        }
        // change of variables term
        for  (j = 0; j < n; j++){
            lh = lh + m*(log(g)-log(gsl_vector_get(d,j)));
        }

    }

    gsl_matrix_free(b);
    gsl_vector_free(B);
    gsl_vector_free(d);
    gsl_matrix_free(pmat);


    return lh;

}






