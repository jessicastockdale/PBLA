/// PBLA log likelihood function, for TdC data analysis ///
// susceptible type (age group) determines infection rate beta
// assume the initial infective is the initial removal


double log_lh(struct individual data[], double beta1, double beta2, double beta3, double g, int N, int n) {

    // set up
    double loglh=0.0;
    double denom, chiphi, psi, b;
    int j, i;
    double B=16*beta1 + 30*beta2 + 168*beta3; // just set the number in each group manually here
    double d=g+B;

    // All parameters must be non-negative
    if (beta1<0 || beta2<0 || beta3<0 || g<0){loglh = -10000000;}
    else{

    // Now calculate the log ikelihood:
    for (j=1; j<n; j++){
            chiphi = 0.0;
            psi = 1.0;

             // find relevant beta
            if (data[j].t==1){b = beta1;}
            else if (data[j].t==2){b = beta2;}
            else{b = beta3;}

            for (i=0; i<n; i++){
                    if (i!=j){
                            if (data[j].r<data[i].r){
                                denom = (1 - (b/(2*(d+b)))*exp(-d*(data[i].r-data[j].r)));
                                psi = psi*(1 - (b/(2*(d+b)))*exp(-d*(data[i].r-data[j].r)));
                            }
                            else{
                                denom = ( d/(d+b) + (b/(2*(d+b)))*exp(-d*(data[j].r-data[i].r)));
                                psi = psi*( d/(d+b) + (b/(2*(d+b)))*exp(-d*(data[j].r-data[i].r)));
                            }
                        chiphi = chiphi + 0.5*exp(-d*fabs(data[i].r-data[j].r))*d/((b+d)*denom);
                    }
            }
            chiphi =  chiphi*b;
            loglh = loglh + log(chiphi) + log(psi);

    }

    loglh = loglh + n*(log(g) - log(d));
    }
    return loglh;
}


/// Prior function ///
// define a gamma prior function which returns the density;
double gamma_prior(double value, double par1, double par2) {
	double f = gsl_ran_gamma_pdf(value, par1, 1.0/par2);
	return f;
}


