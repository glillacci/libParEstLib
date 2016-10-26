/*
 *  dcms.c
 *
 *  This file contains the functions related to the
 *  distribution comparison model selection (DCMS) method.
 *
 *  Created by Gabriele Lillacci in November 2011.
 *	Latest revision: November 2011.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA. 
 */


#include <parestlib.h>


/**
 Supporting function for kolmogorov_cdf_marsaglia, not to be used directly.
 */
void mMultiply (double *A, double *B, double *C, int m)
{
    int i,j,k; double s;
    for(i=0;i<m;i++) for(j=0; j<m; j++)
    {s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s;}
}


/**
 Supporting function for kolmogorov_cdf_marsaglia, not to be used directly.
 */
void mPower (double *A, int eA, double *V, int *eV, int m, int n)
{
    double *B;int eB,i;
    if(n==1) {for(i=0;i<m*m;i++) V[i]=A[i];*eV=eA; return;}
    mPower(A,eA,V,eV,m,n/2);
    B=(double*)malloc((m*m)*sizeof(double));
    mMultiply(V,V,B,m); eB=2*(*eV);
    if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
    else {mMultiply(A,B,V,m); *eV=eA+eB;}
    if(V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;}
    free(B);
}


/**
 This function evaluates the cumulative distribution function of the Kolmorogov distribution
 using the numerical procedure described in Marsaglia 2003.
 */
double kolmogorov_cdf_marsaglia (int n, double d)
{ 
    int k,m,i,j,g,eH,eQ;
    double h,s,*H,*Q;
    
    //OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
    s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
    
    k=(int)(n*d)+1; m=2*k-1; h=k-n*d;
    H=(double*)malloc((m*m)*sizeof(double));
    Q=(double*)malloc((m*m)*sizeof(double));
    for(i=0;i<m;i++) for(j=0;j<m;j++)
        if(i-j+1<0) H[i*m+j]=0; else H[i*m+j]=1;
    for(i=0;i<m;i++) {H[i*m]-=pow(h,i+1); H[(m-1)*m+i]-=pow(h,(m-i));}
    H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
    for(i=0;i<m;i++) for(j=0;j<m;j++)
            if(i-j+1>0) for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
    eH=0; mPower(H,eH,Q,&eQ,m,n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++) {s=s*i/n; if(s<1e-140) {s*=1e140; eQ-=140;}}
    s*=pow(10.,eQ); free(H); free(Q); return s;
}

/**
 Supporting function for kolmogorov_cdf_inverse, not to be used directly.
 */
double kolmogorov_cdf_roots (double x, void * params)
{
    double * pars = (double *) params;
    return kolmogorov_cdf_marsaglia ((int) pars[0], x) - pars[1];
}


/**
 This function inverts the Kolmogorov distribution using the Brent-Dekker method.
 */
double kolmogorov_cdf_inverse (int M, double q)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0.5;
    double x_lo = 0.0, x_hi = 1.0;
    gsl_function F;
    
    // Set initial guesses using the DKW bounds
    x_hi = sqrt(-1/(2*(double)M)*log((1-q)/2));
    
    F.function = &kolmogorov_cdf_roots;
    double pars[2];
    pars[0] = (double) M;
    pars[1] = q;
    F.params = pars;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    
 
    printf ("\n\nusing %s method\n", 
            gsl_root_fsolver_name (s));
    
    printf ("%5s [%9s, %9s] %9s %9s\n",
            "iter", "lower", "upper", "root", 
            "err(est)");
    
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.001);
        if (status == GSL_SUCCESS)
            printf ("Converged:\n");
        
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                iter, x_lo, x_hi,
                r, 
                x_hi - x_lo);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free (s);
    
    return r;
}


/**
 This function inverts the Kolmogorov distribution with respect to the number
 of samples M.
 */
int kolmogorov_cdf_inverse_M (double epsilon, double q)
{
	// Set initial guess for M using the DKW bounds
    int M = ceil (-1/(2*epsilon*epsilon)*log((1-q)/2));
    double kM = -1.0, kM1 = 1.0;
    
    do
    {        
        // Evaluate the Kolmo cdf corresponding to M and M-1
        kM = kolmogorov_cdf_marsaglia (M, epsilon)-q;
        kM1 = kolmogorov_cdf_marsaglia (M+1, epsilon)-q;
        
        printf("M = %d -- kM = %f -- kM1 = %f\n", M, kM, kM1);
        
        M--;
    }
    while (kM>=0 || kM1<=0);
    
    return M+2;
}


/**
 This function evaluates the bound on the critical number of simulations S obtained
 using the DKW inequality
 */
int critical_S_DKW (double epsilon, double alpha, int M)
{
	return ceil (-log (alpha/2)/(2*pow (epsilon - sqrt (((double) -1/(2*M))*log (alpha/2)), 2)));
}
