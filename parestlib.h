/*
 *  parestlib.h
 *  ParEstLib
 *
 *  Header file for ParEstLib
 *
 *  Created by Gabriele Lillacci in April 2009.
 *	Latest revision: October 2011.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA. 
 */

#ifndef _PARESTLIB_H_
#define _PARESTLIB_H_


/*
 System libraries includes
 */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>


/* 
 GSL library includes
 */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_roots.h>


/* 
 Custom library includes
 */
#include <stochmod.h>

/*
 Macro definitions
 */

#define EPS 2.2204e-16
#define PRECISION EPS
#define INITSTEP EPS
#define MAXTIME 3


/*
 Data type definitions
 */

// This structure defines a Lyapunov equation
typedef struct {
	gsl_matrix * A;
	gsl_matrix * Q;
	size_t DIM;
} LyapunovEquation;

// This structure defines an empirical cumulative distribution function
typedef struct {
	gsl_vector * xdata;
	gsl_vector * ydata;
	size_t N;
} ecdf;

// This structure defines a discretely observed sample path of a chemical reaction network
typedef struct {
	size_t nspecies;
	size_t ntimes;
	gsl_vector * times;
	gsl_matrix * counts;
} rxn_sample_path;

// This structure defines a discretely observed ensemble of trajectories of a chemical reaction network
typedef struct {
	size_t nreplic;
	size_t nspecies;
	size_t ntimes;
	rxn_sample_path ** data;
} rxn_ensemble;


/*
 Exported functions prototype declarations == KALMANFILTER.C
 */
int hekf_time_update_x (double tCurrent, double tNext, const gsl_vector * xPlus, gsl_vector * xMinus,
						int (* ptrToFunc) (double, const double[], double[], void *),
						int (* ptrToJac) (double, const double[], double *, double[], void *), void * params);
int hekf_time_update_p (double tCurrent, double tNext, const gsl_matrix * pPlus, gsl_matrix * pMinus, LyapunovEquation * eqn);
int hekf_compute_gain (gsl_matrix * L, const gsl_matrix * pMinus, const gsl_matrix * H, const gsl_matrix * R);
int hekf_measurement_update_x (double tCurrent, const gsl_vector * xMinus, gsl_vector * xPlus, const gsl_matrix * L, 
							   const gsl_vector * mCurrent, int (* ptrToOutput) (double, const double[], double[], void *), void * params);
int hekf_measurement_update_p (gsl_matrix * pPlus, const gsl_matrix * pMinus, const gsl_matrix * L,
							   const gsl_matrix * H, const gsl_matrix * R);


/*
 Exported functions prototype declarations == UTILITY.C
 */
int easy_matrix_product (gsl_matrix *, const gsl_matrix *, const gsl_matrix *);
int easy_matrix_inverse (gsl_matrix *, const gsl_matrix *);
int matlab_semicolon (gsl_matrix *, double, double);
int concatenate_arrays (double[], const double[], int, const double[], int);
int print_matrix (double[], int, int);
int print_gsl_matrix (const gsl_matrix *);
int print_gsl_vector (const gsl_vector * v);
int print_gsl_subvector (const gsl_vector * v, size_t start, size_t nel);
int row_number (int, int, int);
int col_number (int, int, int);
double mean (gsl_matrix * mat);
int gsl_matrix_horzcat (gsl_matrix * cat, const gsl_matrix * mat1, const gsl_matrix * mat2);
int gsl_matrix_vertcat (gsl_matrix * cat, const gsl_matrix * mat1, const gsl_matrix * mat2);
void gsl_vector_set_array (gsl_vector * v, double elem[]);
void gsl_matrix_set_array (gsl_matrix * mat, double elem[]);
gsl_vector * gsl_vector_realloc (gsl_vector * oldvec, size_t newsize);


/*
 Exported functions prototype declarations == MEASUREMENT.C
 */
int generate_measurements (int (* ptrToFunc) (double, const double[], double[], void *),
						   int (* ptrToJac) (double, const double[], double *, double[], void *),
						   int (* ptrToOutput) (double, const double[], double[], void *), 
						   const gsl_matrix * time, gsl_matrix * sig, const gsl_matrix * R, const gsl_matrix * xZero, void * params);


/*
 Exported functions prototype declarations == VARIANCE.C
 */
int chi2test (gsl_vector * vhp, gsl_matrix * vhi,  const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M,
			  const gsl_matrix * R, double p, int (* ptrToFunc) (double, const double[], double[], void *),
			  int (* ptrToJac) (double, const double[], double *, double[], void *), int (* ptrToOutput) (double, const double[], double[], void *),
			  void * params);
int chi2cost (gsl_vector * ave, gsl_vector * var,  const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M,
			  const gsl_matrix * R, int (* ptrToFunc) (double, const double[], double[], void *),
			  int (* ptrToJac) (double, const double[], double *, double[], void *), int (* ptrToOutput) (double, const double[], double[], void *),
			  void * params);


/*
 Exported functions prototype declarations == MINIM.C
 */


/*
 Exported functions prototype declarations == DERIV.C
 */
int central_diff_sv (double * D, double * E, double (* func) (const gsl_vector *, void *), size_t i, const gsl_vector * x, double h, void * params);
int adaptive_gradient_sv (gsl_vector * df, gsl_vector * err, const gsl_vector * x, double (* func) (const gsl_vector *, void *), double tol, void * params);


/*
 Exported functions prototype declarations == ECDF.C
 */
void ecdf_free (ecdf * dist);
ecdf * ecdf_create_from_data (const gsl_vector * data);
double ecdf_eval (const ecdf * dist, double x);
double ecdf_sample (const ecdf * dist, const gsl_rng * r);


/*
 Exported functions prototype declarations == KDMO.C
 */
double ksdist_two_sample_nr_presort (const gsl_vector * data1, const gsl_vector * data2);
double ksdist_one_sample_nr (const gsl_vector * data, double (* exact)(), size_t count, ...);
double ksdist_two_sample_nr (const gsl_vector * data1, const gsl_vector * data2);
double ksdist_two_sample_2 (const gsl_vector * data1, const gsl_vector * data2);
double ksdist_two_sample (const ecdf * ecdf1, const ecdf * ecdf2);
double ksdist (const gsl_vector * data, double (* exact)(), size_t count, ...);
double kscost (const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M, 
			   const gsl_matrix * exact, int (* ptrToFunc) (double, const double[], double[], void *),
			   int (* ptrToJac) (double, const double[], double *, double[], void *),
			   int (* ptrToOutput) (double, const double[], double[], void *), void * params);


/*
 Exported functions prototype declarations == STOCHSIM.C
 */
int ssa_direct_evolve (stochmod * model, const gsl_vector * params, gsl_vector * X, double * t, const gsl_rng * r);
rxn_sample_path * rxn_sample_path_alloc (size_t nspecies, size_t ntimes);
void rxn_sample_path_free (rxn_sample_path * rsp);
int ssa_direct_trajectory (stochmod * model, const gsl_vector * params, gsl_vector * times, gsl_vector * X0,
						   rxn_sample_path * path, const gsl_rng * r);
rxn_ensemble * rxn_ensemble_alloc (size_t nreplic, size_t nspecies, size_t ntimes);
void rxn_ensemble_free (rxn_ensemble * ens);
int ssa_direct_ensemble (stochmod * model, const gsl_vector * params, gsl_vector * times, const gsl_vector * X0,
						 rxn_ensemble * ens, const gsl_rng * r);
int rxn_ensemble_counts (const rxn_ensemble * ens, gsl_matrix * counts, const gsl_matrix * out);


/*
 Exported functions prototype declarations == DCMS.C
 */
double kolmogorov_cdf_marsaglia (int n, double d);
double kolmogorov_cdf_inverse (int M, double q);
int kolmogorov_cdf_inverse_M (double epsilon, double q);
int critical_S_DKW (double epsilon, double alpha, int M);


/*
 Exported functions prototype declarations == FCS.C
 */
int fcs_seek_to_keyword (FILE * data, const char * keyword);
int fread_floats_swap (float * f, size_t count, FILE * data);
int fcs3_read_int_kw (FILE * data, const char * keyword);
int fcs2_read_int_kw (FILE * data, const char * keyword);
int fread_floats_noswap (float * f, size_t count, FILE * data);


#endif
