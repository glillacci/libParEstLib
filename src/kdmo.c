/*
 *  kdmo.c
 *
 *  This file contains the functions that implement the Kolmogorov distance
 *	matching optimization procedure.
 *
 *  Created by Gabriele Lillacci in April 2010.
 *	Latest revision: March 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
 */

#include "../parestlib.h"


/**
 This function computes the Kolmogorov distance between two data samples.
 For this version of the routine, the vectors data1 and data2 MUST BE SORTED IN ASCENDING ORDER!!
 This routine is adapted from Numerical Recipes Third Edition.
 */
double ksdist_two_sample_nr_presort (const gsl_vector * data1, const gsl_vector * data2)
{
	// Initialize necessary variables
	size_t j1 = 0, j2 = 0, n1 = data1->size, n2 = data2->size;
	double D = 0.0, d1, d2, dt, en1 = n1, en2 = n2, fn1 = 0.0, fn2 = 0.0;

	// Perform the computation
	while (j1 < n1 && j2 < n2)
	{
		if ((d1 = gsl_vector_get (data1, j1)) <= (d2 = gsl_vector_get (data2, j2)))
		{
			do
				fn1 = ++j1/en1;
			while (j1 < n1 && d1 == gsl_vector_get (data1, j1));
		}
		if (d2 <= d1)
		{
			do
				fn2 = ++j2/en2;
			while (j2 < n2 && d2 == gsl_vector_get (data2, j2));
		}
		if ((dt = fabs (fn2 - fn1)) > D)
			D = dt;
	}

	// Return
	return D;
}


/**
 This function computes the Kolmogorov distance between a data sample and an exact CDF.
 This routine is adapted from Numerical Recipes Third Edition.
 */
double ksdist_one_sample_nr (const gsl_vector * data, double (* exact)(), size_t count, ...)
{
	// Initialize necessary variables
	size_t j = 0, n = data->size;
	double D = 0.0, dt, en = (double) n, ff, fn, fo = 0.0;

	// Create a sorted copy of data
	gsl_vector * datas = gsl_vector_alloc (data->size);
	gsl_vector_memcpy (datas, data);
	gsl_sort_vector (datas);

	// Create a vector for the values of the exact cdf
	gsl_vector * ex = gsl_vector_alloc (datas->size);

	// Recover the distribution parameters from the variadic list
	double pars[count];

	va_list ap;
	va_start (ap, count);

	for	(size_t i = 0; i < count; i++) {
		pars[i] = va_arg (ap, double);
	}

	va_end (ap);

	// Evaluate the exact distribution for the data points contained in data. We have 3 cases.
	switch (count) {
	case 1:
		// The exact cdf is a one-parameter distribution
		for (size_t i = 0; i < ex->size; i++) {
			gsl_vector_set (ex, i, (* exact) (gsl_vector_get (datas, i), pars[0]));
		}
		break;
	case 2:
		// The exact cdf is a two-parameter distribution
		for (size_t i = 0; i < ex->size; i++) {
			gsl_vector_set (ex, i, (* exact) (gsl_vector_get (datas, i), pars[0], pars[1]));
		}
		break;
	case 3:
		// The exact cdf is a three-parameter distribution
		for (size_t i = 0; i < ex->size; i++) {
			gsl_vector_set (ex, i, (* exact) (gsl_vector_get (datas, i), pars[0], pars[1], pars[2]));
		}
		break;
	}

	// Compute the Kolmogorov distance
	for (j = 0; j < n; j++)
	{
		// Compute the data cdf
		fn = (j+1)/en;
		// Extract value of exact cdf for current data point
		ff = gsl_vector_get (ex, j);
		// Find the maximum distance for current data point
		dt = (fabs (fo - ff) > fabs (fn - ff)) ? fabs (fo - ff) : fabs (fn - ff);
		// Update global distance if necessary
		D = (dt > D) ? dt : D;
		// Update fo
		fo = fn;
	}

	// Free the manually allocated resources
	gsl_vector_free (datas);
	gsl_vector_free (ex);

	// Return
	return D;
}


/**
 This function computes the Kolmogorov distance between two data samples.
 This routine is adapted from Numerical Recipes Third Edition.
 */
double ksdist_two_sample_nr (const gsl_vector * data1, const gsl_vector * data2)
{
	// Initialize necessary variables
	size_t j1 = 0, j2 = 0, n1 = data1->size, n2 = data2->size;
	double D = 0.0, d1, d2, dt, en1 = n1, en2 = n2, fn1 = 0.0, fn2 = 0.0;

	// Create sorted copies of data1 and data2
	gsl_vector * data1s = gsl_vector_alloc (data1->size);
	gsl_vector * data2s = gsl_vector_alloc (data2->size);
	gsl_vector_memcpy (data1s, data1);
	gsl_vector_memcpy (data2s, data2);
	gsl_sort_vector (data1s);
	gsl_sort_vector (data2s);

	// Perform the computation
	while (j1 < n1 && j2 < n2)
	{
		if ((d1 = gsl_vector_get (data1s, j1)) <= (d2 = gsl_vector_get (data2s, j2)))
		{
			do
				fn1 = ++j1/en1;
			while (j1 < n1 && d1 == gsl_vector_get (data1s, j1));
		}
		if (d2 <= d1)
		{
			do
				fn2 = ++j2/en2;
			while (j2 < n2 && d2 == gsl_vector_get (data2s, j2));
		}
		if ((dt = fabs (fn2 - fn1)) > D)
			D = dt;
	}

	// Free the manually allocated resources
	gsl_vector_free (data1s);
	gsl_vector_free (data2s);

	// Return
	return D;
}


/**
 This function computes the Kolmogorov distance between two data samples.
 Here data2 is treated as the "reference" sample, while data1
 represents the "estimated" sample.
 */
double ksdist_two_sample_2 (const gsl_vector * data1, const gsl_vector * data2)
{
	// Create the empirical cdfs of the two samples
	ecdf * ecdf1 = ecdf_create_from_data (data1);
	ecdf * ecdf2 = ecdf_create_from_data (data2);

	// Initialize computation of the Kolmogorov distance
	double D = 0, A = 0, B = 0, C = 0;
	size_t N = ecdf1->N;

	// Scan the data vector
	for (size_t	i = 1; i < N; i++) {
		// Compute distance for current data point
		A = fabs (ecdf_eval (ecdf2, gsl_vector_get (ecdf1->xdata, i)) - gsl_vector_get (ecdf1->ydata, i-1));
		B = fabs (gsl_vector_get (ecdf1->ydata, i) - ecdf_eval (ecdf2, gsl_vector_get (ecdf1->xdata, i)));
		C = (A >= B ? A : B);
		// Update overall distance
		D = (D >= C ? D : C);
	}

	/*
	// Scan the data1 vector
	for (size_t i = 0; i < N; i++)
	{
		// Compute distance for current point
		A = fabs (ecdf_eval (ecdf2, gsl_vector_get (ecdf1->xdata, i)) - gsl_vector_get (ecdf1->ydata, i));
		// Update distance
		D = (D < A ? A : D);
	}
	*/

	// Destroy the two ecdfs
	ecdf_free (ecdf1);
	ecdf_free (ecdf2);

	// Return
	return D;
}


/**
 This function computes the Kolmogorov distance between two univariate empirical cumulative
 distribution functions (ecdf). Here ecdf2 is treated as the "reference" cdf, while ecdf1
 represents the "estimated" cdf.
 */
double ksdist_two_sample (const ecdf * ecdf1, const ecdf * ecdf2)
{
	// Initialize computation of the Kolmogorov distance
	double D = 0, A = 0, B = 0, C = 0;
	size_t N = ecdf1->N;

	// Scan the data vector
	for (size_t	i = 1; i < N; i++) {
		// Compute distance for current data point
		A = fabs (ecdf_eval (ecdf2, gsl_vector_get (ecdf1->xdata, i)) - gsl_vector_get (ecdf1->ydata, i-1));
		B = fabs (gsl_vector_get (ecdf1->ydata, i) - ecdf_eval (ecdf2, gsl_vector_get (ecdf1->xdata, i)));
		C = (A >= B ? A : B);
		// Update overall distance
		D = (D >= C ? D : C);
	}

	// Return
	return D;
}


/**
 This function computes the Kolmogorov distance of the empirical cumulative distribution function
 of the univariate data points in data from an exact cdf.
 */
double ksdist (const gsl_vector * data, double (* exact)(), size_t count, ...)
{
	// Create a copy of the data vector and sort it
	gsl_vector * data2 = gsl_vector_alloc (data->size);
	gsl_vector_memcpy (data2, data);
	gsl_sort_vector (data2);

	// Recover the distribution parameters from the variadic list
	double pars[count];

	va_list ap;
	va_start (ap, count);

	for	(size_t i = 0; i < count; i++) {
		pars[i] = va_arg (ap, double);
	}

	va_end (ap);

	// Generate the exact distribution data points. We have 3 (or more...) cases.
	gsl_vector * ex = gsl_vector_alloc (data2->size);

	switch (count) {
		case 1:
			// The exact cdf is a one-parameter distribution
			for (size_t i = 0; i < ex->size; i++) {
				gsl_vector_set (ex, i, (* exact) (gsl_vector_get (data2, i), pars[0]));
			}
			break;
		case 2:
			// The exact cdf is a two-parameter distribution
			for (size_t i = 0; i < ex->size; i++) {
				gsl_vector_set (ex, i, (* exact) (gsl_vector_get (data2, i), pars[0], pars[1]));
			}
			break;
		case 3:
			// The exact cdf is a three-parameter distribution
			for (size_t i = 0; i < ex->size; i++) {
				gsl_vector_set (ex, i, (* exact) (gsl_vector_get (data2, i), pars[0], pars[1], pars[2]));
			}
			break;
	}

	// Initialize computation of the Kolmogorov distance
	double D = 0, A = 0, B = 0, C = 0;
	size_t N = ex -> size;

	// Scan the data vector
	for (size_t	i = 1; i <= N; i++) {
		// Compute distance for current data point
		A = gsl_vector_get (ex, i-1) - ((double) (i-1)/N);
		B = ((double) i/N) - gsl_vector_get (ex, i-1);
		C = (A >= B ? A : B);
		// Update overall distance
		D = (D >= C ? D : C);
	}

	// Free all the manually allocated resources
	gsl_vector_free (data2);
	gsl_vector_free (ex);

	// Return
	return D;
}


double kscost (const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M,
			   const gsl_matrix * exact, int (* ptrToFunc) (double, const double[], double[], void *),
			   int (* ptrToJac) (double, const double[], double *, double[], void *),
			   int (* ptrToOutput) (double, const double[], double[], void *), void * params)
{
	// Detect the dimensions of the problem
	size_t S = time->size;					// Number of samples
	size_t N = M->size1;					// Number of outputs
	size_t DIM = xZero->size + theta->size;	// Dimension of the system

	// *** 1 *** Evaluate the "expected measurement" by solving the ODE system

	// Allocate the expected measurement matrix
	gsl_matrix * M2 = gsl_matrix_alloc (M->size1, M->size2);

	// Set up ODE solve environment
	gsl_odeiv_system sys = {ptrToFunc, ptrToJac, DIM, params};
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (PRECISION, 0.0);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM);

	// Set initial time and initial step size
	double t = gsl_vector_get (time, 0);
	double k = INITSTEP;

	// Declare array y and initialize to initial conditions
	double y[DIM];
	int i;
	for (i = 0; i < xZero->size; i++)
	{
		y[i] = gsl_vector_get (xZero, i);
	}
	for (i = 0; i < theta->size; i++)
	{
		y[i+xZero->size] = gsl_vector_get (theta, i);
	}

	// Declare array for the "output" to store the results for each desired time step
	double h[N];

	for (i = 0; i < S; i++)
	{
		// Loop to obtain a result for every specified time step
		while (t < gsl_vector_get (time, i))
		{
			// Evolve the system to the next time step
			int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t,  gsl_vector_get (time, i), &k, y);

			// Break loop if an error occurs
			if (status != GSL_SUCCESS)
		        break;
		}
		// Evaluate the output of the system at the current time step
		ptrToOutput (t, y, h, params);

		// Store the computed output in appropriate column of M2
		gsl_vector_view m2v = gsl_matrix_column (M2, i);
		gsl_vector_set_array (&m2v.vector, h);
	}
	// At this point every row of M2 contains one expected output


	// *** 2 *** Compute the costs in terms of means and variances

	// Compute the residuals, the result is stored in M2
	gsl_matrix_sub (M2, M);
	gsl_matrix_scale (M2, -1.0);
	double dist = 0;

	for (i = 0; i < N; i++)
	{
		gsl_vector_const_view data = gsl_matrix_const_row (M2, i);
		gsl_vector_const_view ex = gsl_matrix_const_row (exact, i);

		dist += ksdist_two_sample_2 (&data.vector, &ex.vector);
	}

	// Free all the manually allocated resources
	gsl_matrix_free (M2);
	gsl_odeiv_step_free (s);
	gsl_odeiv_control_free (c);
	gsl_odeiv_evolve_free (e);

	return dist;
}
