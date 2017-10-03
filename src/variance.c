/*
 *  variance.c
 *
 *  This file contains the functions that implement the chi2 variance test and the
 *	computation of the costs for the moment matching optimization.
 *
 *  This file is part of libParEstLib.
 *  Copyright 2011-2017 Gabriele Lillacci.
 *
 *  libParEstLib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  libParEstLib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libParEstLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../parestlib.h"

int chi2test (gsl_vector * vhp, gsl_matrix * vhi,  const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M,
			  const gsl_matrix * R, double p, int (* ptrToFunc) (double, const double[], double[], void *),
			  int (* ptrToJac) (double, const double[], double *, double[], void *), int (* ptrToOutput) (double, const double[], double[], void *),
			  void * params)
{
	// Detect the dimensions of the problem
	size_t S = time->size;					// Number of samples
	size_t N = R->size1;					// Number of outputs
	size_t DIM = xZero->size + theta->size;	// Dimension of the system

	// Set the confidence level for chi2 test
	// double d = 1-p;

	// *** 1 *** Evaluate the "expected measurement" by solving the ODE system

	// Allocate the expected measurement matrix
	gsl_matrix * M2 = gsl_matrix_alloc (M->size1, M->size2);

	// Set up ODE solve environment
	gsl_odeiv_system sys = {ptrToFunc, ptrToJac, DIM, params};
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (PRECISION, PRECISION);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM);

	// Set initial time and initial step size
	double t = gsl_vector_get (time, 0);
	double k = 1e-10;

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


	// *** 2 *** Compute the variances

	// Compute the difference, the result is stored in M2
	gsl_matrix_sub (M2, M);

	for (i = 0; i < N; i++) {
		gsl_matrix_view out = gsl_matrix_submatrix (M2, i, 0, 1, S);
		gsl_vector_set (vhp, i, gsl_stats_variance ((&out.matrix)->data, 1, S));
	}

	// Free all the manually allocated resources
	gsl_matrix_free (M2);
	gsl_odeiv_step_free (s);
	gsl_odeiv_control_free (c);
	gsl_odeiv_evolve_free (e);

	return GSL_SUCCESS;
}


int chi2cost (gsl_vector * ave, gsl_vector * var,  const gsl_vector * time, const gsl_vector * theta, const gsl_vector * xZero, const gsl_matrix * M,
			  const gsl_matrix * R, int (* ptrToFunc) (double, const double[], double[], void *),
			  int (* ptrToJac) (double, const double[], double *, double[], void *), int (* ptrToOutput) (double, const double[], double[], void *),
			  void * params)
{
	// Detect the dimensions of the problem
	size_t S = time->size;					// Number of samples
	size_t N = R->size1;					// Number of outputs
	size_t DIM = xZero->size + theta->size;	// Dimension of the system

	// *** 1 *** Evaluate the "expected measurement" by solving the ODE system

	// Allocate the expected measurement matrix
	gsl_matrix * M2 = gsl_matrix_alloc (M->size1, M->size2);

	// Set up ODE solve environment
	gsl_odeiv_system sys = {ptrToFunc, ptrToJac, DIM, params};
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-8, 1e-8);
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

	// Compute the difference, the result is stored in M2
	gsl_matrix_sub (M2, M);

	for (i = 0; i < N; i++)
	{
		gsl_vector_const_view out = gsl_matrix_const_row (M, i);
		gsl_vector_const_view out2 = gsl_matrix_const_row (M2, i);
        double Z = gsl_stats_mean((&out2.vector)->data, (&out2.vector)->stride, (&out2.vector)->size);
        double Z1 = gsl_stats_mean((&out.vector)->data, (&out.vector)->stride, (&out.vector)->size);
        double Y = gsl_stats_variance ((&out2.vector)->data, (&out2.vector)->stride, (&out2.vector)->size);
        double Y1 = gsl_matrix_get (R, i, i);
        double Z2 = 2.58; //// **** NORMALIZATION value for mean, if noise is not zero mean

        gsl_vector_set (ave, i, (Z + Z2)/Z1);
		gsl_vector_set (var, i, (Y - Y1)/Y1);
	}

	// Free all the manually allocated resources
	gsl_matrix_free (M2);
	gsl_odeiv_step_free (s);
	gsl_odeiv_control_free (c);
	gsl_odeiv_evolve_free (e);

	return GSL_SUCCESS;
}
