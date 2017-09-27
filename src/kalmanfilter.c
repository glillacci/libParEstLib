/*
 *  kalmanfilter.c
 *
 *  This file contains the functions that implement the single steps of the
 *	Hybrid extended Kalman filter.
 *
 *  Created by Gabriele Lillacci in July 2008.
 *	Latest revision: May 2009.
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


/*
 Internal functions == NOT TO BE USED DIRECTLY
 */


/*
 Problem definition for the solution of the Lyapunov equation. This is an auxiliary function
 for hekf_time_update_p and is not to be used directly.
 */
int funcL (double t, const double y[], double dydt[], void * eqn)
{
	// Recover problem dimension from Lyapunov Equation structure
	size_t DIM = ((LyapunovEquation *) eqn)->DIM;

	// Define view for matrix P
	gsl_matrix_const_view P_view = gsl_matrix_const_view_array (y, DIM, DIM);

	// Symmetrize P by setting P = (P + P^T)/2
	gsl_matrix * P = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix_transpose_memcpy (P, &P_view.matrix);
	gsl_matrix_add (P, &P_view.matrix);
	gsl_matrix_scale (P, 0.5);

	// Get the A and Q matrices from the eqn struct
	gsl_matrix * A = ((LyapunovEquation *) eqn)->A;
	gsl_matrix * Q = ((LyapunovEquation *) eqn)->Q;

	// Declare matrices for A*P and P*A^T
	gsl_matrix * AP = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix * PAT = gsl_matrix_alloc (DIM, DIM);

	// Compute A*P (A times P)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, A, P, 0.0, AP);

	// Compute P*A^T (P times transpose of A)
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, P, A, 0.0, PAT);

	// Add elements of A*P to P*A^T. The sum is stored in the matrix PAT
	gsl_matrix_add (PAT, AP);

	// Add elements of Q to matrix PAT. The sum is stored in the matrix PAT
	gsl_matrix_add (PAT, Q);

	// Set the values of the array dydt using gsl_matrix_view
	gsl_matrix_view dydt_mat = gsl_matrix_view_array (dydt, DIM, DIM);
	gsl_matrix_memcpy (&dydt_mat.matrix, PAT);

	// Free allocated matrices
	gsl_matrix_free (P);
	gsl_matrix_free (PAT);
	gsl_matrix_free (AP);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/*
 Jacobian definition for the solution of the Lyapunov equation. This is an auxiliary function
 for hekf_time_update_p and is not to be used directly.
 */
int jacL (double t, const double y[], double * dfdy, double dfdt[],  void * eqn)
{
	// Recover problem dimension from Lyapunov Equation structure
	size_t DIM = ((LyapunovEquation *) eqn)->DIM;

	// Get the matrix A from the eqn struct
	gsl_matrix * A = ((LyapunovEquation *) eqn)->A;

	// Create A transpose
	gsl_matrix * AT = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix_transpose_memcpy (AT, A);

	// Add the matrices A and AT. The sum is stored in AT
	gsl_matrix_add (AT, A);

	// Set the values of the array dydt[] using gsl_matrix_view
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, DIM, DIM);
	gsl_matrix_memcpy (&dfdy_mat.matrix, AT);

	// Set dfdt to 0
	int i;
	for (i = 0; i < DIM*DIM; i++)
	{
		dfdt[i] = 0.0;
	}

	// Free allocated matrices
	gsl_matrix_free(AT);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/*
 Exported functions
 */


/**
 This function projects the state ahead from the current time step tCurrent
 to the next time step tNext.
 */
int hekf_time_update_x (double tCurrent, double tNext, const gsl_vector * xPlus, gsl_vector * xMinus,
						int (* ptrToFunc) (double, const double[], double[], void *),
						int (* ptrToJac) (double, const double[], double *, double[], void *), void * params)
{
	// Detect problem dimension
	size_t DIM = xPlus->size;

	// Set up ODE solver environment
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (PRECISION, PRECISION);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM);
	gsl_odeiv_system sys = {ptrToFunc, ptrToJac, DIM, params};

	// Declare variables for time and step size
	double t = tCurrent, h = INITSTEP;

	// Declare array y and set the initial conditions
	double y[DIM];
	int i;
	for (i = 0; i < DIM; i++)
	{
		y[i] = gsl_vector_get (xPlus, i);
	}

	// Get the starting time
	clock_t tic = clock();

	// Main solution loop
	while (t < tNext)
	{
		// Monitor the elapsed time since the function has started solving the system
		clock_t toc = clock();
		// If more than MAXTIME seconds have passed, abort the computation
		if ((((double) (toc - tic)) / CLOCKS_PER_SEC) > MAXTIME)
		{
			printf ("\n>> error in function hekf_time_udpate_x: ");
			printf ("the solution of the ODE system is taking too long, function exiting...\n");
			return GSL_EFAILED;
		}

		// Evolve the system to the next time step
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tNext, &h, y);

		// Break loop if an error occurs
		if (status != GSL_SUCCESS)
		{
			printf ("\n>> error in function hekf_time_udpate_x: ");
			printf ("error in the GSL library solver, function exiting...\n");
			return GSL_EFAILED;
		}
	}

	// Store the values of y into xMinus
	gsl_vector_set_array(xMinus, y);

	// Free all the manually allocated resources
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/*	This function projects the state ahead from the current time step tCurrent
 to the next time step tNext, using the second order approximation.
int hekf_time_update_x_second_order (double tCurrent, double tNext, const double xPlus[], double xMinus[],
									 const double pPlus[], int (* ptrToFunc2o) (double, const double[], double[], void *),
									 int (* ptrToJac) (double, const double[], double *, double[], void *), void * params)
{
	// Set up the params2 struct for the ODE systems
	struct {
		parameters * p;
		const double * pP;
	} params2;

	params2.p = params;
	params2.pP = pPlus;

	// Set up ODE solver environment
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (PRECISION, PRECISION);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM);
	gsl_odeiv_system sys = {ptrToFunc2o, ptrToJac, DIM, &params2};

	// Declare variables for time and step size
	double t = tCurrent, h = INITSTEP;

	// Declare array y and set the initial conditions
	double y[DIM];
	int i;
	for (i = 0; i < DIM; i++)
	{
		y[i] = xPlus[i];
	}

	// Get the starting time
	clock_t tic = clock();

	// Main solution loop
	while (t < tNext)
	{
		// Monitor the elapsed time since the function has started solving the system
		clock_t toc = clock();
		// If more than 30 seconds have passed, abort the computation
		if ((((double) (toc - tic)) / CLOCKS_PER_SEC) > MAXTIME)
		{
			printf ("\n>> error in function hekf_time_udpate_x: ");
			printf ("the solution of the ODE system is taking too long, function exiting...\n");
			return GSL_EFAILED;
		}

		// Evolve the system to the next time step
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tNext, &h, y);

		// Break loop if an error occurs
		if (status != GSL_SUCCESS)
		{
			printf ("\n>> error in function hekf_time_udpate_x: ");
			printf ("error in the GSL library solver, function exiting...\n");
			return GSL_EFAILED;
		}
	}

	// Store the values of y into xMinus
	for (i = 0; i < DIM; i++)
	{
		xMinus[i] = y[i];
	}

	// Free all the manually allocated resources
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
} */


/**
 This function projects the error covariance ahead from the current time step tCurrent
 to the next time step tNext.
 */
int hekf_time_update_p (double tCurrent, double tNext, const gsl_matrix * pPlus, gsl_matrix * pMinus, LyapunovEquation * eqn)
{
	// Detect problem dimension
	size_t DIM = ((LyapunovEquation *) eqn)->DIM;

	// Set up ODE solver environment
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM*DIM);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (PRECISION, PRECISION);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM*DIM);
	gsl_odeiv_system sys = {funcL, jacL, DIM*DIM, eqn};

	// Declare variables for time and initial step size
	double t = tCurrent, h = INITSTEP;

	// Declare array y and set the initial conditions
	double y[DIM*DIM];
	int i;
	for (i = 0; i <  DIM*DIM; i++)
	{
        y[i] = pPlus->data[i];
	}

	// Get the starting time
	clock_t tic = clock();

	// Main solution loop
	while (t < tNext)
	{
		// Monitor the elapsed time since the function has started solving the Lyapunov equation
		clock_t toc = clock();
		// If more than MAXTIME seconds have passed, abort the computation
		if ((((double) (toc - tic)) / CLOCKS_PER_SEC) > MAXTIME)
		{
			printf ("\n>> error in function hekf_time_udpate_p: ");
			printf ("the solution of the lyapunov equation is taking too long, function exiting...\n");
			return GSL_EFAILED;
		}

		// Evolve the system to the next time step
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tNext, &h, y);

	    // Break loop if an error occurs
		if (status != GSL_SUCCESS)
		{
			printf ("\n>> error in function hekf_time_udpate_x: ");
			printf ("error in the GSL library solver, function exiting...\n");
			return GSL_EFAILED;
		}
	}

	// Assign the values of y[] to pMinus
	gsl_matrix_set_array (pMinus, y);

	// Free the allocated resources
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/**
 This function computes the optimal filter gain.
 */
int hekf_compute_gain (gsl_matrix * L, const gsl_matrix * pMinus, const gsl_matrix * H, const gsl_matrix * R)
{
	// Detect problem dimensions
	size_t DIM = pMinus->size1;
	size_t OUTDIM = H->size1;

	// Allocate necessary matrices
	gsl_matrix * W = gsl_matrix_alloc (OUTDIM, OUTDIM);
	gsl_matrix * T1 = gsl_matrix_alloc (OUTDIM, DIM);
	gsl_matrix * T2 = gsl_matrix_alloc (DIM, OUTDIM);

	// Compute H * P (result stored in T1)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, H, pMinus, 0.0, T1);
	// Compute H * P * H^T (i.e. T1 * H^T, result stored in W)
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, T1, H, 0.0, W);

	// Compute (H * Pminus * H^T + R) (result stored in W)
	gsl_matrix_add (W, R);
	// Invert W
	easy_matrix_inverse (W, W);

	// Compute H^T * W^(-1) (result stored in T2)
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, H, W, 0.0, T2);
	// Compute L = P * H^T * W^(-1) = P * T2 (result stored in L)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, pMinus, T2, 0.0, L);

	// Free all the manually allocated resources
	gsl_matrix_free (W);
	gsl_matrix_free (T1);
	gsl_matrix_free (T2);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/**
 This function corrects the predicted state estimate with the information from the output.
 */
int hekf_measurement_update_x (double tCurrent, const gsl_vector * xMinus, gsl_vector * xPlus, const gsl_matrix * L,
							   const gsl_vector * mCurrent, int (* ptrToOutput) (double, const double[], double[], void *), void * params)
{
	// Detect problem dimensions
	size_t OUTDIM = L->size2;

	// Allocate necessary gsl_matrix objects
	gsl_vector * h = gsl_vector_alloc (OUTDIM);

	// Evaluate the output corresponding to xMinus (expected measure)
	ptrToOutput (tCurrent, xMinus->data, h->data, params);

	// Compute the difference h(x)-M (result stored in h, att. the signs are opposite!!)
	gsl_vector_sub (h, mCurrent);

	// Multiply h(x)-y by L and change sign (result stored in X, now the signs match again!!)
	gsl_blas_dgemv (CblasNoTrans, -1.0, L, h, 0.0, xPlus);

	// Finally, add xMinus and the computation of xPlus is complete
	gsl_vector_add (xPlus, xMinus);

	// Free the manually allocated resources
	gsl_vector_free (h);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}


/**
 This function corrects the predicted error covariance with the information from the output.
 */
int hekf_measurement_update_p (gsl_matrix * pPlus, const gsl_matrix * pMinus, const gsl_matrix * L,
							   const gsl_matrix * H, const gsl_matrix * R)
{
	// Detect problem dimensions
	size_t DIM = L->size1;
	size_t OUTDIM = L->size2;

	// Allocate the necessary gsl_matrix objects
	gsl_matrix * I = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix_set_identity (I);
	gsl_matrix * T1 = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix * T2 = gsl_matrix_alloc (DIM, DIM);
	gsl_matrix * T3 = gsl_matrix_alloc (OUTDIM, DIM);
	gsl_matrix * T4 = gsl_matrix_alloc (DIM, DIM);

	// Compute -L*H (result stored in T1)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, L, H, 0.0, T1);
	// Compute I-L*H = -L*H+I = T1 + I (result stored in T1)
	gsl_matrix_add (T1, I);
	// Compute pMinus * (I-L*H)^T = pMinus * T1^T (result stored in T2)
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, pMinus, T1, 0.0, T2);
	// Compute (I-L*H) * pMinus * (I-L*H)^T = T1 * T2 (partial result stored in pPlus)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, T1, T2, 0.0, pPlus);

	// Compute R * L^T (result stored in T3)
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, R, L, 0.0, T3);
	// Compute L * R * L^T = L * T3 (result stored in T4)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, T3, 0.0, T4);
	// Add T4 to pPlus to get the final value!!
	gsl_matrix_add (pPlus, T4);

	// Free the manually allocated resources
	gsl_matrix_free (I);
	gsl_matrix_free (T1);
	gsl_matrix_free (T2);
	gsl_matrix_free (T3);
	gsl_matrix_free (T4);

	// Signal that the computation was done correctly
	return GSL_SUCCESS;
}
