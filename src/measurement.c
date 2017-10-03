/*
 *  measurement.c
 *
 *  This file contains the functions that implement the generation of the measurement
 *	signals.
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


/*	This function generates the measurement signal(s) at the required time steps by adding
 white Gaussian noise to the solution of the specified ODE system. */
int generate_measurements (int (* ptrToFunc) (double, const double[], double[], void *),
						   int (* ptrToJac) (double, const double[], double *, double[], void *),
						   int (* ptrToOutput) (double, const double[], double[], void *),
						   const gsl_matrix * time, gsl_matrix * sig, const gsl_matrix * R, const gsl_matrix * xZero, void * params)
{
	// ** 0 ** Detect problem dimensions
	size_t DIM = xZero->size1;
	size_t OUTDIM = R->size1;
	size_t N = time->size2;

	// ** 1 ** Solve the ODE system

	// Declare ODE system
	gsl_odeiv_system sys = {ptrToFunc, ptrToJac, DIM, params};

	// Declare step type
	const gsl_odeiv_step_type * T = gsl_odeiv_step_bsimp;
	gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, DIM);

	// Initialize the adaptive step size control objects
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-15, 1e-15);
	gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (DIM);

	// Set initial time and initial step size
	double t = gsl_matrix_get (time, 0, 0);
	double k = 1e-10;

	// Declare array y and initialize to initial conditions
	double y[DIM];
	int l;
	for (l = 0; l < DIM; l++)
	{
		y[l] = gsl_matrix_get (xZero, l, 0);
	}

	// Declare array for the "output" to store the results for each desired time step
	double h[OUTDIM];

	// Initialize counters
	int i, j;

	for (i = 0; i < N; i++)
	{
		// Loop to obtain a result for every specified time step
		while (t < gsl_matrix_get (time, 0, i))
		{
			// Evolve the system to the next time step
			int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t,  gsl_matrix_get (time, 0, i), &k, y);

			// Break loop if an error occurs
			if (status != GSL_SUCCESS)
		        break;
		}
		// Evaluate the output of the system at the current time step
		ptrToOutput (t, y, h, params);

		// Store the computed output in the array sig
		for (j = 0; j < OUTDIM; j++)
		{
			gsl_matrix_set (sig, j, i, h[j]);
		}

		// At this point every row of sig contains one output
	}

	// ** 2 ** Add gaussian noise to the output

	// Initialize the random number generator
	const gsl_rng_type * TT;
	gsl_rng * r;

	gsl_rng_env_setup();

	TT = gsl_rng_taus;
	r = gsl_rng_alloc (TT);

	// Read the state of the random number generator from file
	FILE * rng_state = fopen ("rng_state.dat", "rb");
	if (rng_state == NULL)
	{
		// File cannot be opened, or does not exist. Set the random number generator
		// state to the default value.
		gsl_rng_set (r, 0);
	}
	else
	{
		// Load the state of the generator
		if (gsl_rng_fread (rng_state, r) == GSL_EFAILED)
		{
			// State cannot be read. Exit.
			printf("Error reading random number generator state. Exiting...");
			return -1;
		}
	}
	fclose (rng_state);

	// Now we can finally get the random numbers from the generator
	for (j = 0; j < OUTDIM; j++)
	{
		for (i = 0; i < N; i++)
		{
			gsl_matrix_set (sig, j, i, gsl_matrix_get (sig, j, i) + gsl_ran_gaussian_ziggurat (r, sqrt (gsl_matrix_get (R, j, j))));
		}
	}

	// Write the new state of the random number generator to file
	rng_state = fopen ("rng_state.dat", "wb");
	if (rng_state == NULL)
	{
		// File cannot be opened. Exit.
		printf("Error opening random number generator state. Exiting...");
		return -1;
	}
	else
	{
		// Write the state of the generator
		if (gsl_rng_fwrite (rng_state, r) == GSL_EFAILED)
		{
			// State cannot be read. Exit.
			printf("Error reading random number generator state. Exiting...");
			return -1;
		}
	}
	fclose (rng_state);

	// Free the allocated resources
	gsl_rng_free (r);
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	// Signal that the computation was done successfully
	return GSL_SUCCESS;
}
