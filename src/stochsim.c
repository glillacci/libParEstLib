/*
 *  stochsim.c
 *
 *  This file contains functions to implement simulation algorithms for stochastic
 *  chemical kinetics models.
 *
 *  Created by Gabriele Lillacci in July 2010.
 *	Latest revision: July 2010.
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
 This function evolves the stochatic model described by the structure model to the next
 reaction by using the SSA direct method.
 */
int ssa_direct_evolve (stochmod * model, const gsl_vector * params, gsl_vector * X, double * t, const gsl_rng * r)
{
	// Check the parameters
	if (gsl_vector_isnonneg (params) != 1) {
		printf(">>> ssa_direct_evolve: parameters vector contains unfeasible values\n");
		print_gsl_vector (params);
	}

	// Allocate necessary variables
	gsl_vector * prop = gsl_vector_alloc (model->nrxns);

	// Evaluate the reaction propensities at the current state
	model->propensity (X, params, prop);

	// Compute the sum a0
	double a0 = gsl_blas_dasum (prop);

	// Generate the two uniform random numbers
	double r1 = gsl_rng_uniform (r);
	double r2 = gsl_rng_uniform (r);

	//print_gsl_vector(X);
	//print_gsl_vector(prop);
	//printf("r1 = %f, r2 = %f, a0 = %f, r2*a0 = %f\n", r1, r2, a0, r2*a0);

	// Compute the time to the next reaction
	double tau = (1/a0)*log(1/r1);

	// Compute the index of the next reaction
	size_t mu = 0;
	double cumsum = 0; // cumsum has to be the sum from the first reaction to the current MINUS ONE!!!
	while (cumsum <= r2*a0)
	{
		cumsum += gsl_vector_get (prop, mu);
		mu++;
		//printf("mu = %d, cumsum = %f\n", (int) mu, cumsum);
	}
	//printf("Reaction that fired = %d\n", (int) mu);

	// Adjust the time
	*t = *t + tau;

	// Update the state according to the reaction that fired
	model->update (X, mu-1);

	// Free the manually allocated resources
	gsl_vector_free (prop);

	// Signal that computation was done correctly
	return GSL_SUCCESS;
}


/**
 This function allocates memory for a new rxn_sample_path object.
 */
rxn_sample_path * rxn_sample_path_alloc (size_t nspecies, size_t ntimes)
{
	// Allocate memory for the structure
	rxn_sample_path * rsp = (rxn_sample_path *) malloc (sizeof (rxn_sample_path));

	// Check if allocation was successful
	if (rsp != NULL)
	{
		// We can proceed and allocate the matrix and vector
		rsp->ntimes = ntimes;
		rsp->nspecies = nspecies;

		rsp->times = gsl_vector_alloc (ntimes);
		rsp->counts = gsl_matrix_alloc (nspecies, ntimes);

		return rsp;
	}
	else
	{
		// We are out of memory
		printf (">> error in function rxn_sample_path_alloc: could not allocate rxn_sample_path object...");
		return NULL;
	}
}


/**
 This function destroys a rxn_sample_path object.
 */
void rxn_sample_path_free (rxn_sample_path * rsp)
{
	// Destroy the vector and matrix
	gsl_vector_free (rsp->times);
	gsl_matrix_free (rsp->counts);

	// Free the struct
	free (rsp);
}


/**
 This function generates a discretely observed trajectory from the stochastic model described
 by the structure model using the SSA direct method.
 */

int ssa_direct_trajectory (stochmod * model, const gsl_vector * params, gsl_vector * times, gsl_vector * X0,
						   rxn_sample_path * path, const gsl_rng * r)
{
	// Intialization
	double t = gsl_vector_get (times, 0);
	double tf;

	// Transfer the times vector into the path struct
	gsl_vector_memcpy (path->times, times);

	// The vector X0 must be the first column of the matrix count
	gsl_vector_view col = gsl_matrix_column (path->counts, 0);
	gsl_vector_memcpy (&col.vector, X0);

	// Loop for the desired time points
	for (size_t i = 1; i < times->size; i++)
	{
		// Get next time point
		tf = gsl_vector_get (path->times, i);

		// Evolve the system until the time point is hit
		while (t < tf)
		{
			ssa_direct_evolve (model, params, X0, &t, r);
		}

		// At this point the state vector contains what we want
		gsl_vector_view col = gsl_matrix_column (path->counts, i);
		gsl_vector_memcpy (&col.vector, X0);
	}

	// Signal that computation was done correctly
	return GSL_SUCCESS;
}



/**
 This function allocates memory for a new rxn_ensemble object.
 */
rxn_ensemble * rxn_ensemble_alloc (size_t nreplic, size_t nspecies, size_t ntimes)
{
	// Allocate memory for the struct
	rxn_ensemble * ens = (rxn_ensemble *) malloc (sizeof (rxn_ensemble));

	// Check if allocation was successful
	if (ens != NULL)
	{
		// Set the struct fields
		ens->nreplic = nreplic;
		ens->nspecies = nspecies;
		ens->ntimes = ntimes;

		// Allocate the double array of pointers to rxn_sample path
		ens->data = (rxn_sample_path **) malloc (ens->nreplic * sizeof (rxn_sample_path *));

		// Proceed with the allocation
		for (size_t i = 0; i < ens->nreplic; i++)
		{
			ens->data[i] = rxn_sample_path_alloc (ens->nspecies, ens->ntimes);
		}

		// Return the pointer
		return ens;
	}
	else
	{
		// We are out of memory
		printf (">> error in function rxn_ensemble_alloc: could not allocate rxn_ensemble object...");
		return NULL;
	}
}


/**
 This function destroys a rxn_ensemble object.
 */
void rxn_ensemble_free (rxn_ensemble * ens)
{
	// Destroy all the single rxn_sample_path objects
	for (size_t i = 0; i < ens->nreplic; i++)
	{
		rxn_sample_path_free (ens->data[i]);
	}

	// Destroy the data structure
	free (ens->data);

	// Destroy the struct
	free (ens);
}


/**
 This function extracts the counts of desired outputs from a rxn_ensemble data structure. The counts
 are returned in the form of a gsl_matrix object, with the following format:
			---time 1---	---time 2---	...		---time K---
			o1 o2 ... oP	o1 o2 ... oK			o1 o2 ... oK
 sim 1		x  x  ... x
 sim 2      x  x  ... x
 ...        ...........
 sim N      x  x  ... x
 */
int rxn_ensemble_counts (const rxn_ensemble * ens, gsl_matrix * counts, const gsl_matrix * out)
{
	// Allocate necessary objects
	gsl_matrix * outs = gsl_matrix_alloc (out->size1, ens->ntimes);
	gsl_matrix * colt = gsl_matrix_alloc (1, outs->size1);

	// Declare necessary objects
	gsl_matrix_view col;
	gsl_matrix_view sro;

	for (size_t i = 0; i < ens->nreplic; i++)
	{
		// Multiply the counts matrix by the output matrix to generate the outputs
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, out, ens->data[i]->counts, 0.0, outs);

		// Scan the columns of the matrix outs
		for (size_t j = 0; j < outs->size2; j++) {
			// Create a submatrix view of the outs matrix
			col = gsl_matrix_submatrix (outs, 0, j, outs->size1, 1);
			gsl_matrix_transpose_memcpy(colt, &col.matrix);

			// Create a submatrix view of the counts matrix
			sro = gsl_matrix_submatrix (counts, i, outs->size1*j, 1, outs->size1);

			// Copy the content into the appropriate location of counts
			gsl_matrix_memcpy (&sro.matrix, colt);
		}
	}

	// Free the manually allocated resources
	gsl_matrix_free (outs);
	gsl_matrix_free (colt);

	// Signal that computation was done successfully
	return GSL_SUCCESS;
}


/**
 This function generates a discretely observed trajectory ensemble from the stochastic model
 described by the structure model using the SSA direct method.
 */
int ssa_direct_ensemble (stochmod * model, const gsl_vector * params, gsl_vector * times, const gsl_vector * X0,
						 rxn_ensemble * ens, const gsl_rng * r)
{
	// Create a copy of the initial conditions vector
	gsl_vector * X0cpy = gsl_vector_alloc (X0->size);

	// Loop for number of replicates
	for (size_t i = 0; i < ens->nreplic; i++)
	{
		// Reset the initial conditions vector to the correct quantities
		gsl_vector_memcpy (X0cpy, X0);

		// Generate a sample path
		ssa_direct_trajectory (model, params, times, X0cpy, ens->data[i], r);
	}

	// Free the manually allocated resources
	gsl_vector_free (X0cpy);

	// Signal that computation was done correctly
	return GSL_SUCCESS;
}
