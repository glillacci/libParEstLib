/*
 *  deriv.c
 *
 *  This file contains the functions that implement various types of numerical
 *	derivative computations.
 *
 *  Created by Gabriele Lillacci in May 2009.
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


/**
 This function computes the numerical i-th partial derivative of a function func: Rn --> R. The derivative is
 computed at the point x using a central difference three-point formula with step size h. An estimate of the error is
 computed by subtracting the three-point derivative value and the corresponding five-point formula value.
*/
int central_diff_sv (double * D, double * E, double (* func) (const gsl_vector *, void *), size_t i, const gsl_vector * x, double h, void * params)
{
	// Get the current point for derivative computation
	double p = gsl_vector_get (x, i);

	// Make sure that h is machine representable
	double temp = h + p;
	h = temp - p;

	// Allocate necessary vector objects
	gsl_vector * Xp = gsl_vector_alloc (x->size);

	// Add h to the i-th component of Xp
	gsl_vector_memcpy (Xp, x);
	gsl_vector_set (Xp, i, gsl_vector_get (Xp, i) + h);
	// Evaluate F at x + h*ei
	double Fp = func (Xp, params);

	// Add another h to the i-th component of Xp (so that it's now the same as x + 2*h*ei)
	gsl_vector_set (Xp, i, gsl_vector_get (Xp, i) + h);
	// Evaluate F at x + 2*h*ei
	double Fpp = func (Xp, params);

	// Subtract 3*h to the i-th component of Xp (so that it's now the same as x - h*ei)
	gsl_vector_set (Xp, i, gsl_vector_get (Xp, i) - 3*h);
	// Evaluate F at x - h*ei
	double Fm = func (Xp, params);

	// Subtract another h to the i-th component of Xp (so that it's now the same as x - 2*h*ei)
	gsl_vector_set (Xp, i, gsl_vector_get (Xp, i) - h);
	// Evaluate F at x - 2*h*ei
	double Fmm = func (Xp, params);

	// Compute D using the three-point central difference
	double D3 = (Fp - Fm)/(2*h);
	// Compute D using the five-point central difference
	double D5 = (- Fpp + 8*Fp - 8*Fm + Fmm)/(12*h);

	// Set the value of D
	* D = D3;
	// Set the value of E
	* E = fabs ((D5-D3)/D3);

	// Free all the manually allocated resources
	gsl_vector_free (Xp);

	// Signal that the computation was done successfully
	return GSL_SUCCESS;
}


/**
 This function computes the numerical gradient of a function func: Rn --> R. The gradient is
 evaluated at the point x using a central difference three-point formula. An estimate of the error is
 computed by subtracting the three-point derivative value and the corresponding five-point formula value.
 The step size is adaptively selected so that the relative accuracy in smaller than tol.
 */
int adaptive_gradient_sv (gsl_vector * df, gsl_vector * err, const gsl_vector * x, double (* func) (const gsl_vector *, void *), double tol, void * params)
{
	// Initialize the variables D and E
	double D=0, E=0;

	// Initialize error count
	size_t errc=0;

	// Cycle for function component
	size_t i;
	for (i = 0; i < df->size; i++)
	{
		// Select initial h = sqrt(eps) * |x_i|
		double h = sqrt (2.2204e-16) * fabs (gsl_vector_get (x, i));

		// Compute the derivative & relative error
		central_diff_sv (&D, &E, func, i, x, h, params);

		// Check if the error tolerance is satisfied
		if (E <= tol)
		{
			// Error is lower that tol --> set the appropriate component of the gradient and we are done
			gsl_vector_set (df, i, D);
			gsl_vector_set (err, i, E);
		}
		else
		{
			// The derivative is not accurate enough --> use the adaptive approach
			size_t count = 0;

			// Choose an initially big step size (not smaller than 1e-3)
			h = (0.1 * fabs (gsl_vector_get (x, i)) < 1e-3) ? 1e-3 : 0.1 * fabs (gsl_vector_get (x, i));

			// Compute the derivative and error
			central_diff_sv (&D, &E, func, i, x, h, params);

			// A maximum of 50 iterations is allowed
			while (E >= tol && count <= 50)
			{
				// Save the last value of the derivative and of the error
				double D2 = D;
				double E2 = E;

				// Compute the new derivative and error
				central_diff_sv (&D, &E, func, i, x, h, params);

				// Compare the current error and the previous error
				if (E > E2 && count > 3)
				{
					// If the current error if bigger than the previous error (and already 3 iterations were completed), the
					// optimal step-size h has been overstepped, so we need to go back the previous derivative and exit.
					D = D2;
					E = E2;
					break;
				}

				// Set the new h, by making sure the reduction is not too big
				double r = (cbrt (0.9 * tol/E) < 0.2) ? 0.2 : cbrt (0.9 * tol/E);
				h = h * r;

				// Move to next iteration
				count++;
			}

			// The while loop ended, so at this point we either have a good derivative or we have too many iterations!!
			if (E <= tol)
			{
				// The derivative is good --> set it and exit
				gsl_vector_set (df, i, D);
				gsl_vector_set (err, i, E);
			}
			else
			{
				// The derivative is not good --> warn the user and continue
				printf (">> warning... cannot achieve desired accuracy in gradient component %d, error is %g \n", (int) i, E);
				gsl_vector_set (df, i, D);
				gsl_vector_set (err, i, E);
				errc++;
			}
		}
	}

	// Check how many errors we have and return accordingly
	if (errc==0)
		return GSL_SUCCESS;
	else
		return GSL_EFAILED;
}
