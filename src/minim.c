/*
 *  minim.c
 *
 *  This file contains the functions that implement various types of minimization,
 *	including constrained minimization and quadratic programming.
 *
 *  Created by Gabriele Lillacci in June 2009.
 *	Latest revision: June 2009.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA. 
 */

#include "parestlib.h"


/**
 This function solves a convex quadratic program subject to equality contraints,
 of the following form:
 
					min F(x)  = 0.5 x^T*G*x + x^T*c
			 subject to Aeq*x = beq
 
 where G is positive definite and A is full rank. This function uses a nullspace method,
 so G is allowed to be singular.
 */
int quadprog_equality_nullspace (gsl_vector * sol, const gsl_matrix * G, const gsl_vector * c, const gsl_matrix * Aeq, const gsl_vector * beq)
{
	// Allocate necessary objects
	// gsl_matrix * U = gsl_matrix_alloc 
	// Find an orthonormal basis for the nullspace of Aeq
	
	
	return GSL_SUCCESS;
}