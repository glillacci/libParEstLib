/*
 *  minim.c
 *
 *  This file contains the functions that implement various types of minimization,
 *	including constrained minimization and quadratic programming.
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
