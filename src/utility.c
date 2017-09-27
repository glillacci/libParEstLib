/*
 *  utility.c
 *
 *  This file implements some functions that do silly things in an easy way!!
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

#include "parestlib.h"


/**
 * Matrix product for dummies!!
 */
int easy_matrix_product (gsl_matrix * p, const gsl_matrix * a, const gsl_matrix * b)
{
	// Compute the product using BLAS library
	int c;
	return gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, a, b, 0.0, p);
}


// Matrix inverse for dummies!!
/* REMARK
 Matrix inverse is an operation to be avoided whenever possible, but
 having it handy might be useful anyway. */
int easy_matrix_inverse (gsl_matrix * inv, const gsl_matrix * a)
{
	// Allocate the permutation p
	gsl_permutation * p = gsl_permutation_alloc (a->size1);
	
	// Create another copy of a to prevent its destruction by the LU decomposition
	gsl_matrix * aa = gsl_matrix_alloc(a->size1, a->size2);
	gsl_matrix_memcpy (aa, a);
	
	// Get the LU decomposition of aa
	int s;
	gsl_linalg_LU_decomp (aa, p, &s);
	
	// Get the inverse
	int r = gsl_linalg_LU_invert (aa, p, inv);
	
	// Free the allocated resources
	gsl_permutation_free (p);
	gsl_matrix_free (aa);
	
	// Return
	return r;
}


// MATLAB semicolon operator for making equally spaced arrays
int matlab_semicolon (gsl_matrix * array, double start, double step)
{
	int i;
	size_t samples = array->size2;
	
	for (i = 0; i < samples; i++)
	{
		gsl_matrix_set (array, 0, i, start+i*step);
	}
	
	return GSL_SUCCESS;
}


// This function concatenates array1 and array2 into array
int concatenate_arrays (double array[], const double array1[], int L1, const double array2[], int L2)
{
	int i;
	for (i=0; i < L1; i++)
	{
		array[i] = array1[i];
	}
	
	for (i=0; i < L2; i++)
	{
		array[i+L1] = array2[i];
	}
	
	return GSL_SUCCESS;
}


// This function prints a matrix to console
int print_matrix (double mat[], int rows, int cols)
{
	printf ("[%dx%d DOUBLE ARRAY]\n", rows, cols);

	int i, j;
	for (i=0; i < rows; i++)
	{
		for (j=0; j < cols; j++)
		{
			printf ("%12.6f \t ", mat[cols*i+j]);
		}
		printf ("\n");
	}
	printf ("\n");
	
	return i*j;
}


// This function prints a gsl_matrix to console
int print_gsl_matrix (const gsl_matrix * mat)
{
	printf ("\n[%dx%d GSL_MATRIX]\n", (int) mat->size1, (int) mat->size2);
	
	int i, j;
	for (i=0; i < mat->size1; i++)
	{
		for (j=0; j < mat->size2; j++)
		{
			printf (" %12.6e \t ", gsl_matrix_get (mat, i, j));
		}
		printf ("\n");
	}
	printf("\n");
	
	return i*j;
}


// This function prints a gsl_vector to console
int print_gsl_vector (const gsl_vector * v)
{
	printf ("\n[%d-GSL_VECTOR]\n", (int) v->size);
	
	int i;
	for (i=0; i < v->size; i++)
	{
		printf (" %12.6e \n", gsl_vector_get (v, i));
	}
	printf("\n");
	
	return i;
}

// This function prints a part of a gsl_vector object to console
int print_gsl_subvector (const gsl_vector * v, size_t start, size_t nel)
{
	printf ("\n[%d-GSL_VECTOR]\n", (int) v->size);
	
	if (start > 0)
		printf (" ...\n");
	
	for (size_t i = start; i < start + nel; i++)
	{
		printf (" %d. %12.6f \n", (int) i, gsl_vector_get (v, i));
	}
	
	if (start + nel < v->size)
		printf (" ...\n");
	
	printf("\n");
	
	return nel;
}


// This function returns the row number converted from linear indexing
int row_number (int linearIndex, int rows, int cols)
{
	return ceil (linearIndex/cols);
}

// This function returns the column number converted from linear indexing
int col_number (int linearIndex, int rows, int cols)
{
	return linearIndex - (row_number (linearIndex,rows,cols) *cols);
}

// This function computes the average of the elements of a gsl_matrix
double mean (gsl_matrix * mat)
{
	double sum = 0;
	int i, j, count = 0;
	
	for (i = 0; i < mat->size1; i++)
	{
		for (j = 0; j < mat->size2; j++)
		{
			sum += gsl_matrix_get (mat, i, j);
			count++;
		}
	}
	
	return sum/count;
}


// This function concatenates two gsl_matrix objects **horizontally**
int gsl_matrix_horzcat (gsl_matrix * cat, const gsl_matrix * mat1, const gsl_matrix * mat2)
{
	// Check that the dimensions are compatible
	if (mat1->size1 != mat2->size1)
	{
		printf ("\n>> error in function gsl_matrix_horzcat: ");
		printf ("the two matrices must have the same number of rows...\n");
		return GSL_EFAILED;
	}
	
	// Perform the concatenation
	int i, j;
	for (i = 0; i < mat1->size1; i++)
	{
		for (j = 0; j < mat1->size2; j++)
		{
			gsl_matrix_set (cat, i, j, gsl_matrix_get (mat1, i, j));
		}
	}
	
	for (i = 0; i < mat2->size1; i++)
	{
		for (j = 0; j < mat2->size2; j++)
		{
			gsl_matrix_set (cat, i, j+mat1->size2, gsl_matrix_get (mat2, i, j));
		}
	}
	
	return GSL_SUCCESS;
}


// This function concatenates two gsl_matrix objects **vertically**
int gsl_matrix_vertcat (gsl_matrix * cat, const gsl_matrix * mat1, const gsl_matrix * mat2)
{
	// Check that the dimensions are compatible
	if (mat1->size2 != mat2->size2)
	{
		printf ("\n>> error in function gsl_matrix_vertcat: ");
		printf ("the two matrices must have the same number of columns...\n");
		return GSL_EFAILED;
	}
	
	// Perform the concatenation
	int i, j;
	for (i = 0; i < mat1->size1; i++)
	{
		for (j = 0; j < mat1->size2; j++)
		{
			gsl_matrix_set (cat, i, j, gsl_matrix_get (mat1, i, j));
		}
	}
	
	for (i = 0; i < mat2->size1; i++)
	{
		for (j = 0; j < mat2->size2; j++)
		{
			gsl_matrix_set (cat, i+mat1->size1, j, gsl_matrix_get (mat2, i, j));
		}
	}
	
	return GSL_SUCCESS;
}


/** 
 This function sets the elements of a gsl_vector object to a specified double array.
 */
void gsl_vector_set_array (gsl_vector * v, double elem[])
{
	// Get a view for the array elem[]
	gsl_vector_view ev = gsl_vector_view_array (elem, v->size);
	
	// Set the elements of v
	gsl_vector_memcpy (v, &ev.vector);
}


/** 
 This function fills up a gsl_matrix with the elements of a constant double array in row major order.
 */
void gsl_matrix_set_array (gsl_matrix * mat, double elem[])
{
	// Get a view for the array elem[]
	gsl_matrix_view ev = gsl_matrix_view_array (elem, mat->size1, mat->size2);
	
	// Set the elements of mat
	gsl_matrix_memcpy (mat, &ev.matrix);
}


/**
 This function changes the size of a gsl_vector object. If the new size is larger, trailing zeros are added.
 */
gsl_vector * gsl_vector_realloc (gsl_vector * oldvec, size_t newsize)
{
	// Create the new vector
	gsl_vector * newvec = gsl_vector_calloc (newsize);

	// Transfer the elements
	for (size_t i = 0; i < ((newsize <= oldvec->size) ? newsize : oldvec->size); i++)
	{
		gsl_vector_set (newvec, i, gsl_vector_get (oldvec, i));
	}

	// Destroy the old vector
	gsl_vector_free (oldvec);

	return newvec;
}
