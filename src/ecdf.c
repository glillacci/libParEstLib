/*
 *  ecdf.c
 *
 *  This file contains the functions that allow to compute and manipulate
 *	the empirical cumulative distribution functions objects (ecdf).
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
 This function allocates a new ecdf struct and fills it up with the empirical cumulative
 distribution function of the data points in the vector data.
 The resulting object must be manually destroyed with the function ecdf_free.
 */
ecdf * ecdf_create_from_data (const gsl_vector * data)
{
	// Create a temporary copy of the data vector and sort it
	gsl_vector * tmpdata = gsl_vector_alloc (data->size);
	gsl_vector_memcpy (tmpdata, data);
	gsl_sort_vector (tmpdata);

	// Get the number of data points
	size_t npoints = tmpdata->size;

	// Allocate two empty vectors and set them to all zeros
	gsl_vector * tmpdata2 = gsl_vector_calloc (npoints);
	gsl_vector * freq = gsl_vector_calloc (npoints);

	// Copy the first element of tmpdata to tmpdata2
	gsl_vector_set (tmpdata2, 0, gsl_vector_get (tmpdata, 0));

	// Initialize indices
	size_t i = 1; // Index to scan the tmpdata vector
	size_t j = 0; // Index for tmpdata2 vector
	size_t f = 1; // Frequency count

	// Scan tmpdata for multiple occurrencies of the same value
	while (i < tmpdata->size)
	{
		if (gsl_vector_get (tmpdata2, j) == gsl_vector_get (tmpdata, i))
		{
			// The next element of tmpdata is equal to the last element of tmpdata2
			i++;
			f++;
		}
		else
		{
			// The next element of tmpdata is different from the last element of tmpdata2
			// Set the frequency of the last element
			gsl_vector_set (freq, j, f);
			// Copy the next element
			j++;
			gsl_vector_set (tmpdata2, j, gsl_vector_get (tmpdata, i));
			// Increase i
			i++;
			// Reset f
			f = 1;
		}
	}

	// Process the last data point
	gsl_vector_set (tmpdata2, j, gsl_vector_get (tmpdata, i-1));
	gsl_vector_set (freq, j, f);

	// The index j is now set to the last element that was transfered
	gsl_vector_view tmpdata2view = gsl_vector_subvector (tmpdata2, 0, j+1);
	gsl_vector_view	freqview = gsl_vector_subvector (freq, 0, j+1);

	// Allocate memory for the new ecdf object.
	ecdf * dist = (ecdf *) malloc (sizeof (ecdf));

	// If malloc fails, warn the user and return NULL
	if (dist == NULL) {
		printf(">>> warning: can not allocate memory for ecdf object...");
		return NULL;
	}
	// Otherwise...
	else
	{
		// Set the number of points in the ecdf
		dist->N = (&tmpdata2view.vector)->size;

		// Allocate the xdata and ydata vectors
		dist->xdata = gsl_vector_calloc (dist->N);
		dist->ydata = gsl_vector_calloc (dist->N);

		// Copy the sorted and trimmed data vector into the ecdf xdata
		gsl_vector_memcpy (dist->xdata, &tmpdata2view.vector);

		// Fill up the ydata vector
		// First point
		gsl_vector_set (dist->ydata, 0, gsl_vector_get(&freqview.vector, 0)/((double) npoints));
		// The other points
		for (size_t i = 1; i < dist->N; i++) {
			gsl_vector_set (dist->ydata, i, gsl_vector_get (dist->ydata, i-1) + gsl_vector_get(&freqview.vector, i)/((double) npoints));
		}

		// Free all the manually allocated resources
		gsl_vector_free (freq);
		gsl_vector_free (tmpdata2);
		gsl_vector_free (tmpdata);

		// Return the pointer
		return dist;
	}
}


/**
 This function evaluates the ECDF contained in dist at the point x.
 */
double ecdf_eval (const ecdf * dist, double x)
{
	double ans = 0.0;

	if (x < gsl_vector_get (dist->xdata, 0))
	{
		// Evaluation point is smaller than first data point
		ans = 0.0;
	}
	else if (x >= gsl_vector_get (dist->xdata, (dist->N)-1))
	{
		// Evaluation point is greater or equal than last data point
		ans = 1.0;
	}
	else {
		// Find the index of the first value of xdata not smaller than x
		size_t i = 0;
		while (gsl_vector_get (dist->xdata, i) <= x) {
			i++;
		}

		// Return the corresponding value of ydata
		ans = gsl_vector_get (dist->ydata, i-1);
	}

	// Return
	return ans;
}


/**
 This function generates one random sample from the ecdf dist.
 */
double ecdf_sample (const ecdf * dist, const gsl_rng * r)
{
	// Generate a standard uniform random number
	double u = gsl_rng_uniform (r);

	// Find the first value of f not smaller than u
	size_t i = 0;
	while (gsl_vector_get(dist->ydata, i) < u) {
		i++;
	}

	// Return the corresponding value of x
	return gsl_vector_get (dist->xdata, i);
}




/**
 This function frees an ecdf object.
 */
void ecdf_free (ecdf * dist)
{
	// Destroy the xdata and ydata vector
	gsl_vector_free (dist->xdata);
	gsl_vector_free (dist->ydata);

	// Free the ecdf struct
	free (dist);
}
