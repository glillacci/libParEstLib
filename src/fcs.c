/*
 *  fcs.c
 *
 *  This file contains functions to read data from FCS files.
 *
 *  Created by Gabriele Lillacci in April 2012.
 *	Latest revision: April 2012.
 *
 *
 *	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
 *	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
 *	is given to the original author, all derivative works are distributed under the same license or a compatible one,
 *	and this software and its derivatives are not used for commercial purposes.
 *	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
 *	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.
 */

#include <parestlib.h>

/**
 This function scans the TEXT segment of an FCS file, and positions the file cursor to
 the specified keyword in the FCS file. The return value is 0 if the seeking was performed
 successfully, or a non-zero value otherwise.
 */
int fcs_seek_to_keyword (FILE * data, const char * keyword)
{
	// Initialize reading buffer and return value
	char buf[256];
	int R = 1;

	// Go back to the starting position on the file
	rewind (data);

	// Keep filling up buffer with data from the file
	while (fgets (buf, sizeof (buf), data) != NULL)
	{
		// Check if the buffer contains the desired keyword
		char * off = strstr (buf, keyword);
		if (off != NULL)
		{
			// If it does, set the cursor to the beginning of the keyword
			R = fseek (data, off - buf - sizeof (buf) + 1, SEEK_CUR);
			break;
		}
	}

	// Return
	return R;
}


/**
 This function reads the value of an integer-valued keyword in an FCS version 3 file.
 The return value is the value of the keyword if the read was successful, or -1 otherwise.
 */
int fcs3_read_int_kw (FILE * data, const char * keyword)
{
	// Initialize the return value
	int R = -1;

	// Initialize a read buffer
	char buf[256];

	// Position the cursor at the desired keyword
	if (fcs_seek_to_keyword (data, keyword) == 0)
	{
		// Read the value
		if (fscanf (data, "%[^\|]\|%d\|", buf, &R) != 2)
			R = -1;
	}

	// Return
	return R;
}


/**
 This function reads the value of an integer-valued keyword in an FCS version 2 file.
 The return value is the value of the keyword if the read was successful, or -1 otherwise.
 */
int fcs2_read_int_kw (FILE * data, const char * keyword)
{
	// Initialize the return value
	int R = -1;

	// Initialize a read buffer
	char buf[256];

	// Position the cursor at the desired keyword
	if (fcs_seek_to_keyword (data, keyword) == 0)
	{
		// Read the value
		if (fscanf (data, "%[^\\]\\%d\\", buf, &R) != 2)
			R = -1;
	}

	// Return
	return R;
}


/**
 This function reads IEEE-754 single precision floating point numbers (float in C) from the
 specified file. The byte order is swapped to change the endianness of the data.
 */
int fread_floats_swap (float * f, size_t count, FILE * data)
{
	// Initialize the read buffer
	char bbuf[4];
	char bbuf2[4];

	// Initialize a float
	float g;

	// Loop for the number of floats to read
	for (size_t i = 0; i < count; i++)
	{
		// Fill up the buffer
		fread (bbuf, 1, 4, data);

		// Swap the bytes to change the endianness
		bbuf2[0] = bbuf[3];
		bbuf2[1] = bbuf[2];
		bbuf2[2] = bbuf[1];
		bbuf2[3] = bbuf[0];

		// Copy the byte array into the float
		memcpy (&g, bbuf2, 4);

		// Place the received float into the output array
		f[i] = g;
	}

	return (int) count;
}


/**
 This function reads IEEE-754 single precision floating point numbers (float in C) from the
 specified file. The byte order is not swapped, so endianness is preserved.
 */
int fread_floats_noswap (float * f, size_t count, FILE * data)
{
	// Initialize the read buffer
	char bbuf[4];

	// Initialize a float
	float g;

	// Loop for the number of floats to read
	for (size_t i = 0; i < count; i++)
	{
		// Fill up the buffer
		fread (bbuf, 1, 4, data);

		// Copy the byte array into the float
		memcpy (&g, bbuf, 4);

		// Place the received float into the output array
		f[i] = g;
	}

	return (int) count;
}
