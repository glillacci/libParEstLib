#	Autoconf script for
#
#	PARESTLIB
#
# 
#	Created by Gabriele Lillacci in December 2011.
#	Latest revision: December 2011.
#
#
#	This free software is available under the Creative Commons Attribution Non-Commercial Share Alike License.
#	You are permitted to use, share, copy, redistribute and adapt this software as long as appropriate credit
#	is given to the original author, all derivative works are distributed under the same license or a compatible one,
#	and this software and its derivatives are not used for commercial purposes.
#	For more information, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or contact
#	Creative Commons, 171 2nd Street, Suite 300, San Francisco, California, 94105, USA.

 
AC_PREREQ(2.59)
AC_INIT(parestlib, 1.0)

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(foreign -Wall -Werror)
AM_PROG_AR

AC_CONFIG_MACRO_DIR([m4])

# Check for C compiler and set it to C99 mode
AC_PROG_CC_C99

# Initialize libtool
LT_INIT()

# Check for required libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AC_CHECK_LIB([stochmod],[syncirc_mod_setup])

# Specify output files
AC_CONFIG_HEADER(config.h)
AC_CONFIG_FILES(Makefile src/Makefile)

AC_OUTPUT