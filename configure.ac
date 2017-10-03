#	Autoconf script for libParEstLib
#
# This file is part of libParEstLib.
# Copyright 2011-2017 Gabriele Lillacci.
#
# libParEstLib is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libParEstLib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libParEstLib.  If not, see <http://www.gnu.org/licenses/>.


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
