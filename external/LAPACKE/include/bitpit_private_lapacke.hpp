/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
 \*---------------------------------------------------------------------------*/

#ifndef __BITPIT_PRIVATE_LAPACKE_HPP__
#define __BITPIT_PRIVATE_LAPACKE_HPP__

/*! \file */

/*!
 * There is a bug in old gcc versions (maybe versions prior to the 4.8.3) when
 * including the file "complex.h" in C++ code[1]. A workaround for this bug
 * is to include the "<complex>" header and define some lapacke macros before
 * including the lapacke header.
 *
 * [1] https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59087
 */
#if defined(__clang__)
#    include <lapacke.h>
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#    include <lapacke.h>
#elif defined(__GNUC__) || defined(__GNUG__)
#    if __GNUC__ == 4 && __GNUC_MINOR__ <= 8
#        include <complex>
#        define lapack_complex_float std::complex<float>
#        define lapack_complex_double std::complex<double>
#        include <lapacke.h>
#    else
#        include <lapacke.h>
#    endif
#else
#    include <lapacke.h>
#endif

#endif
