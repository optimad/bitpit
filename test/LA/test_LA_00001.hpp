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

// ========================================================================== //
//                         LINEAR ALGEBRA PACKAGE                             //
//                                                                            //
// Demo for linear algebra routines.                                          //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author            : Alessandro Alaia                                       //
// Data              : Sept 26, 2014                                          //
// Version           : v2.0                                                   //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __BITPIT_TEST_LA_00001_HPP__
# define __BITPIT_TEST_LA_00001_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <iostream>

// CC_lin
# include "LinearAlgebra.hpp"

// Others
// none

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //

// Double vectors
typedef vector<double>                  dvector1D;
typedef vector<dvector1D>               dvector2D;
typedef vector<dvector2D>               dvector3D;
typedef vector<dvector3D>               dvector4D;

// Integer vectors
typedef vector<int>                     ivector1D;
typedef vector<ivector1D>               ivector2D;
typedef vector<ivector2D>               ivector3D;
typedef vector<ivector3D>               ivector4D;

// ========================================================================== //
// FUNCTIONS PROTOTYPES                                                       //
// ========================================================================== //
void Test_000(                                                                // Basic matrix template demo
    void                                                                      // (input) none
);
void Test_001(                                                                // Basic matrix operations demo
    void                                                                      // (input) none
);
void Test_002(                                                                // Basic matrix manipulation demo
    void                                                                      // (input) none
);
void Test_003(                                                                // Basic linear system solver demo
    void                                                                      // (input) none
);

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
// none

# endif
