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
# ifndef __EXAMPLES_HPP__
# define __EXAMPLES_HPP__

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
