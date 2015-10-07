// ================================================================================== //
//                         OPERATORS - EXAMPLES OF USAGE -                            //
//                                                                                    //
// Examples of Operators usage.                                                       //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //
# ifndef __OPERATORS_EXAMPLES_HH__
# define __OPERATORS_EXAMPLES_HH__

// ================================================================================== //
// INCLUDES                                                                           //
// ================================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <iostream>

// BitPib Library
# include "Operators.hpp"

// ================================================================================== //
// NAMESPACES                                                                         //
// ================================================================================== //
using namespace std;

// ================================================================================== //
// TYPES DEFINITIONS                                                                  //
// ================================================================================== //

// boolean vectors
typedef vector< bool >                 bvector1D;
typedef vector< bvector1D >            bvector2D;
typedef vector< bvector2D >            bvector3D;
typedef vector< bvector3D >            bvector4D;

// characters vectors
typedef vector< char >                 cvector1D;
typedef vector< cvector1D >            cvector2D;
typedef vector< cvector2D >            cvector3D;
typedef vector< cvector3D >            cvector4D;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// ================================================================================== //
// FUNCTIONS PROTOTYPES                                                               //
// ================================================================================== //
void vectorOperators_Ex(                                                              // Demo for operators between vectors
    void                                                                              // (input) none
);
void vectorMathFunct_Ex(                                                              // Demo for math operators between vectors
    void                                                                              // (input) none
);
void arrayOperators_Ex(                                                               // Demo for operators between vectors
    void                                                                              // (input) none
);
void arrayMathFunct_Ex(                                                               // Demo for operators between vectors
    void                                                                              // (input) none
);
void stringOperators_Ex(                                                              // Demo for operators between strings
    void                                                                              // (input) none
);

// ================================================================================== //
// TEMPLATES                                                                          //
// ================================================================================== //
// none

# endif