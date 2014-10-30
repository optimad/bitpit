// ========================================================================== //
//             - Class_VolTri - Examples of usage                             //
//                                                                            //
// Examples of usage for Class_VolTri unstructured mesh manager.              //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v2.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __EXAMPLES_HPP__
# define __EXAMPLES_HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include<iostream>

// CC_lib
# include "Class_VolTri.hpp"

// Others
// none

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITION                                                           //
// ========================================================================== //

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

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// FUNCTIONS PROTOTYPES                                                       //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Test_000(                                                                // Class_VolTri selection tools demo
    void                                                                      // (input) none
);

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
// none

# endif