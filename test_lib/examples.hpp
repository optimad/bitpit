// ========================================================================== //
//                - Class_SurfTri - Examples of usage                         //
//                                                                            //
// Example of usage for Class_SurfTri (grid manager for unstructured meshes)  //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __EXAMPLES_HPP__
# define __EXAMPLES_HPP__

// ========================================================================== //
// INCLUDE                                                                    //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <sstream>
# include <fstream>
# include <iostream>
# include <chrono>

// CC_lib
# include "Class_SurfTri.hpp"

// Others
// none


// ========================================================================== //
// NAME SPACES                                                                //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
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
// FUNCTION PROTOTYPES                                                        //
// ========================================================================== //
void Demo_CleaningTools(                                                      // Class_SurfTri cleaning tools - demo
    void                                                                      // (input) none
);
void Demo_CleaningTools2(                                                     // Class_SurfTri cleaning tools - demo 2
    void                                                                      // (input) none
);
void Demo_GenerationTools(                                                    // Class_SurfTri generation tools - demo
    void                                                                      // (input) none
);


# endif