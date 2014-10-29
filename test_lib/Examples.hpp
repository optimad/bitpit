// ========================================================================== //
//                - SORTING ALGORITHMS - EXAMPLES OF USAGE -                  //
//                                                                            //
// Functions for data sorting.                                                //
// ========================================================================== //
// INFO                                                                       //
// Author    : Alessandro Alaia                                               //
// Version   : v2.0                                                           //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __EXAMPLES_HH__
# define __EXAMPLES_HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <sstream>
# include <iostream>
# include <chrono>

// Classes
// none

// Sort algorithms
# include "SortAlgorithms.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
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
void Test000(                                                                 // LIFO stack demo
    void                                                                      // none
);
void Test001(                                                                 // FIFO stack demo
    void                                                                      // none
);
void Test002(                                                                 // minPQUEUE demo
    void                                                                      // none
);
void Test003(                                                                 // maxPQUEUE demo
    void                                                                      // none
);
void Test004(                                                                 // kd tree demo
    void                                                                      // none
);
void Test005(                                                                 // sort-algorithms demo
    void                                                                      // none
);

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
// none

# endif