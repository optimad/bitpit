/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
//                  - SORTING ALGORITHMS, EXAMPLES OF USAGE -                 //
//                                                                            //
// Functions for data sorting.                                                //
// ========================================================================== //
// INFO                                                                       //
// Author    : Alessandro Alaia                                               //
// Version   : v2.0                                                           //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <sstream>
# include <iostream>
# include <chrono>
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// Sort algorithms
# include "bitpit_common.hpp"
# include "bitpit_common.hpp"
# include "bitpit_SA.hpp"

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
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
int subtest_001(
    void
) {

// ========================================================================== //
// int subtest_001(                                                           //
//     void)                                                                  //
//                                                                            //
// LIFO stack demo.                                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                            N = 10;
bitpit::LIFOStack<ivector1D>           stack(5);

// Counters
int                            i;

// ========================================================================== //
// OPENING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "================== TEST 000: LIFO stack demo ==================" << endl;

}

// ========================================================================== //
// INSERT ELEMENTS INTO THE LIFO STACK                                        //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    ivector1D        dummy(2, 0);

    // PUSH elements into the LIFO stack ------------------------------------ //
    for (i = 0; i < N; i++) {
        dummy[0] = i+1;    dummy[1] = -dummy[0];
        cout << " * pushing element: " << dummy << endl;
        stack.push(dummy);
        stack.display(cout);
    } //next i

    // POP elements from the LIFO stack ------------------------------------- //
    for (i = 0; i < N; i++) {
        dummy = stack.pop();
        cout << " * popped element: " << dummy << endl;
        stack.display(cout);
    } //next i
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "====================== TEST 000: done!! =======================" << endl;

}

return 0; };

// -------------------------------------------------------------------------- //
int subtest_002(
    void
) {

// ========================================================================== //
// int subtest_002(                                                           //
//     void)                                                                  //
//                                                                            //
// Demo for MinPQueue (min heap queue)                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
vector< array<int,2> >       map(20), *map_ = &map;
bitpit::MinPQueue<double, string>    heap(2, true, map_);

// Counters
// none

// ========================================================================== //
// OPENING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "=================== TEST 001: min heap queue ==================" << endl;
}



// ========================================================================== //
// BUILD MIN-HEAP                                                             //
// ========================================================================== //
{


    // Scope variables ------------------------------------------------------ //
    int            i, N = 10;
    double         dummy_key;
    stringstream   sdummy;
    string         dummy_label;

    for(i=0; i<20; ++i){
        map[i].fill(-1) ;
    };

    // Display min heap ----------------------------------------------------- //
    heap.display(cout);

    // Create min-heap ------------------------------------------------------ //
    cout << " - Creating min-heap by element insertion" << endl;
    for (i = 0; i < N; i++) {

        // Sorting key
        dummy_key = ((double) (N-i));

        // Create label
        sdummy << "i" << i;
        dummy_label = sdummy.str();
        sdummy.str("");
        map[heap.heap_size][0] = i;
        map[i][1] = i;
        heap.insert(dummy_key, dummy_label);

    } //next i
    heap.display(cout);

    // Insert new element(s) ------------------------------------------------ //
    cout << " - Inserting a new element with key: -0.1" << endl;
    dummy_key = -0.1;
    dummy_label = "i10";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 11.0" << endl;
    dummy_key = 11.0;
    dummy_label = "i11";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 9.0" << endl;
    dummy_key = 9.0;
    dummy_label = "i12";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 2.1" << endl;
    dummy_key = 2.1;
    dummy_label = "i13";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);

    // Modify an old element ------------------------------------------------ //
    cout << " - Modifying element with index 4" << endl;
    dummy_key = 4.1;
    dummy_label = "i8";
    heap.modify(4, dummy_key, dummy_label);
    heap.display(cout);
    cout << "mapper: " << map << endl;

    // Extract data from min-heap ------------------------------------------- //
    cout << " - Extracting data from min-heap" << endl;
    while (heap.heap_size > 0) {
        heap.extract(dummy_key, dummy_label);
        cout << "key: " << dummy_key << ", label: " << dummy_label << endl;
    } //next heap item
    heap.display(cout);
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== TEST 001: done!! =====================" << endl;
}

return 0; }

// -------------------------------------------------------------------------- //
int subtest_003(
    void
) {

// ========================================================================== //
// int subtest_003(                                                           //
//     void)                                                                  //
//                                                                            //
// Demo for MaxPQueue (max heap queue)                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
vector< array<int,2> >       map(20), *map_ = &map;
bitpit::MaxPQueue<double, string>    heap(2, true, map_);

// Counters
// none

// ========================================================================== //
// OPENING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "=================== TEST 002: min heap queue ==================" << endl;
}

// ========================================================================== //
// BUILD MIN-HEAP                                                             //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    int            i, N = 10;
    double         dummy_key;
    stringstream   sdummy;
    string         dummy_label;
    

    for( i=0; i<20; ++i){
        map[i].fill(-1) ;
    };

    // Display min heap ----------------------------------------------------- //
    heap.display(cout);

    // Create min-heap ------------------------------------------------------ //
    cout << " - Creating min-heap by element insertion" << endl;
    for (i = 0; i < N; i++) {

        // Sorting key
        dummy_key = ((double) i);

        // Create label
        sdummy << "i" << i;
        dummy_label = sdummy.str();
        sdummy.str("");
        map[heap.heap_size][0] = i;
        map[i][1] = i;
        heap.insert(dummy_key, dummy_label);

    } //next i
    heap.display(cout);

    // Insert new element(s) ------------------------------------------------ //
    cout << " - Inserting a new element with key: -0.1" << endl;
    dummy_key = -0.1;
    dummy_label = "i10";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 11.0" << endl;
    dummy_key = 11.0;
    dummy_label = "i11";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 9.0" << endl;
    dummy_key = 9.0;
    dummy_label = "i12";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);
    cout << " - Inserting a new element with key: 2.1" << endl;
    dummy_key = 2.1;
    dummy_label = "i13";
    map[heap.heap_size][0] = heap.heap_size;
    map[heap.heap_size][1] = heap.heap_size;
    heap.insert(dummy_key, dummy_label);
    heap.display(cout);

    // Modify an old element ------------------------------------------------ //
    cout << " - Modifying element with index 4" << endl;
    dummy_key = 4.1;
    dummy_label = "i8";
    heap.modify(4, dummy_key, dummy_label);
    heap.display(cout);
    cout << "mapper: " << map << endl;

    // Extract data from min-heap ------------------------------------------- //
    cout << " - Extracting data from min-heap" << endl;
    while (heap.heap_size > 0) {
        heap.extract(dummy_key, dummy_label);
        cout << "key: " << dummy_key << ", label: " << dummy_label << endl;
    } //next heap item
    heap.display(cout);
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== TEST 002: done!! =====================" << endl;
}

return 0; }

// -------------------------------------------------------------------------- //
int subtest_004(
    void
) {

// ========================================================================== //
// int subtest_004(                                                           //
//     void)                                                                  //
//                                                                            //
// kd-tree demo.                                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                         N = 5120000;
bitpit::KdTree<2, vector<double> >           KD(N);
dvector2D                   X(N, dvector1D(2, 0.0));
dvector1D                   x(2, -0.1);

// Counters
// none

// ========================================================================== //
// OPENING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "=================== TEST 003: kd-tree =========================" << endl;
}


// ========================================================================== //
// BUILD KD-TREE                                                              //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    chrono::time_point<chrono::system_clock> start, end;
    int                                      i;
    std::chrono::milliseconds                elapsed_seconds;


    // Output message ------------------------------------------------------- //
    cout << " - Generating k-d tree by element insertion" << endl;

    // Create vertex list --------------------------------------------------- //
    std::srand(0);
    for (i = 0; i < N; i++) {
        X[i][0] = ((double) rand())/((double) RAND_MAX);
        X[i][1] = ((double) rand())/((double) RAND_MAX);
    } //next i

    // Build kd tree -------------------------------------------------------- //
    start = std::chrono::system_clock::now();
    for (i = 0; i < N; i++) {
        KD.insert(&X[i]);
    } //next i
    end = chrono::system_clock::now();
    elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start);
    cout << "elapsed time: " << elapsed_seconds.count() << " ms" << endl;

    // output kd-tree ------------------------------------------------------- //
//     for (i = 0; i < KD.n_nodes; i++) {
//         cout << "j: " << i << endl;
//         cout << "target: " << *(KD.nodes[i].object_) << endl;
//         cout << "children: (L) " << KD.nodes[i].lchild_ << " (R) " << KD.nodes[i].rchild_ << endl;
//     } //next i
}

// ========================================================================== //
// INSERT ELEMENT                                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << " - Inserting a new element" << endl;

    // Inserting a new element ---------------------------------------------- //
    KD.insert(&x);
    cout << "# of nodes: " << KD.n_nodes << endl;

    // output kd-tree ------------------------------------------------------- //
//     for (i = 0; i < KD.n_nodes; i++) {
//         cout << "j: " << i << endl;
//         cout << "target: " << *(KD.nodes[i].object_) << endl;
//         cout << "children: (L) " << KD.nodes[i].lchild_ << " (R) " << KD.nodes[i].rchild_ << endl;
//     } //next i

}

// ========================================================================== //
// CHECK FOR EXISTING ELEMENTS                                                //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int            i;

    // Output message ------------------------------------------------------- //
    cout << " - Checking existing elements" << endl;

    // Create vertex list --------------------------------------------------- //

    // Build kd tree -------------------------------------------------------- //
     for (i = N-1; i >= 0; i--) {
         if( KD.exist(&X[i]) < 0 ){
             cout << "error: " << i << "  " << X[i] << endl ;
         };
//ht         cout << "elem found: " << KD.exist(&X[i]) << endl;
     } //next i

}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== TEST 003: done!! =====================" << endl;
}

return 0; }


// -------------------------------------------------------------------------- //
int subtest_005(
    void
) {

// ========================================================================== //
// int subtest_005(                                                           //
//     void)                                                                  //
//                                                                            //
// Sorting routines demo.                                                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err       : int, error flag:                                             //
//               err = 0  --> no error(s)                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                selection;

// Counters
// none

// ========================================================================== //
// OPENING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "=================== TEST 004: kd-tree =========================" << endl;
}

// ========================================================================== //
// EXTRACTION WITHOUT REPLACEMENT                                             //
// ========================================================================== //
selection = 0;
if (selection == 0) {

    // Scope variables ------------------------------------------------------ //
    int               n = 10, N = 9, out;
    ivector1D         list;

    // Initialize variables ------------------------------------------------- //
    std::srand(0);

    // Output message ------------------------------------------------------- //
    cout << endl;
    cout << "Test: random extraction wihtout replacement" << endl;
    cout << "Subtest: random permutation of integers in [0, " << N << "]" << endl;

    // Extract n integers in [0, N] ----------------------------------------- //
    bitpit::utils::extractWithoutReplacement(n, N, list);
    cout << "Extracted: " << list << endl;
    cout << endl;

    // Output message ------------------------------------------------------- //
    cout << "Test: random extraction wihtout replacement" << endl;
    cout << "Subtest: expected value of 1 extraction in [0, " << N << "]" << endl;

    // Extract n integers in [0, N] ----------------------------------------- //
    out = 0;
    for (int I_ = 0; I_ < 100000; I_++) {
        bitpit::utils::extractWithoutReplacement(1, N, list);
        out += list[0];
    } //next I
    cout << "Expectation: " << ((double) out)/(100000.0) << endl;
    cout << endl;
}

// ========================================================================== //
// CLOSING MESSAGE                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== TEST 004: done!! =====================" << endl;
}

return 0; }

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }

        status = subtest_002();
        if (status != 0) {
            return (20 + status);
        }

        status = subtest_003();
        if (status != 0) {
            return (30 + status);
        }

        status = subtest_004();
        if (status != 0) {
            return (40 + status);
        }

        status = subtest_005();
        if (status != 0) {
            return (50 + status);
        }
    } catch (const std::exception &exception) {
        cout << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
