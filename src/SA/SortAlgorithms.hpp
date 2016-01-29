/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by by the Free Software Foundation.
 *
 *  BitPit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

// ========================================================================== //
//                         - SORTING ALGORITHMS -                             //
//                                                                            //
// Functions for data sorting.                                                //
// ========================================================================== //
// INFO                                                                       //
// Author    : Alessandro Alaia                                               //
// Version   : v2.0                                                           //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __SORT_ALG_HH__
# define __SORT_ALG_HH__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <iostream>

// Classes
// none

// BitPit
# include "Operators.hpp"

namespace bitpit{

// KdTree ------------------------------------------------------------------- //
template <class T, class T1 = int>
class KdNode {

    // Members ========================================================== //
    public:
    int        lchild_;                                                   // pointer to left child
    int        rchild_;                                                   // pointer to left child
    T         *object_;                                                   // pointer to object
    T1         label;

    // Constructor ====================================================== //
    public:
    KdNode(                                                               // default constructor for KdNode variables
        void                                                              // (input) none
    );

    // Destructor ======================================================= //
    public:
    ~KdNode(                                                              // default destructor for KdNode variables
        void                                                              // (input) none
    );
};

template <int d, class T, class T1 = int>
class KdTree {

    // Members ============================================================== //
    public:
    int                                 MAXSTK;                               // max stack size
    int                                 n_nodes;                              // number of nodes
    std::vector< KdNode<T, T1> >        nodes;                                // kd-tree nodes

    // Constructors ========================================================= //
    public:
    KdTree(                                                                   // Default constructor for KdTree
        int                 stack_size = 10                                   // (input/optional stack size)
    );

    // Destructors ========================================================== //
    public:
    ~KdTree(                                                                  // default destructor for KdTree variables
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    int exist(                                                                // Check if element exist in the kd-tree
        T           *                                                         // (input) pointer to element to be tested
    );
    int exist(                                                                // Check if element exist in the kd-tree
        T           *,                                                        // (input) pointer to element to be tested
        T1          &                                                         // (input/output) label of the kd node matching test object
    );

    template <class T2>
    int hNeighbor(                                                           // Check if a kd-node exists in the h-neighborhood of a given item
        T           *,                                                        // (input) pointer to element to be tested
        T2           ,                                                        // (input) radius of ball
        bool         ,
        int         n = 0,                                                    // (input/optional) root for binary search algorithm
        int         l = 0                                                     // (input/optional) level of root on binary tree
    );

    template <class T2>
    int hNeighbor(                                                           // Check if a kd-node exists in the h-neighborhood of a given item
        T           *,                                                        // (input) pointer to element to be tested
        T1          &,                                                      // (input/output) label of the kd node matching test object
        T2           ,                                                        // (input) radius of ball
        int         n = 0,                                                    // (input/optional) root for binary search algorithm
        int         l = 0                                                     // (input/optional) level of root on binary tree
    );
    void insert(                                                              // Insert new element in the kd-tree
        T           *                                                         // (input) pointer to element to be inserted
    );
    void insert(                                                              // Insert new element in the kd-tree
        T           *,                                                        // (input) pointer to element to be inserted
        T1          &                                                         // (input) label of the new element
    );
    private:
    void increaseStack(                                                       // Increase stack size
        void                                                                  // )input) none
    );
    void decreaseStack(                                                       // Decrease stack size
        void                                                                  // )input) none
    );

};

// min PQUEUE --------------------------------------------------------------- //
template <class T, class T1 = T>
class MinPQueue {

    // Members ============================================================== //
    public:
    int                                          MAXSTK;                      // Maximal stack size between resize
    int                                          heap_size;                   // number of elements in stack
    std::vector< T >                             keys;                        // stack
    std::vector< T1 >                            labels;                      // labels associated to keys
    std::vector< std::array<int,2> >            *map;                         // pointer to mapper

    private:
    bool                                        use_labels;                   // flag for key labelling

    // Constructor ========================================================== //
    public:
    MinPQueue(                                                                // Default constructor for min priority queue
        bool                                    a = false,                    // (input/optional) flag for key labelling
        std::vector< std::array<int,2> >       *b = NULL                      // (input/optional) pointer to user-defined map
    );
    MinPQueue(                                                                // Default constructor for min priority queue
        int                                      ,                            // (input) stack size
        bool                                    a = false,                    // (input/optional) flag for key labelling
        std::vector< std::array<int,2> >       *b = NULL                      // (input/optional) pointer to user-defined map
    );

    // Destructor =========================================================== //
    public:
    ~MinPQueue(                                                               // Standard destructor for min priority queues
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear min heap content
        void                                                                  // (input) none
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &                                                     // (input/output) root value
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &,                                                    // (input/output) root value
        T1              &                                                     // (input/output) root label
    );
    void insert(                                                              // Insert a new key
        T               &                                                     // (input) new key to be inserted
    );
    void insert(                                                              // Insert a new key
        T               &,                                                    // (input) new key to be inserted
        T1              &                                                     // (input) label
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &                                                     // (input) new key value
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &,                                                    // (input) new key value
        T1              &                                                     // (input) label attached to the new key
    );
    void buildHeap(                                                          // Build min-heap
        void                                                                  // (input) none
    );
    void display(                                                             // Display min-heap content
        std::ostream    &                                                     // (input) output stream
    );
    private:
    void increaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void decreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (input) none
    );
    void heapify(                                                             // Restore min heap condition on spacified element
        int                                                                   // (input) position of element in stack
    );

};

// max PQUEUE --------------------------------------------------------------- //
template <class T, class T1 = T>
class MaxPQueue {

    // Members ============================================================== //
    public:
    int                                         MAXSTK;                       // Maximal stack size between resize
    int                                         heap_size;                    // number of elements in stack
    std::vector< T >                            keys;                         // stack
    std::vector< T1 >                           labels;                       // labels associated to keys
    std::vector< std::array<int,2> >           *map;                          // pointer to mapper
    private:
    bool                                        use_labels;                   // flag for key labelling

    // Constructor ========================================================== //
    public:
    MaxPQueue(                                                                // Default constructor for min priority queue
        bool                                    a = false,                    // (input/optional) flag for key labelling
        std::vector< std::array<int,2> >       *b = NULL                      // (input/optional) pointer to user-defined map
    );
    MaxPQueue(                                                                // Default constructor for min priority queue
        int              ,                                                    // (input) stack size
        bool                                    a = false,                    // (input/optional) flag for key labelling
        std::vector< std::array<int,2> >       *b = NULL                      // (input/optional) pointer to user-defined map
    );

    // Destructor =========================================================== //
    public:
    ~MaxPQueue(                                                               // Standard destructor for min priority queues
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear max heap content
        void                                                                  // (input) none
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &                                                     // (input/output) root value
    );
    void extract(                                                             // Extract root from the min-heap data structure
        T               &,                                                    // (input/output) root value
        T1              &                                                     // (input/output) root label
    );
    void insert(                                                              // Insert a new key
        T               &                                                     // (input) new key to be inserted
    );
    void insert(                                                              // Insert a new key
        T               &,                                                    // (input) new key to be inserted
        T1              &                                                     // (input) label
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &                                                     // (input) new key value
    );
    void modify(                                                              // Modify key value
        int              ,                                                    // (input) index of key to be modified
        T               &,                                                    // (input) new key value
        T1              &                                                     // (input) label attached to the new key
    );
    void buildHeap(                                                          // Build min-heap
        void
    );
    void display(                                                             // Display min-heap content
        std::ostream         &                                                     // (input) output stream
    );
    private:
    void increaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void decreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (input) none
    );
    void heapify(                                                             // Restore min heap condition on spacified element
        int                                                                   // (input) position of element in stack
    );

};

// LIFO stack --------------------------------------------------------------- //
template <class T>
class LIFOStack {

    // Members ============================================================== //
    public:
    int                 MAXSTK;                                               // Maximal stack size between resize
    int                 TOPSTK;                                               // Current stack size
    std::vector<T>      STACK;                                                // LIFO stack

    // Constructor ========================================================== //
    public:
    LIFOStack(                                                                // Standard constructor for LIFO stack
        void                                                                  // (input) none
    );
    LIFOStack(                                                                // Custom constructor #1 for LIFO stack
        int                                                                   // (input) maximal stack size
    );
    LIFOStack(                                                                // Custom constructor #2 for LIFO stack
        std::vector<T>  &                                                     // (input) items to be added in the LIFO stack
    );

    // Destructor =========================================================== //
    public:
    ~LIFOStack(                                                               // Standard destructor for LIFO stack
        void                                                                  // (input) none
    );

    // Methods ============================================================== //
    public:
    void clear(                                                               // Clear stack content
        void                                                                  // (input) none
    );
    void increaseSTACK(                                                       // Increase stack size
        void                                                                  // (input) none
    );
    void decreaseSTACK(                                                       // Decrease stack size
        void                                                                  // (output) none
    );
    T pop(                                                                    // Pop last item from stack
        void                                                                  // (input) none
    );
    void push(                                                                // Pusk item into the stack
        T                                                                     // (input) item to be pushed into the stack list
    );
    void push(                                                                // Pusk items into the stack
        std::vector<T>  &                                                     // (input) items to be pushed into the stack list
    );
    void display(                                                             // Display LIFO stack infos
        std::ostream    &                                                     // (input) output stream
    );
};

// ========================================================================== //
// FUNCTIONS PROTOTYPES                                                       //
// ========================================================================== //
//None

}

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "LIFOStack.tpp"
# include "PQueue.tpp"
# include "KdTree.tpp"



# endif
