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

namespace bitpit{

/*!
 \ingroup   SortAlgorithms
 \{
 */

/*!

    \class KdNode
    \brief class for kd-tree node.

    Store the informatino of a node in kd-tree data structure.
    
    Template parameters are:
    - T, container used to store node coordinates
    - T1, label associated to the node.

    Template parameters can be any type fulfilling the following requirements:
    1. operator* (dereferencing operator) must be defined for class T
    2. T1 must be a copy-constructible type.

*/

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR KDNODE                                        //
// ========================================================================== //

// Constructor(s) =========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class KdNode.
    Initialize an empty node in the kd-tree.
*/
template<class T, class T1>
KdNode<T, T1>::KdNode(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET DEFAULT VALUES                                                         //
// ========================================================================== //
lchild_ = -1;
rchild_ = -1;

return; }

// Destructor(s) ============================================================ //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class KdNode.
    Clear KdNode content and release memory.
*/
template<class T, class T1>
KdNode<T, T1>::~KdNode(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET DEFAULT VALUES                                                         //
// ========================================================================== //
lchild_ = -1;
rchild_ = -1;

return; }

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR KDTREE                                        //
// ========================================================================== //

/*!

    \ingroup   SortAlgorithms
    \class KdTree
    \brief class for kd-tree data structure.

    Sort vertices in a d-dimensional Euclidean space into a kd-tree structure.
    
    Template parameters are:
    - d, number dimensions (i.e. number of coordinates for vertices)
    - T, container used to store vertex coordinates
    - T1, label type associated to each node in the kd-tree.

    Template parameters can be any type fulfilling the following requirements:
    1. operator* (dereferencing operator) must be defined for class T
    2. T1 must be a copy-constructible type.
    3. operator== must be defined for container T, which returns true
       if two objects of type T have the same content (i.e. two vertices have the
       same coordinates)

*/

// Constructors ============================================================= //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class KdTree.

    Initialize an empty kd-tree structure and reserve memory for the insertion
    of maxstack nodes.

    \param[in] maxstack memory reserve (in terms of number of elements)
*/
template<int d, class T, class T1>
KdTree<d, T, T1>::KdTree(
    int      maxstack
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET DEFAULT VALUES                                                         //
// ========================================================================== //
MAXSTK = maxstack;
n_nodes = 0;

// ========================================================================== //
// RESIZE NODE LIST                                                           //
// ========================================================================== //
increaseStack();

return; }

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class KdTree.

    Clear kd-tree content and release memory.
*/
template<int d, class T, class T1>
KdTree<d, T, T1>::~KdTree(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// DESTRUCT MEMBERS                                                           //
// ========================================================================== //
n_nodes = 0;
MAXSTK = 0;
nodes.clear();

return; }

// Methods ================================================================== //

// -------------------------------------------------------------------------- //
/*!
    Check whether a given vertex already exist in the kd-tree.
    Check is performed via lexicographical comparison of vertex coordinates.

    \param[in] P_ pointer to container storing the coordinates of the vertex
    to be checked.

    \result on output returns the index of the kd-node in the tree having the
    same coordinates of the input vertex. If no node is found in the kd-tree,
    returns -1.
*/
template<int d, class T, class T1>
int KdTree<d, T, T1>::exist(
    T      *P_
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             check = false;
int              index = -1;
int              prev_ = -1, next_ = 0;
int              lev = 0, dim;

// Counters
// none

// ========================================================================== //
// EXIT FOR EMPTY TREE                                                        //
// ========================================================================== //
if (n_nodes == 0) { return(index); };

// ========================================================================== //
// MOVE ON TREE BRANCHES                                                      //
// ========================================================================== //

// Find node on leaf -------------------------------------------------------- //
while ((next_ >= 0) && (!check)) {
    check = ((*P_) == (*(nodes[next_].object_)));
    prev_ = next_;
    dim = lev % d;
    if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
        next_ = nodes[next_].lchild_;
    }
    else {
        next_ = nodes[next_].rchild_;
    }
    lev++;
} //next
if (check) {
    index = prev_;
}

return(index); };

// -------------------------------------------------------------------------- //
/*!
    Check whether a given vertex already exist in the kd-tree.
    Check is performed via lexicographical comparison of vertex coordinates.

    \param[in] P_ pointer to container storing the coordinates of the vertex
    to be checked.
    \param[in,out] label on output stores the label associated with the KdNode
    (if any) whose coordinates match those of the input vector

    \result on output returns the index of the kd-node in the tree having the
    same coordinates of the input vertex. If no node is found in the kd-tree,
    returns -1.
*/
template<int d, class T, class T1>
int KdTree<d, T, T1>::exist(
    T               *P_,
    T1              &label
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             check = false;
int              index = -1;
int              prev_ = -1, next_ = 0;
int              lev = 0, dim;

// Counters
// none

// ========================================================================== //
// EXIT FOR EMPTY TREE                                                        //
// ========================================================================== //
if (n_nodes == 0) { return(index); };

// ========================================================================== //
// MOVE ON TREE BRANCHES                                                      //
// ========================================================================== //

// Find node on leaf -------------------------------------------------------- //
while ((next_ >= 0) && (!check)) {
    check = ((*P_) == (*(nodes[next_].object_)));
    prev_ = next_;
    dim = lev % d;
    if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
        next_ = nodes[next_].lchild_;
    }
    else {
        next_ = nodes[next_].rchild_;
    }
    lev++;
} //next
if (check) {
    index = prev_;
    label = nodes[prev_].label;
}

return(index); };

// -------------------------------------------------------------------------- //
/*!
    Given an input vertex P, returns the index of the first node encountered
    in the kd-tree which is in the 1-ball centered on P and having a radius of h.
    The 1-ball is defined as:
    B1(x; h) = {y: norm_1(y-x) <= h}
    \param[in] P_ pointer to container storing the coordinates of P
    \param[in] h 1-ball radius
    \param[in] debug (default = false) flag for running the search in debug mode
    \param[in] next_ (default = 0) index of element in the kd tree used as starting
    point by the kd-tree search algorithm
    \param[in] lev   (default = 0) level in kd-tree of node used as starting point
    by the kd-tree search algorithm.

    \result on output returns the index of the first kd-node encountered in the tree
    which lies in the 1-ball centered on P. If no node is found, returns -1.
*/
template<int d, class T, class T1 >
template< class T2>
int KdTree<d, T, T1>::hNeighbor(
    T               *P_,
    T2               h,
    bool             debug,
    int              next_,
    int              lev
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             check = false;
int              index_l = -1, index_r = -1;
int              prev_ = next_;
int              dim;

// Counters
// none

// ========================================================================== //
// EXIT FOR EMPTY TREE                                                        //
// ========================================================================== //
if (n_nodes == 0) { return(-1); };

// ========================================================================== //
// MOVE ON TREE BRANCHES                                                      //
// ========================================================================== //

// Check if root is in the h-neighbor of P_ --------------------------------- //
//if (debug) { cout << "visiting: " << prev_ << endl; }
if (norm2((*(nodes[prev_].object_)) - (*P_)) <= h) { return(prev_); }

// Move on next branch ------------------------------------------------------ //
dim = lev % d;
if (((*(nodes[prev_].object_))[dim] >= (*P_)[dim] - h)
 && (nodes[prev_].lchild_ >= 0)) {
// if (nodes[prev_].lchild_ >= 0) {
    next_ = nodes[prev_].lchild_;
    index_l = hNeighbor(P_, h, debug, next_, lev+1);
}
if (((*(nodes[prev_].object_))[dim] <= (*P_)[dim] + h)
 && (nodes[prev_].rchild_ >= 0)) {
// if (nodes[prev_].rchild_ >= 0) {
    next_ = nodes[prev_].rchild_;
    index_r = hNeighbor(P_, h, debug, next_, lev+1);
}

//if (debug) { cout << "result is: " << max(index_l, index_r) << endl; }
return(std::max(index_l, index_r)); };

// -------------------------------------------------------------------------- //
/*!
    Given an input vertex P, returns the index of the first node encountered
    in the kd-tree which is in the 1-ball centered on P and having a radius of h.
    The 1-ball is defined as:
    B1(x; h) = {y: norm_1(y-x) <= h}
    \param[in] P_ pointer to container storing the coordinates of P
    \param[in,out] label on output stores the label of the KdNode in the 1-ball
    centered on P (if any).
    \param[in] h 1-ball radius
    \param[in] next_ (default = 0) index of element in the kd tree used as starting
    point by the kd-tree search algorithm
    \param[in] lev   (default = 0) level in kd-tree of node used as starting point
    by the kd-tree search algorithm.

    \result on output returns the index of the first kd-node encountered in the tree
    which lies in the 1-ball centered on P. If no node is found, returns -1.
*/
template<int d, class T, class T1 >
template< class T2>
int KdTree<d, T, T1>::hNeighbor(
    T               *P_,
    T1              &label,
    T2               h,
    int              next_,
    int              lev
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             check = false;
int              index_l = -1, index_r = -1;
int              dim;

// Counters
// none

// // ========================================================================== //
// // EXIT FOR EMPTY TREE                                                        //
// // ========================================================================== //
// if (n_nodes == 0) { return(index); };

// // ========================================================================== //
// // MOVE ON TREE BRANCHES                                                      //
// // ========================================================================== //

// // Find node on leaf -------------------------------------------------------- //
// while ((next_ >= 0) && (!check)) {
    // check = (norm2((*P_) - (*(nodes[next_].object_))) <= h);
    // prev_ = next_;
    // dim = lev % d;
    // if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
        // next_ = nodes[next_].lchild_;
    // }
    // else {
        // next_ = nodes[next_].rchild_;
    // }
    // lev++;
// } //next
// if (check) {
    // index = prev_;
    // label = nodes[prev_].label;
// }

return(std::max(index_l, index_r)); };

// -------------------------------------------------------------------------- //
/*!
    Insert a new vertex in the kd-tree.

    \param[in] P_ pointer to container storing vertex coordinates.
*/
template<int d, class T, class T1>
void KdTree<d, T, T1>::insert(
    T               *P_
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             left;
int              prev_ = -1, next_ = 0;
int              lev = 0, dim;

// Counters
// none

// ========================================================================== //
// EXIT CONDITION FOR EMPTY TREE                                              //
// ========================================================================== //
if (n_nodes == 0) {
    nodes[0].object_ = P_;
    n_nodes++;
    return;
};

// ========================================================================== //
// MOVE ON TREE BRANCHES                                                      //
// ========================================================================== //

// Find node on leaf -------------------------------------------------------- //
while (next_ >= 0) {
    prev_ = next_;
    dim = lev % d;
    if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
        next_ = nodes[next_].lchild_;
        left = true;
    }
    else {
        next_ = nodes[next_].rchild_;
        left = false;
    }
    lev++;
} //next

// Insert new element ------------------------------------------------------- //

// Increase stack size
if (n_nodes+1 > nodes.size()) {
    increaseStack();
}

// Update parent
if (left) {
    nodes[prev_].lchild_ = n_nodes;
}
else {
    nodes[prev_].rchild_ = n_nodes;
}

// Insert children
nodes[n_nodes].object_ = P_;
n_nodes++;

return; };


// -------------------------------------------------------------------------- //
/*!
    Insert a new vertex and the associated label into kd-tree.

    \param[in] P_ pointer to container storing vertex coordinates.
    \param[in] label label associated to P_
*/
template<int d, class T, class T1>
void KdTree<d, T, T1>::insert(
    T               *P_,
    T1              &label
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             left;
int              prev_ = -1, next_ = 0;
int              lev = 0, dim;

// Counters
// none

// ========================================================================== //
// QUIT FOR EMPTY TREE                                                        //
// ========================================================================== //
if (n_nodes == 0) {
    nodes[0].object_ = P_;
    nodes[0].label   = label;
    n_nodes++;
    return;
};

// ========================================================================== //
// MOVE ON TREE BRANCHES                                                      //
// ========================================================================== //

// Find node on leaf -------------------------------------------------------- //
while (next_ >= 0) {
    prev_ = next_;
    dim = lev % d;
    if ((*P_)[dim] <= (*(nodes[next_].object_))[dim]) {
        next_ = nodes[next_].lchild_;
        left = true;
    }
    else {
        next_ = nodes[next_].rchild_;
        left = false;
    }
    lev++;
} //next

// Insert new element ------------------------------------------------------- //

// Increase stack size
if (n_nodes+1 > nodes.size()) {
    increaseStack();
}

// Update parent
if (left) {
    nodes[prev_].lchild_ = n_nodes;
}
else {
    nodes[prev_].rchild_ = n_nodes;
}

// Insert children
nodes[n_nodes].object_ = P_;
nodes[n_nodes].label   = label;
n_nodes++;

return; };

// -------------------------------------------------------------------------- //
/*!
    Whenever the kd-tree reaches its full capacity (i.e. the number of KdNode
    stored in the tree is equal to the memory reserved), increase the memory
    reserve by maxstack. The parameters maxstack is set at construction.
*/
template<int d, class T, class T1>
void KdTree<d, T, T1>::increaseStack(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// INCREASE STACK SIZE                                                        //
// ========================================================================== //
nodes.resize(nodes.size() + MAXSTK);

return; };

// -------------------------------------------------------------------------- //
/*!
    Decrease the memory reserved for kd-tree by maxstack.
    The parameters maxstack is set at construction.
*/
template<int d, class T, class T1>
void KdTree<d, T, T1>::decreaseStack(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// INCREASE STACK SIZE                                                        //
// ========================================================================== //
nodes.resize(max(MAXSTK, nodes.size() - MAXSTK));

return; };
/*!
 \}
 */

}
