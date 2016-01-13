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

/*!

    \class kdnode
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
    Default constructor for class kdnode.
    Initialize an empty node in the kd-tree.
*/
template<class T, class T1>
kdnode<T, T1>::kdnode(
    void
) {

// ========================================================================== //
// template<class T, class T1>                                                //
// kdnode<T, T1>::kdnode(                                                     //
//     void)                                                                  //
//                                                                            //
// Default constructor for class kdnode variables                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
    Default destructor for class kdnode.
    Clear kdnode content and release memory.
*/
template<class T, class T1>
kdnode<T, T1>::~kdnode(
    void
) {

// ========================================================================== //
// template<class T, class T1>                                                //
// kdnode<T, T1>::~kdnode(                                                    //
//     void)                                                                  //
//                                                                            //
// Default destructor for class kdnode variables                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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

    \class kdtree
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
    Default constructor for class kdtree.

    Initialize an empty kd-tree structure and reserve memory for the insertion
    of maxstack nodes.

    \param[in] maxstack memory reserve (in terms of number of elements)
*/
template<int d, class T, class T1>
kdtree<d, T, T1>::kdtree(
    int      maxstack
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// kdtree<d, T, T1>::kdtree(                                                  //
//     int      maxstack)                                                     //
//                                                                            //
// Default constructor for kdtree variables.                                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
IncreaseStack();

return; }

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class kdtree.

    Clear kd-tree content and release memory.
*/
template<int d, class T, class T1>
kdtree<d, T, T1>::~kdtree(
    void
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// kdtree<d, T, T1>::~kdtree(                                                 //
//     void)                                                                  //
//                                                                            //
// Default destructor for kdtree variables.                                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
int kdtree<d, T, T1>::exist(
    T      *P_
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// int kdtree<d, T, T1>::exist(                                               //
//     T               *P_)                                                   //
//                                                                            //
// Check if a object is already included in the k-d tree.                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_        : T          , pointer to object to be checked                 //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index     : int, index of kdnode containing object (returns -1 if no     //
//               node is found)                                               //
// ========================================================================== //

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
    \param[in,out] label on output stores the label associated with the kdnode
    (if any) whose coordinates match those of the input vector

    \result on output returns the index of the kd-node in the tree having the
    same coordinates of the input vertex. If no node is found in the kd-tree,
    returns -1.
*/
template<int d, class T, class T1>
int kdtree<d, T, T1>::exist(
    T               *P_,
    T1              &label
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// int kdtree<d, T, T1>::exist(                                               //
//     T               *P_,                                                   //
//     T1              &label)                                                //
//                                                                            //
// Check if a object is already included in the k-d tree.                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_        : T          , pointer to object to be checked                 //
// - label     : T1, label of kd-node matching P_                             //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index     : int, index of kdnode containing object (returns -1 if no     //
//               node is found)                                               //
// ========================================================================== //

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
int kdtree<d, T, T1>::h_neighbor(
    T               *P_,
    T2               h,
    bool             debug,
    int              next_,
    int              lev
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// int kdtree<d, T, T1>::h_neighbor(                                          //
//     T               *P_,                                                   //
//     T                h,                                                    //
//     int              next_,                                                //
//     int              lev)                                                  //
//                                                                            //
// Check if there exist a kdnode in the h-neighborhood of the given item P    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_        : T            pointer to object to be checked                 //
// - h         : T, radius of the ball centered at P_                         //
// - next_     : int, root for binary search                                  //
// - lev       : int, level of root on binary tree                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index     : int, index of kdnode in the h-neighborhood of the given item //
//               (returns -1 if no kd-node is found)                          //
// ========================================================================== //

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
if (norm_2((*(nodes[prev_].object_)) - (*P_)) <= h) { return(prev_); }

// Move on next branch ------------------------------------------------------ //
dim = lev % d;
if (((*(nodes[prev_].object_))[dim] >= (*P_)[dim] - h)
 && (nodes[prev_].lchild_ >= 0)) {
// if (nodes[prev_].lchild_ >= 0) {
    next_ = nodes[prev_].lchild_;
    index_l = h_neighbor(P_, h, debug, next_, lev+1);
}
if (((*(nodes[prev_].object_))[dim] <= (*P_)[dim] + h)
 && (nodes[prev_].rchild_ >= 0)) {
// if (nodes[prev_].rchild_ >= 0) {
    next_ = nodes[prev_].rchild_;
    index_r = h_neighbor(P_, h, debug, next_, lev+1);
}

//if (debug) { cout << "result is: " << max(index_l, index_r) << endl; }
return(max(index_l, index_r)); };

// -------------------------------------------------------------------------- //
/*!
    Given an input vertex P, returns the index of the first node encountered
    in the kd-tree which is in the 1-ball centered on P and having a radius of h.
    The 1-ball is defined as:
    B1(x; h) = {y: norm_1(y-x) <= h}
    \param[in] P_ pointer to container storing the coordinates of P
    \param[in,out] label on output stores the label of the kdnode in the 1-ball
    centered on P (if any).
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
int kdtree<d, T, T1>::h_neighbor(
    T               *P_,
    T1              &label,
    T2               h,
    int              next_,
    int              lev
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// int kdtree<d, T, T1>::h_neighbor(                                          //
//     T               *P_,                                                   //
//     T1              &label,                                                //
//     T2               h)                                                    //
//                                                                            //
// Check if a object is already included in the k-d tree.                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_        : T         *, pointer to object to be checked                 //
// - label     : T1, label of kd-node matching P_                             //
// - h         : T, radius of the ball centered at P_                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - index     : int, index of kdnode containing object (returns -1 if no     //
//               node is found)                                               //
// ========================================================================== //

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
    // check = (norm_2((*P_) - (*(nodes[next_].object_))) <= h);
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

return(max(index_l, index_r)); };

// -------------------------------------------------------------------------- //
/*!
    Insert a new vertex in the kd-tree.

    \param[in] P_ pointer to container storing vertex coordinates.
*/
template<int d, class T, class T1>
void kdtree<d, T, T1>::insert(
    T               *P_
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// void kdtree<d, T, T1>::insert(                                             //
//     T               *P_)                                                   //
//                                                                            //
// Insert a new element in kd-tree.                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_         : T           *, pointer to object                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
    IncreaseStack();
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
void kdtree<d, T, T1>::insert(
    T               *P_,
    T1              &label
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// void kdtree<d, T, T1>::insert(                                             //
//     T               *P_,                                                   //
//     T1              &label)                                                //
//                                                                            //
// Insert a new element in kd-tree.                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P_         : T           *, pointer to object                            //
// - label      : T1, label of new object                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
    IncreaseStack();
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
    Whenever the kd-tree reaches its full capacity (i.e. the number of kdnode
    stored in the tree is equal to the memory reserved), increase the memory
    reserve by maxstack. The parameters maxstack is set at construction.
*/
template<int d, class T, class T1>
void kdtree<d, T, T1>::IncreaseStack(
    void
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// void kdtree<d, T, T1>::IncreaseStack(                                      //
//     void)                                                                  //
//                                                                            //
// Increase stack size.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
void kdtree<d, T, T1>::DecreaseStack(
    void
) {

// ========================================================================== //
// template<int d, class T, class T1>                                         //
// void kdtree<d, T, T1>::DecreaseStack(                                      //
//     void)                                                                  //
//                                                                            //
// Decrease stack size.                                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

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
