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
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 *  for more details.
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

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR MinPQueue                                     //
// ========================================================================== //

/*!

    \class MinPQueue
    \brief class for min priority queue.

    Class for priority element insertion and extraction. Elements inserted in
    a min. priority queue are internally sorted on a binary tree to ensure the
    following property (heap property):
    given a parent node in the tree (P), and given its children nodes (L and R),
    P < L and P < R.
    In this way, the root element is the one with smallest value in the tree.

    Each time a new element is inserted or removed from the heap, the tree is updated
    by moving the elements upwards or downwards on each branch to maintain the heap
    property.
    The new position of elements in the tree can be tracked by passing a non-null
    pointer to a mapping (a vector< array<int, 2> >) to the heap constructor.
    At any time, the mapping will store the following information:
        - mapping[i][0] stores the index of the node currently stored in the i-th
        position of the tree.
        - mapping[i][1] stores the current position in the tree of the i-th node

    - Template parameter T can be of any copy-constructible type for which operator<
    is defined
    - Template parameter T1 is used for labelling tree nodes, and can be any
    copy-constructible type.
*/

// Constructors ============================================================= //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class MinPQueue.
    Initialize an empty priority queue.
    If a pointer to a vector<array<int, 2>> is provided, the mapping between
    the original and current element position in the tree is tracked.

    \param[in] flag_label flag for using node labelling (true) or not (false).
    \param[in] map_ (default = NULL) pointer to map between the original->current element positions
*/
template <class T, class T1>
MinPQueue<T, T1>::MinPQueue(
        bool                            flag_label,
        std::vector< std::array<int,2> >*map_
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CREATE LIFO STACK                                                          //
    // ========================================================================== //

    // Max stack dimensions
    MAXSTK = 10;

    // Currenst stack size
    heap_size = 0;

    // Set flag
    use_labels = flag_label;

    // Pointers
    map = map_;

    // Initialize stack
    increaseSTACK();

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class MinPQueue.
    Initialize an empty priority queue.
    The size of memory reserve is specified by maxstack.
    If a pointer to a vector<array<int, 2>> is provided, the mapping between
    the original and current element position in the tree is tracked.

    \param[in] maxstack size of memory reserve.
    \param[in] flag_label flag for using node labelling (true) or not (false).
    \param[in] map_ (default = NULL) pointer to map between the original->current element positions
*/
template <class T, class T1>
MinPQueue<T, T1>::MinPQueue(
        int                             maxstack,
        bool                            flag_label,
        std::vector< std::array<int,2> >*map_
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CREATE LIFO STACK                                                          //
    // ========================================================================== //

    // Max stack dimensions
    MAXSTK = maxstack;

    // Currenst stack size
    heap_size = 0;

    // Flags
    use_labels = flag_label;

    // Pointers
    map = map_;

    // Initialize stack
    increaseSTACK();

    return; 
};

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class MinPQueue.
    Clear queue content and free memory.
*/
template <class T, class T1>
MinPQueue<T, T1>::~MinPQueue(
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
    // DESTROY LIFO STACK                                                         //
    // ========================================================================== //

    // Maximal stack size
    MAXSTK = 0;

    // Current stack dimensions
    heap_size = 0;

    // Destroy items
    keys.clear();
    if (use_labels) { labels.clear(); }
    map = NULL;

    // Flags
    use_labels = false;

    return; 
};

// Methods ================================================================== //

// -------------------------------------------------------------------------- //
/*!
    Clear current content without freeing memory.
*/
template <class T, class T1>
void MinPQueue<T, T1>::clear(
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
    // CLEAR CONTENT                                                              //
    // ========================================================================== //
    heap_size = 0;
    keys.resize(MAXSTK);
    if (use_labels) { labels.resize(MAXSTK); }

    return; 
}

// -------------------------------------------------------------------------- //
/*!
    Increase memory reserve by maxstack. The value of the parameter maxstack
    is set at construction.
*/
template <class T, class T1>
void MinPQueue<T, T1>::increaseSTACK(
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

    // stack
    keys.resize(heap_size + MAXSTK);

    // labels
    if (use_labels) {
        labels.resize(heap_size + MAXSTK);
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Decrease memory reserve by maxstack. The value of the parameter maxstack
    is set at construction.
*/
template <class T, class T1>
void MinPQueue<T, T1>::decreaseSTACK(
        void
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int           n = keys.size();

    // Counters
    // none

    // ========================================================================== //
    // DECREASE STACK SIZE                                                        //
    // ========================================================================== //

    // Stack
    keys.resize(std::max(n - MAXSTK, MAXSTK));

    // labels
    if (use_labels) {
        labels.resize(std::max(n - MAXSTK, MAXSTK));
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Restore heap property after element insertion or extraction

    \param[in] i index of element to be moved down/upwards in the tree
*/
template <class T, class T1 >
void MinPQueue<T, T1>::heapify(
        int                       i
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int             L, R, idummy, imap;
    T               dummy;
    T1              dummy2;

    // Counters
    // none

    // ========================================================================== //
    // RESTORE THE MIN HEAP PROPERTY                                              //
    // ========================================================================== //

    // Index of childrens
    L = 2*i + 1;
    R = 2*(i + 1);

    // Check children's key
    if ((L < heap_size) && (keys[L] < keys[i])) {
        idummy = L;
    }
    else {
        idummy = i;
    }
    if ((R < heap_size) && (keys[R] < keys[idummy])) {
        idummy = R;
    }

    // Move up-heap
    if (idummy != i) {
        dummy = keys[i];
        keys[i] = keys[idummy];
        keys[idummy] = dummy;
        if (use_labels) {
            dummy2 = labels[i];
            labels[i] = labels[idummy];
            labels[idummy] = dummy2;
        }
        if (map != NULL) {
            imap = (*map)[i][0];
            (*map)[i][0] = (*map)[idummy][0];
            (*map)[idummy][0] = imap;
            (*map)[(*map)[i][0]][1] = i;
            (*map)[(*map)[idummy][0]][1] = idummy;
        }
        heapify(idummy);
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Build a min heap, from data currently stored in the underlying container.
*/
template <class T, class T1>
void MinPQueue<T, T1>::buildHeap(
        void
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    int                i;

    // ========================================================================== //
    // BUILD MIN HEAP                                                             //
    // ========================================================================== //
    for (i = (heap_size - 1)/2; i >= 0; i--) {
        heapify(i);
    } //next i

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Extract the next element (i.e. the element on the root of the tree) from heap
    list and reduce heap size by one element. After extraction heap property is
    automatically restored by updating the tree.

    \param[in,out] root value of element at tree root.
    \param[in,out] root_label label associated to root element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MinPQueue<T, T1>::extract(
        T                      &root,
        T1                     &root_label
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int        imap;

    // Counters
    // none

    // ========================================================================== //
    // EXTRACT HEAP MIN ELEMENT                                                   //
    // ========================================================================== //

    // Extract min element
    if (heap_size == 0) {
        return;
    }
    root = keys[0];

    // Update heap data structure
    keys[0] = keys[heap_size-1];
    if (use_labels) {
        root_label = labels[0];
        labels[0] = labels[heap_size-1];
    }
    if (map != NULL) {
        imap = (*map)[0][0];
        (*map)[0][0] = (*map)[heap_size-1][0];
        (*map)[heap_size-1][0] = imap;
        (*map)[imap][1] = heap_size-1;
        (*map)[(*map)[0][0]][1] = 0;
    }

    // Reduce stack dimensions
    heap_size--;
    if (heap_size <= keys.size() - MAXSTK) {
        decreaseSTACK();
    }

    // Restore min-heap condition
    heapify(0);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Extract the next element (i.e. the element on the root of the tree) from heap
    list and reduce heap size by one element. After extraction heap property is
    automatically restored by updating the tree.

    \param[in,out] root value of element at tree root.
*/
template <class T, class T1>
void MinPQueue<T, T1>::extract(
        T                      &root
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int        imap;

    // Counters
    // none

    // ========================================================================== //
    // EXTRACT HEAP MIN ELEMENT                                                   //
    // ========================================================================== //

    // Extract min element
    if (heap_size == 0) {
        return;
    }
    root = keys[0];

    // Update heap data structure
    keys[0] = keys[heap_size-1];
    if (map != NULL) {
        imap = (*map)[0][0];
        (*map)[0][0] = (*map)[heap_size-1][0];
        (*map)[heap_size-1][0] = imap;
        (*map)[imap][1] = heap_size-1;
        (*map)[(*map)[0][0]][1] = 0;
    }

    // Reduce stack dimensions
    heap_size--;
    if (heap_size <= keys.size() - MAXSTK) {
        decreaseSTACK();
    }

    // Restore min-heap condition
    heapify(0);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Insert a new element in the heap. After insertion the heap property is
    automatically restored by updating the tree.

    \param[in] key_new value of the new element.
    \param[in] label_new label associated to the new element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MinPQueue<T, T1>::insert(
        T                      &key_new,
        T1                     &label_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // INSERT A NEW KEY                                                           //
    // ========================================================================== //

    // Insert key
    if (heap_size+1 > keys.size()) {
        increaseSTACK();
    }

    // Add new key
    keys[heap_size] = key_new;

    // Update labels
    if (use_labels) {
        labels[heap_size] = label_new;
    }

    // Update heap size
    heap_size++;

    // Restore max heap condition
    modify(heap_size-1, key_new, label_new);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Insert a new element in the heap. After insertion the heap property is
    automatically restored by updating the tree.

    \param[in] key_new value of the new element.
*/
template <class T, class T1>
void MinPQueue<T, T1>::insert(
        T                      &key_new
        ) {

    // ========================================================================== //
    // template <class T, class T1>                                               //
    // void MinPQueue<T, T1>::insert(                                             //
    //     T                      &key_new)                                       //
    //                                                                            //
    // Insert new key into a min heap                                             //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - key_new   : T1, key to be inserted                                       //
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
    // INSERT A NEW KEY                                                           //
    // ========================================================================== //

    // Insert key
    if (heap_size+1 > keys.size()) {
        increaseSTACK();
    }

    // Add new key
    keys[heap_size] = key_new;

    // Update heap size
    heap_size++;

    // Restore max heap condition
    modify(heap_size-1, key_new);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Modify the value of an existing element in the heap. The heap property is
    automatically restored by moving down/upwards the element in the tree.

    \param[in] i index of element in the tree.
    \param[in] key_new new value of the element.
    \param[in] label_new new label to assign to the element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MinPQueue<T, T1>::modify(
        int                     i,
        T                      &key_new,
        T1                     &label_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                 imap;
    T                   dummy1;
    T1                  dummy2;

    // Counters
    int                 j;

    // ========================================================================== //
    // INCREASE VALUE OF LAST KEY.                                                //
    // ========================================================================== //
    if (key_new > keys[i]) {

        // Update keys
        keys[i] = key_new;
        if (use_labels) {
            labels[i] = label_new;
        }

        // move down-heap
        heapify(i);

    }
    else {

        // Update keys
        keys[i] = key_new;
        if (use_labels) {
            labels[i] = label_new;
        }

        // move up-heap
        j = (i + 1)/2 - 1;
        while ((i > 0) && (keys[j] > keys[i])) {
            dummy1 = keys[i];
            keys[i] = keys[j];
            keys[j] = dummy1;
            if (use_labels) {
                dummy2 = labels[i];
                labels[i] = labels[j];
                labels[j] = dummy2;
            }
            if (map != NULL) {
                imap = (*map)[i][0];
                (*map)[i][0] = (*map)[j][0];
                (*map)[j][0] = imap;
                (*map)[(*map)[i][0]][1] = i;
                (*map)[(*map)[j][0]][1] = j;
            }
            i = j;
            j = (i + 1)/2 - 1;
        } //next parent
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Modify the value of an existing element in the heap. The heap property is
    automatically restored by moving down/upwards the element in the tree.

    \param[in] i index of element in the tree.
    \param[in] key_new new value of the element.
*/
template <class T, class T1>
void MinPQueue<T, T1>::modify(
        int                    i,
        T                      &key_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                 imap;
    T                   dummy1;

    // Counters
    int                 j;

    // ========================================================================== //
    // INCREASE VALUE OF SELECTED KEY.                                            //
    // ========================================================================== //
    if (key_new > keys[i]) {

        // Update keys
        keys[i] = key_new;

        // move down-heap
        heapify(i);

    }
    else {

        // Update keys
        keys[i] = key_new;

        // move up-heap
        j = (i + 1)/2 - 1;
        while ((i > 0) && (keys[j] > keys[i])) {
            dummy1 = keys[i];
            keys[i] = keys[j];
            keys[j] = dummy1;
            if (map != NULL) {
                imap = (*map)[i][0];
                (*map)[i][0] = (*map)[j][0];
                (*map)[j][0] = imap;
                (*map)[(*map)[i][0]][1] = i;
                (*map)[(*map)[j][0]][1] = j;
            }
            i = j;
            j = (i + 1)/2 - 1;
        } //next parent
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Display info about the heap to output stream (e.g. current number of element
    stored in the tree, memory reserved, etc.)

    \param[in,out] out output stream
*/
template <class T, class T1>
void MinPQueue<T, T1>::display(
        std::ostream    &out
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    int           i, j, k;

    // ========================================================================== //
    // DISPLAY MIN-HEAP INFO                                                      //
    // ========================================================================== //

    // General info ------------------------------------------------------------- //
    out << "min heap:" << std::endl;
    out << "heap size:         " << heap_size << std::endl;
    out << "data struct. size: " << keys.size() << std::endl;
    out << "max stack size:    " << MAXSTK << std::endl;

    // heap content ------------------------------------------------------------- //
    out << "min heap data:" << std::endl;
    i = 0;
    k = 0;
    while (i < heap_size) {
        j = 0;
        while ((j < pow(2, k)) && (i < heap_size)) {
            out << "[" << keys[i];
            if (use_labels) {
                out << ", '" << labels[i] << "'";
            }
            out << "] ";
            i++;
            j++;
        } //next j
        out << std::endl;
        k++;
    } //next i

    return; 
}

/*!
 \}
 */

/*!
 \ingroup   SortAlgorithms
 \{
 */

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR MaxPQueue                                     //
// ========================================================================== //

/*!

    \class MaxPQueue
    \brief class for max priority queue.

    Class for priority element insertion and extraction. Elements inserted in
    a max. priority queue are internally sorted on a binary tree to ensure the
    following property (heap property):
    given a parent node in the tree (P), and given its children nodes (L and R),
    P > L and P > R.
    In this way, the root element is the one with the largest value in the tree.

    Each time a new element is inserted or removed from the heap, the tree is updated
    by moving the elements upwards or downwards on each branch to maintain the heap
    property.
    The new position of elements in the tree can be tracked by passing a non-null
    pointer to a mapping (a vector< array<int, 2> >) to the heap constructor.
    At any time, the mapping will store the following information:
        - mapping[i][0] stores the index of the node currently stored in the i-th
        position of the tree.
        - mapping[i][1] stores the current position in the tree of the i-th node

    - Template parameter T can be of any copy-constructible type for which operator<
    is defined
    - Template parameter T1 is used for labelling tree nodes, and can be any
    copy-constructible type.
*/

// Constructors ============================================================= //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class MaxPQueue.
    Initialize an empty priority queue.
    If a pointer to a vector<array<int, 2>> is provided, the mapping between
    the original and current element position in the tree is tracked.

    \param[in] flag_label flag for using node labelling (true) or not (false).
    \param[in] map_ (default = NULL) pointer to map between the original->current element positions
*/
template <class T, class T1>
MaxPQueue<T, T1>::MaxPQueue(
        bool                            flag_label,
        std::vector< std::array<int,2> >*map_
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CREATE LIFO STACK                                                          //
    // ========================================================================== //

    // Max stack dimensions
    MAXSTK = 10;

    // Currenst stack size
    heap_size = 0;

    // Set flag
    use_labels = flag_label;

    // Pointers
    map = map_;

    // Initialize stack
    increaseSTACK();

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class MaxPQueue.
    Initialize an empty priority queue.
    The size of memory reserve is specified by maxstack.
    If a pointer to a vector<array<int, 2>> is provided, the mapping between
    the original and current element position in the tree is tracked.

    \param[in] maxstack size of memory reserve.
    \param[in] flag_label flag for using node labelling (true) or not (false).
    \param[in] map_ (default = NULL) pointer to map between the original->current element positions
*/
template <class T, class T1>
MaxPQueue<T, T1>::MaxPQueue(
        int                             maxstack,
        bool                            flag_label,
        std::vector< std::array<int,2> >*map_
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CREATE LIFO STACK                                                          //
    // ========================================================================== //

    // Max stack dimensions
    MAXSTK = maxstack;

    // Currenst stack size
    heap_size = 0;

    // Flags
    use_labels = flag_label;

    // Pointers
    map = map_;

    // Initialize stack
    increaseSTACK();

    return; 
};

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class MaxPQueue.
    Clear queue content and free memory.
*/
template <class T, class T1>
MaxPQueue<T, T1>::~MaxPQueue(
        void
        ) {

    // ========================================================================== //
    // template <class T, class T1>                                               //
    // MaxPQueue<T, T1>::~MaxPQueue(                                              //
    //     void)                                                                  //
    //                                                                            //
    // Standard destructor for max heap.                                          //
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
    // DESTROY LIFO STACK                                                         //
    // ========================================================================== //

    // Maximal stack size
    MAXSTK = 0;

    // Current stack dimensions
    heap_size = 0;

    // Destroy items
    keys.clear();
    if (use_labels) { labels.clear(); }
    map = NULL;

    // Flags
    use_labels = false;

    return; 
};

// Methods ================================================================== //

// -------------------------------------------------------------------------- //
/*!
    Clear current content without freeing memory.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::clear(
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
    // CLEAR CONTENT                                                              //
    // ========================================================================== //
    heap_size = 0;
    keys.resize(MAXSTK);
    if (use_labels) { labels.resize(MAXSTK); }

    return; 
}

// -------------------------------------------------------------------------- //
/*!
    Increase memory reserve by maxstack. The value of the parameter maxstack
    is set at construction.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::increaseSTACK(
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

    // stack
    keys.resize(heap_size + MAXSTK);

    // labels
    if (use_labels) {
        labels.resize(heap_size + MAXSTK);
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Decrease memory reserve by maxstack. The value of the parameter maxstack
    is set at construction.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::decreaseSTACK(
        void
        ) {

    // ========================================================================== //
    // template <class T, class T1>                                               //
    // void MaxPQueue<T, T1>::decreaseSTACK(                                      //
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
    int           n = keys.size();

    // Counters
    // none

    // ========================================================================== //
    // DECREASE STACK SIZE                                                        //
    // ========================================================================== //

    // Stack
    keys.resize(std::max(n - MAXSTK, MAXSTK));

    // labels
    if (use_labels) {
        labels.resize(std::max(n - MAXSTK, MAXSTK));
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Restore heap property after element insertion or extraction

    \param[in] i index of element to be moved down/upwards in the tree
*/
template <class T, class T1 >
void MaxPQueue<T, T1>::heapify(
        int                       i
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int             L, R, idummy, imap;
    T               dummy;
    T1              dummy2;

    // Counters
    // none

    // ========================================================================== //
    // RESTORE THE MIN HEAP PROPERTY                                              //
    // ========================================================================== //

    // Index of childrens
    L = 2*i + 1;
    R = 2*(i + 1);

    // Check children's key
    if ((L < heap_size) && (keys[L] > keys[i])) {
        idummy = L;
    }
    else {
        idummy = i;
    }
    if ((R < heap_size) && (keys[R] > keys[idummy])) {
        idummy = R;
    }

    // Move up-heap
    if (idummy != i) {
        dummy = keys[i];
        keys[i] = keys[idummy];
        keys[idummy] = dummy;
        if (use_labels) {
            dummy2 = labels[i];
            labels[i] = labels[idummy];
            labels[idummy] = dummy2;
        }
        if (map != NULL) {
            imap = (*map)[i][0];
            (*map)[i][0] = (*map)[idummy][0];
            (*map)[idummy][0] = imap;
            (*map)[(*map)[i][0]][1] = i;
            (*map)[(*map)[idummy][0]][1] = idummy;
        }
        heapify(idummy);
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Build a max heap, from data currently stored in the underlying container.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::buildHeap(
        void
        ) {

    // ========================================================================== //
    // template <class T, class T1>                                               //
    // void MaxPQueue<T, T1>::buildHeap(                                         //
    //     void)                                                                  //
    //                                                                            //
    // Build max heap tree.                                                       //
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
    int                i;

    // ========================================================================== //
    // BUILD MAX HEAP                                                             //
    // ========================================================================== //
    for (i = (heap_size - 1)/2; i >= 0; i--) {
        heapify(i);
    } //next i

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Extract the next element (i.e. the element on the root of the tree) from heap
    list and reduce heap size by one element. After extraction heap property is
    automatically restored by updating the tree.

    \param[in,out] root value of element at tree root.
    \param[in,out] root_label label associated to root element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MaxPQueue<T, T1>::extract(
        T                      &root,
        T1                     &root_label
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int        imap;

    // Counters
    // none

    // ========================================================================== //
    // EXTRACT HEAP MAX ELEMENT                                                   //
    // ========================================================================== //

    // Extract max element
    if (heap_size == 0) {
        return;
    }
    root = keys[0];

    // Update heap data structure
    keys[0] = keys[heap_size-1];
    if (use_labels) {
        root_label = labels[0];
        labels[0] = labels[heap_size-1];
    }
    if (map != NULL) {
        imap = (*map)[0][0];
        (*map)[0][0] = (*map)[heap_size-1][0];
        (*map)[heap_size-1][0] = imap;
        (*map)[imap][1] = heap_size-1;
        (*map)[(*map)[0][0]][1] = 0;
    }

    // Reduce stack dimensions
    heap_size--;
    if (heap_size <= keys.size() - MAXSTK) {
        decreaseSTACK();
    }

    // Restore max-heap condition
    heapify(0);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Extract the next element (i.e. the element on the root of the tree) from heap
    list and reduce heap size by one element. After extraction heap property is
    automatically restored by updating the tree.

    \param[in,out] root value of element at tree root.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::extract(
        T                      &root
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int        imap;

    // Counters
    // none

    // ========================================================================== //
    // EXTRACT HEAP MAX ELEMENT                                                   //
    // ========================================================================== //

    // Extract max element
    if (heap_size == 0) {
        return;
    }
    root = keys[0];

    // Update heap data structure
    keys[0] = keys[heap_size-1];
    if (map != NULL) {
        imap = (*map)[0][0];
        (*map)[0][0] = (*map)[heap_size-1][0];
        (*map)[heap_size-1][0] = imap;
        (*map)[imap][1] = heap_size-1;
        (*map)[(*map)[0][0]][1] = 0;
    }

    // Reduce stack dimensions
    heap_size--;
    if (heap_size <= keys.size() - MAXSTK) {
        decreaseSTACK();
    }

    // Restore max-heap condition
    heapify(0);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Insert a new element in the heap. After insertion the heap property is
    automatically restored by updating the tree.

    \param[in] key_new value of the new element.
    \param[in] label_new label associated to the new element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MaxPQueue<T, T1>::insert(
        T                      &key_new,
        T1                     &label_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // INSERT A NEW KEY                                                           //
    // ========================================================================== //

    // Insert key
    if (heap_size+1 > keys.size()) {
        increaseSTACK();
    }

    // Add new key
    keys[heap_size] = key_new;

    // Update labels
    if (use_labels) {
        labels[heap_size] = label_new;
    }

    // Update heap size
    heap_size++;

    // Restore max heap condition
    modify(heap_size-1, key_new, label_new);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Insert a new element in the heap. After insertion the heap property is
    automatically restored by updating the tree.

    \param[in] key_new value of the new element.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::insert(
        T                      &key_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // INSERT A NEW KEY                                                           //
    // ========================================================================== //

    // Insert key
    if (heap_size+1 > keys.size()) {
        increaseSTACK();
    }

    // Add new key
    keys[heap_size] = key_new;

    // Update heap size
    heap_size++;

    // Restore max heap condition
    modify(heap_size-1, key_new);

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Modify the value of an existing element in the heap. The heap property is
    automatically restored by moving down/upwards the element in the tree.

    \param[in] i index of element in the tree.
    \param[in] key_new new value of the element.
    \param[in] label_new new label to assign to the element (available only if
    flag_label is set to 'true' at heap construction)
*/
template <class T, class T1>
void MaxPQueue<T, T1>::modify(
        int                     i,
        T                      &key_new,
        T1                     &label_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                 imap;
    T                   dummy1;
    T1                  dummy2;

    // Counters
    int                 j;

    // ========================================================================== //
    // INCREASE VALUE OF LAST KEY.                                                //
    // ========================================================================== //
    if (key_new < keys[i]) {

        // Update keys
        keys[i] = key_new;
        if (use_labels) {
            labels[i] = label_new;
        }

        // move down-heap
        heapify(i);

    }
    else {

        // Update keys
        keys[i] = key_new;
        if (use_labels) {
            labels[i] = label_new;
        }

        // move up-heap
        j = (i + 1)/2 - 1;
        while ((i > 0) && (keys[j] < keys[i])) {
            dummy1 = keys[i];
            keys[i] = keys[j];
            keys[j] = dummy1;
            if (use_labels) {
                dummy2 = labels[i];
                labels[i] = labels[j];
                labels[j] = dummy2;
            }
            if (map != NULL) {
                imap = (*map)[i][0];
                (*map)[i][0] = (*map)[j][0];
                (*map)[j][0] = imap;
                (*map)[(*map)[i][0]][1] = i;
                (*map)[(*map)[j][0]][1] = j;
            }
            i = j;
            j = (i + 1)/2 - 1;
        } //next parent
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Modify the value of an existing element in the heap. The heap property is
    automatically restored by moving down/upwards the element in the tree.

    \param[in] i index of element in the tree.
    \param[in] key_new new value of the element.
*/
template <class T, class T1>
void MaxPQueue<T, T1>::modify(
        int                    i,
        T                      &key_new
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int                 imap;
    T                   dummy1;

    // Counters
    int                 j;

    // ========================================================================== //
    // INCREASE VALUE OF SELECTED KEY.                                            //
    // ========================================================================== //
    if (key_new < keys[i]) {

        // Update keys
        keys[i] = key_new;

        // move down-heap
        heapify(i);

    }
    else {

        // Update keys
        keys[i] = key_new;

        // move up-heap
        j = (i + 1)/2 - 1;
        while ((i > 0) && (keys[j] < keys[i])) {
            dummy1 = keys[i];
            keys[i] = keys[j];
            keys[j] = dummy1;
            if (map != NULL) {
                imap = (*map)[i][0];
                (*map)[i][0] = (*map)[j][0];
                (*map)[j][0] = imap;
                (*map)[(*map)[i][0]][1] = i;
                (*map)[(*map)[j][0]][1] = j;
            }
            i = j;
            j = (i + 1)/2 - 1;
        } //next parent
    }

    return; 
};

// -------------------------------------------------------------------------- //
/*!
    Display info about the heap to output stream (e.g. current number of element
    stored in the tree, memory reserved, etc.)

    \param[in,out] out output stream
*/
template <class T, class T1>
void MaxPQueue<T, T1>::display(
        std::ostream    &out
        ) {

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    int           i, j, k;

    // ========================================================================== //
    // DISPLAY MAX-HEAP INFO                                                      //
    // ========================================================================== //

    // General info ------------------------------------------------------------- //
    out << "max heap:" << std::endl;
    out << "heap size:         " << heap_size << std::endl;
    out << "data struct. size: " << keys.size() << std::endl;
    out << "max stack size:    " << MAXSTK << std::endl;

    // heap content ------------------------------------------------------------- //
    out << "max heap data:" << std::endl;
    i = 0;
    k = 0;
    while (i < heap_size) {
        j = 0;
        while ((j < pow(2, k)) && (i < heap_size)) {
            out << "[" << keys[i];
            if (use_labels) {
                out << ", '" << labels[i] << "'";
            }
            out << "] ";
            i++;
            j++;
        } //next j
        out << std::endl;
        k++;
    } //next i

    return; 
}
/*!
 \}
 */
}
