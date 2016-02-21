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
// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                   //
// ========================================================================== //

/*!
    \class LIFOStack
    \brief class for Last In First Out stack

    Manage insertion and extraction from/to a  list of object
    using the L.I.F.O. pragma.

    Class template parameter, T, can be of any copy-constructible type.
*/

// Constructors ============================================================= //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class LIFOStack. Initalize an empty stack
*/
template <class T>
LIFOStack<T>::LIFOStack(
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
// CREATE LIFO STACK                                                          //
// ========================================================================== //

// Max stack dimensions
MAXSTK = 10;

// Initialize stack
increaseSTACK();

// Currenst stack size
TOPSTK = 0;

return; };

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class LIFOStack.
    Initialize an empty stack, and set the memory reserve to maxstack.
    Once the max capacity is reached, insertion of another elements will cause
    the memory to be further increased by maxstack.

    \param[in] maxstack memory reserved for stack.
*/
template <class T>
LIFOStack<T>::LIFOStack(
    int         maxstack
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

// Initialize stack
increaseSTACK();

// Currenst stack size
TOPSTK = 0;

return; };

// -------------------------------------------------------------------------- //
/*!
    Constructor #2 for LIFOStack.
    Initialize a LIFO stack from a vector of elements.
    The memory reserve is set to the size of the input vector and elements
    are inserted in order from items.begin() to items.end()

    \param[in] items list of objects to be inserted into the stack
*/
template <class T>
LIFOStack<T>::LIFOStack(
    std::vector<T>  &items
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
MAXSTK = items.size();

// Initialize stack
increaseSTACK();

// Currenst stack size
TOPSTK = 0;

// Put items into the stack
push(items);

return; };

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
/*!
    Default destructor for class LIFOStack.

    Clear stack content and release memory.
*/
template <class T>
LIFOStack<T>::~LIFOStack(
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
TOPSTK = 0;

// Destroy items
STACK.clear();

return; };

// Methods ================================================================== //

// -------------------------------------------------------------------------- //
/*!
    Clear stack content without freeing memory.
*/
template <class T>
void LIFOStack<T>::clear(
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
TOPSTK = 0;
STACK.resize(MAXSTK);

return; }

// -------------------------------------------------------------------------- //
/*!
    Increase memory reserved for stack by maxstack. The parameter maxstack
    is assigned at class declaration.
*/
template <class T>
void LIFOStack<T>::increaseSTACK(
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
STACK.resize(STACK.size() + MAXSTK);

return; };

// -------------------------------------------------------------------------- //
/*!
    Decrease memory reserved for stack by maxstack. The parameter maxstack
    is assigned at class declaration.
*/
template <class T>
void LIFOStack<T>::decreaseSTACK(
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
// DECREASE STACK SIZE                                                        //
// ========================================================================== //
STACK.resize(std::max((int) (STACK.size() - MAXSTK), MAXSTK));

return; };

// -------------------------------------------------------------------------- //
/*!
    Pop the first item (i.e. the last added item) from stack, 
    reducing its size by one element.

    \result returns the first item in the stack
*/
template <class T>
T LIFOStack<T>::pop(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
T       item;

// Counters
// none

// ========================================================================== //
// POP LAST ELEMENT FROM THE LIFO STACK                                       //
// ========================================================================== //

// Extract last element from the stack list --------------------------------- //
item = STACK[TOPSTK - 1];
TOPSTK--;

// Resize stack list -------------------------------------------------------- //
if (TOPSTK <= (long) STACK.size() - MAXSTK) {
    decreaseSTACK();
}

return(item); };

// -------------------------------------------------------------------------- //
/*!
    Push item into stack. This item will be the first to be returned if pop()
    is called.

    \param[in] item item to be pushed into the list.
*/
template <class T>
void LIFOStack<T>::push(
    T item
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// POP LAST ELEMENT FROM THE LIFO STACK                                       //
// ========================================================================== //

// Resize stack list -------------------------------------------------------- //
if (TOPSTK >= (long) STACK.size()) {
    increaseSTACK();
}

// Extract last element from the stack list --------------------------------- //
STACK[TOPSTK] = item;
TOPSTK++;

return; };

// -------------------------------------------------------------------------- //
/*!
    Push a list of items into stack.
    Items are inserted in the same order they appear in the input vector, i.e.
    from items.begin() to items.end().

    \param[in] items vector storing the list of items to be pushed into
    the stack
*/
template <class T>
void LIFOStack<T>::push(
    std::vector<T> &items
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            n_items;

// Counters
int            i;

// ========================================================================== //
// PUSH ITEMS INTO THE STACK                                                  //
// ========================================================================== //
n_items = items.size();
for (i = 0; i < n_items; i++) {
    push(items[i]);
} //next i

return; };

// -------------------------------------------------------------------------- //
/*!
    Display stack info to output stream (e.g. number of elements currently
    stored into stack, memory currently reserved for stack, etc.)

    \param[in,out] out output stream
*/
template <class T>
void LIFOStack<T>::display(
    std::ostream &out
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// DISPLAY LIFO INFO                                                          //
// ========================================================================== //
out << "LIFO stack:" << std::endl;
out << "  n. of elements:      " << TOPSTK          << std::endl;
out << "  stack buffer:        " << MAXSTK          << std::endl;
out << "  data struct. size:   " << STACK.size()    << std::endl;
out << "  data: " << std::endl;
out << "[";
if (TOPSTK == 0) return;
for (int i = 0; i < TOPSTK-1; i++) {
    out << STACK[i] << ", ";
} //next i
out << STACK[TOPSTK-1] << "]" << std::endl;

return; };
/*!
 \}
 */

}
