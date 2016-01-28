#define __BITPIT_UTILS_SRC__

#include <iostream>

#include "utils.hpp"
#include "utils.tpp"

namespace bitpit {

namespace utils {

template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
template bool addToOrderedVector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

/*!
    \ingroup commonUtils
    \{
*/

/*!
    Extract n integers in the interval [0,m] without replacement.
    if n = m+1, returns a random permutation of {0, 1, 2, ..., m}

    \param[in] n number of extraction
    \param[in] m upper bound of extraction interval
    \param[in,out] list vector with size n, storing extracted values
*/
void extractWithoutReplacement(
    int                 n,
    int                 m,
    std::vector<int>    &list
) {

// ========================================================================== //
// void Extract_wo_Repl(                                                      //
//     int                 n,                                                 //
//     int                 m,                                                 //
//     vector<int>        &list)                                              //
//                                                                            //
// Extract n integers in the interval [0, m] without replacement.             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - n           : int, number of samples to be extracted                     //
// - m           : int, upper bound of extraction interval                    //
// - list        : ivector1D, list of extracted values                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                 N;
std::vector<int>    set(m+1, -1);

// Counters
int                 i, index;

// ========================================================================== //
// INITIALIZE VARIABLES                                                       //
// ========================================================================== //
if (n > m+1) { std::cout << "error" << std::endl; return; }

// Resize input variables
list.resize(n);

// Initialize extraction set
N = m;
for (i = 0; i < m+1; i++) {
    set[i] = i;
} //next i

// ========================================================================== //
// EXTRACT INTEGERS WITHOUT REPLACEMENT                                       //
// ========================================================================== //
for (i = 0; i < n; i++) {
    index = (int) round(((double) N) * ((double) rand())/((double) RAND_MAX));
    list[i] = set[index];
    set[index] = set[N];
    N--;
} //next i


return; }

/*!
    \}
*/

}

}
