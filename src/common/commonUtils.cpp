/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#define __BITPIT_COMMON_UTILS_SRC__

#include <iostream>
#include <cmath>

#include "commonUtils.hpp"

namespace bitpit {

namespace utils {

template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
template bool addToOrderedVector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

template std::vector<long>::const_iterator findInOrderedVector<>(const long&, const std::vector<long>&, std::less<long>);
template std::vector<unsigned long>::const_iterator findInOrderedVector<>(const unsigned long&, const std::vector<unsigned long>&, std::less<unsigned long>);

/*!
* Extract n integers in the interval [0,m] without replacement.
* if n = m+1, returns a random permutation of {0, 1, 2, ..., m}
*
* \param[in] n is the number of extraction
* \param[in] m is the upper bound of extraction interval
* \param[in,out] list is a vector with size n, storing extracted values
*/
void extractWithoutReplacement(int n, int m, std::vector<int> &list)
{
    // Initialize variables
    if (n > m+1) {
        std::cout << "error" << std::endl;
        return;
    }

    // Resize input variables
    list.resize(n);

    // Initialize extraction set
    int N = m;
    std::vector<int> set(m+1, -1);
    for (int i = 0; i < m+1; i++) {
        set[i] = i;
    }

    // Extract integers without replacement
    for (int i = 0; i < n; i++) {
        int index = (int) round(((double) N) * ((double) rand())/((double) RAND_MAX));
        list[i] = set[index];
        set[index] = set[N];
        N--;
    }
}

}

}
