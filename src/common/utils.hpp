#ifndef __BITPIT_UTILS_HPP__
#define __BITPIT_UTILS_HPP__

/*! \file */

#include <array>
#include <functional>
#include <vector>

namespace bitpit {

namespace utils {

template <typename T, typename Comparator = std::less<T> >
bool addToOrderedVector(const T &value, std::vector<T> &list, Comparator comparator = Comparator());

#ifndef __BITPIT_UTILS_SRC__
extern template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
#endif

void extractWithoutReplacement(                                               // Extract integers without replacement
    int                         ,                                             // (input) number of integers to be extracted
    int                         ,                                             // (input) upper bound of extraction interval
    std::vector<int>           &                                              // (input/output) list of extracted value
);

}

}

#endif
