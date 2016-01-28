#define __BITPIT_UTILS_SRC__

#include "utils.hpp"
#include "utils.tpp"

namespace bitpit {

namespace utils {

template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
template bool addToOrderedVector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

}

}
