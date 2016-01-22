#define __BITPIT_UTILS_SRC__

#include "utils.hpp"
#include "utils.tpp"

namespace utils {

template bool add_to_ordered_vector<>(const long&, std::vector<long>&, std::less<long>);
template bool add_to_ordered_vector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

}
