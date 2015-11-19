//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#define __PATCHMAN_UTILS_SRC__

#include "utils.hpp"
#include "utils.tpp"

namespace pman {

namespace utils {

template bool add_to_ordered_vector<>(const long&, std::vector<long>&, std::less<long>);

}

}
