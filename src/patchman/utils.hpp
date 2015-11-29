//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_UTILS_HPP__
#define __PATCHMAN_UTILS_HPP__

/*! \file */

#include <array>
#include <functional>
#include <vector>

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace pman {

namespace utils {

template <typename T, typename Comparator = std::less<T> >
bool add_to_ordered_vector(const T &value, std::vector<T> &list, Comparator comparator = Comparator());

#ifndef __PATCHMAN_UTILS_SRC__
extern template bool add_to_ordered_vector<>(const long&, std::vector<long>&, std::less<long>);
#endif

std::array<double, 3> cross_3D(std::array<double, 3> &x, std::array<double, 3> &y);
void normalize_3D(std::array<double, 3> &x);
void transpose_3D(std::array<std::array<double, 3>, 3> &A);


}

}

#endif
