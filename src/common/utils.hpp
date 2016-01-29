//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITP_MESH_UTILS_HPP__
#define __BITP_MESH_UTILS_HPP__

/*! \file */

#include <array>
#include <functional>
#include <vector>

namespace utils {

std::array<double, 3> cross_3D(std::array<double, 3> &x, std::array<double, 3> &y);
void normalize_3D(std::array<double, 3> &x);
void transpose_3D(std::array<std::array<double, 3>, 3> &A);

}

#endif
