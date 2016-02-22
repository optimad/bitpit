//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PIERCED_VECTOR_HPP__
#define __PATCHMAN_PIERCED_VECTOR_HPP__

#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

#include <bitpit_containers.hpp>

extern template class bitpit::PiercedVector<Cell>;
extern template class bitpit::PiercedVector<Interface>;
extern template class bitpit::PiercedVector<Vertex>;

#endif
