//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_POSITIONAL_PIERCED_VECTOR_HPP__
#define __PATCHMAN_POSITIONAL_PIERCED_VECTOR_HPP__

#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

#include <bitpit_containers.hpp>

extern template class Positionalbitpit::PiercedVector<Cell>;
extern template class Positionalbitpit::PiercedVector<Interface>;
extern template class Positionalbitpit::PiercedVector<Vertex>;

#endif
