//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PIERCED_VECTOR_HPP__
#define __PATCHMAN_PIERCED_VECTOR_HPP__

#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

#include "piercedVector.tpp"

extern template class PiercedVector<Cell>;
extern template class PiercedVector<Interface>;
extern template class PiercedVector<Vertex>;

#endif
