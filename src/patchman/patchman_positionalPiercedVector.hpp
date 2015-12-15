//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_POSITIONAL_PIERCED_VECTOR_HPP__
#define __PATCHMAN_POSITIONAL_PIERCED_VECTOR_HPP__

#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

#include "positionalPiercedVector.tpp"

namespace pman {

extern template class PositionalPiercedVector<Cell>;
extern template class PositionalPiercedVector<Interface>;
extern template class PositionalPiercedVector<Vertex>;

}

#endif
