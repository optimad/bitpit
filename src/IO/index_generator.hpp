/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BITPIT_INDEX_GENERATOR_HPP__
#define __BITPIT_INDEX_GENERATOR_HPP__

#include <cassert>
#include <deque>
#include <iostream>
#include <limits>

#include "bitpit_common.hpp"

namespace bitpit {

template<typename id_t = long>
class IndexGenerator {

static_assert(std::is_integral<id_t>::value, "Index has to be an integer!");

public:
    typedef id_t id_type;

    static const id_type NULL_ID;

    IndexGenerator();

    id_type generate();
    bool isAssigned(id_type id);
    void setAssigned(id_type id);
    void trash(id_type id);

    id_type getLatest();
    id_type getHighest();

    void reset();

    void dump(std::ostream &stream) const;
    void restore(std::istream &stream);

private:
    id_type m_latest;
    id_type m_highest;
    std::deque<id_type> m_trash;

    int getBinaryArchiveVersion() const;

};

#ifndef __BITPIT_COMMON_UTILS_SRC__
extern template class IndexGenerator<int>;
extern template class IndexGenerator<long>;
#endif

}

// Include the implementation
#include "index_generator.tpp"

#endif
