/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#ifndef __BITPIT_COMMUNICATIONS_TAGS_HPP__
#define __BITPIT_COMMUNICATIONS_TAGS_HPP__

#include <memory>
#include <mpi.h>
#include <iostream>
#include <unordered_map>

#include <bitpit_IO.hpp>

namespace bitpit {

class CommunicationTags : public IndexGenerator<int> {

friend class PatchKernel;

public:
    static CommunicationTags & instance();

    int generate(MPI_Comm communicator);
    bool isAssigned(id_type id, MPI_Comm communicator);
    void setAssigned(int id, MPI_Comm communicator);
    void trash(int id, MPI_Comm communicator);

    id_type getLatest() = delete;
    id_type getHighest() = delete;

    void reset() = delete;

private:
    static std::unique_ptr<CommunicationTags> m_instance;

	CommunicationTags();

    CommunicationTags(CommunicationTags const &other) = delete;
    CommunicationTags& operator=(CommunicationTags const &other) = delete;

};

/*!
    \brief The namespace 'communications' contains routines for interacting
    with the communications.
*/
namespace communications {

    CommunicationTags & tags();

}

}

#endif
