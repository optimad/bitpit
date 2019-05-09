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

#include "communications_tags.hpp"

namespace bitpit {

/*!
	\class CommunicationTags
	\ingroup communications

	\brief The CommunicationTags oversee the handling of the communication
	tags.
*/

/*
    Initialize communication tags instance.
*/
std::unique_ptr<CommunicationTags> CommunicationTags::m_instance = nullptr;

/*!
    Returns an instance of the communication tag manager.

    \result An instance of the communication tag manager.
*/
CommunicationTags & CommunicationTags::instance()
{
    if (!m_instance) {
        m_instance = std::unique_ptr<CommunicationTags>(new CommunicationTags());
    }

    return *m_instance;
}

/*!
	Creates a new communication tag manager.
*/
CommunicationTags::CommunicationTags()
{
}

/*!
    Generates a unique tag.

    If the trash is empty a new tag is generated, otherwise a tag taken
    from the trash is recycled.

    The tag is shared among the processes that belong to the specified
    communicator.

    \param communicator is the MPI communicator
    \return A new unique tag.
*/

int CommunicationTags::generate(MPI_Comm communicator)
{
    // Early return if the communicator is not valid
    if (communicator == MPI_COMM_NULL) {
        return IndexGenerator<int>::generate();
    }

    // The master generates a unique tag
    int tag;

    int rank;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) {
        tag = IndexGenerator<int>::generate();
    }

    // Share the tag among the processes
    MPI_Bcast(&tag, 1, MPI_INT, 0, communicator);

    return tag;
}

/*!
    Checks if an tag is currently assigned.

    \param communicator is the MPI communicator
    \return True if the tag has already been assigned, false otherwise.
*/
bool CommunicationTags::isAssigned(int tag, MPI_Comm communicator)
{
    // The master checks if the tag is assigned
    bool assigned;

    int rank;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) {
        assigned = IndexGenerator<int>::isAssigned(tag);
    }

    // Share the result among the processes
    MPI_Bcast(&assigned, 1, MPI_C_BOOL, 0, communicator);

    return assigned;
}

/*!
    Marks the specified tag as currently assigned.

    \param tag is the tag to be marked as assigned
    \param communicator is the MPI communicator
*/
void CommunicationTags::setAssigned(int tag, MPI_Comm communicator)
{
    // Only the master needs to set the tag as assigned
    int rank;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) {
        IndexGenerator<int>::setAssigned(tag);
    }
}

/*!
    Trashes a tag.

    A trashed tag is a tag no more used that can be recycled.

    \param tag is the tag that will be trashed
    \param communicator is the MPI communicator
*/
void CommunicationTags::trash(int tag, MPI_Comm communicator)
{
    // Wait all the processes
    //
    // We need to avoid that a tag is trashed while some processes are still
    // using it.
    MPI_Barrier(communicator);

    // Only the master needs to trash the tag
    int rank;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) {
        IndexGenerator<int>::trash(tag);
    }
}

// Communication manager global functions
namespace communications {

    /*!
        Returns the communication tag manager.

        \result The communication tag manager.
    */
    CommunicationTags & tags()
    {
        return CommunicationTags::instance();
    }

}

}
