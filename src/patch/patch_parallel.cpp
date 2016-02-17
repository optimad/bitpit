/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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
#if ENABLE_MPI==1

#include <mpi.h>

#include "patch.hpp"

namespace bitpit {

/*!
	\ingroup patch
	@{

	\class Patch
*/

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void Patch::setCommunicator(MPI_Comm *communicator)
{
	// Free previous communicator
	if (m_communicator) {
		MPI_Comm_free(m_communicator);
		m_communicator = nullptr;
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	if (communicator) {
		MPI_Comm_dup(*communicator, m_communicator);
	}

	// Get MPI information
	if (m_communicator) {
		MPI_Comm_size(*m_communicator, &m_nProcessors);
		MPI_Comm_rank(*m_communicator, &m_rank);
	} else {
		m_rank        = 0;
		m_nProcessors = 1;
	}

	// Set parallel data for the VTK output
	setParallel(m_nProcessors, m_rank);
}

/*!
	Unsets the MPI communicator to be used for parallel communications.
*/
void Patch::unsetCommunicator()
{
	setCommunicator(nullptr);
}

/*!
	Gets the MPI communicator associated to the patch

	\return The MPI communicator associated to the patch.
*/
MPI_Comm & Patch::getCommunicator() const
{
	return *m_communicator;
}

/*!
	Gets the MPI rank associated to the patch

	\return The MPI rank associated to the patch.
*/
int Patch::getRank() const
{
	return m_rank;
}

/*!
	Gets the MPI processors in the communicator associated to the patch

	\return The MPI processors in the communicator associated to the patch
*/
int Patch::getProcessorCount() const
{
	return m_nProcessors;
}

/*!
	@}
*/

}

#endif
