/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#include "patch_kernel.hpp"
#include "patch_manager.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	@{
*/

/*!
	\class PatchManager

	\brief The PatchManager oversee the handling of the patches.
*/

int const PatchManager::NULL_PATCH_ID = IndexGenerator::NULL_ID;

/*
    Initialize logger manager instance.
*/
std::unique_ptr<PatchManager> PatchManager::m_manager = nullptr;

/*!
    Returns an instance of the patch manager.

    \result An instance of the patch manager.
*/
PatchManager & PatchManager::manager()
{
    if (!m_manager) {
        m_manager = std::unique_ptr<PatchManager>(new PatchManager());
    }

    return *m_manager;
}

/*!
	Registers a patch in the manager

	\param patch is a pointer to the patch to be registered
	\param id is the id that will be assigned to the patch
	\result The id assigned to the patch.
*/
int PatchManager::registerPatch(PatchKernel *patch, int id)
{
	if (id >= 0) {
		if (m_idGenerator.isIdAssigned(id)) {
			throw std::runtime_error ("A patch with the same id already exists");
		}
	} else {
		id = m_idGenerator.generateId();
	}

	patch->setId(id);
	m_patchIds[patch] = id;
	m_patchOrder.push_back(patch);

	return id;
}

/*!
	Un-registers a patch in the manager

	\param patch is a pointer to the patch to be un-registered
*/
void PatchManager::unregisterPatch(PatchKernel *patch)
{
	auto iterator = m_patchIds.find(patch);
	if (iterator == m_patchIds.end()) {
		throw std::runtime_error ("The patch to be unregistered does not exist");
	}

	int id = iterator->second;
	m_idGenerator.trashId(id);

	m_patchIds.erase(iterator);
	for (auto itr = m_patchOrder.begin(); itr != m_patchOrder.end(); ++itr) {
		if (*itr == patch) {
			m_patchOrder.erase(itr);
			break;
		}
	}
}

/*!
	Creates a new patch manager.
*/
PatchManager::PatchManager()
{
}

/*!
 *  Write the patch manager data to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PatchManager::dump(std::ostream &stream)
{
	m_idGenerator.dump(stream);
}

/*!
 *  Restore the patch manager data from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void PatchManager::restore(std::istream &stream)
{
	m_idGenerator.restore(stream);
}

/*!
 *  Write the registered patches and the patch manager data to the specified
 *  stream.
 *
 *  \param stream is the stream to write to
 */
void PatchManager::dumpAll(std::ostream &stream)
{
	for (PatchKernel *patch : m_patchOrder) {
		patch->dump(stream);
	}

	dump(stream);
}

/*!
 *  Restore the registered patches and the patch manager data from the
 *  specified stream.
 *
 *  \param stream is the stream to read from
 */
void PatchManager::restoreAll(std::istream &stream)
{
	m_idGenerator.reset();

	for (PatchKernel *patch : m_patchOrder) {
		patch->restore(stream);
	}

	restore(stream);
}

/*!
    \ingroup patchkernel
    @{
*/

// Patch manager global functions
namespace patch {

    // Generic global functions

    /*!
        Returns the logger manager.

        \result The logger manager.
    */
    PatchManager & manager()
    {
        return PatchManager::manager();
    }

}

/*!
    @}
*/

}
