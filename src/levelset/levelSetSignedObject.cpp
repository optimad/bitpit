/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif

# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetBoundedObject.hpp"
# include "levelSetKernel.hpp"
# include "levelSetObject.hpp"
# include "levelSetSignedObject.hpp"

namespace bitpit {

/*!
 * \ingroup levelset
 * \class LevelSetSignStorage
 * \brief The class LevelSetSignStorage allows to store the levelset sign
 * on the whole mesh.
 */

const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_UNDEFINED = -2;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_NEGATIVE  = -1;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_ZERO      =  0;
const LevelSetSignStorage::Sign LevelSetSignStorage::SIGN_POSITIVE  =  1;

/*!
 * Constructor.
 *
 * \param kernel is the kernel
 */
LevelSetSignStorage::LevelSetSignStorage(PiercedKernel<long> *kernel) : LevelSetExternalPiercedStorageManager(kernel) {

    m_signs = addStorage<Sign>(getStorageCount(), 1, PiercedSyncMaster::SYNC_MODE_JOURNALED);

}

/*!
 * Get the sign of the specified cell.
 *
 * \param itr is an iterator pointing to the cell
 * \result The sign of the specified cell.
 */
LevelSetSignStorage::Sign LevelSetSignStorage::at(const KernelIterator &itr) const
{
    std::size_t rawId = itr.getRawIndex();

    return m_signs->rawAt(rawId);
}

/*!
 * Get a reference to the sign of the specified cell.
 *
 * \param itr is an iterator pointing to the cell
 * \result A reference to the sign of the specified cell.
 */
LevelSetSignStorage::Sign & LevelSetSignStorage::at(const KernelIterator &itr)
{
    std::size_t rawId = itr.getRawIndex();

    return m_signs->rawAt(rawId);
}

/*!
 * Set the sign of all the stored cells.
 *
 * \param sign is the sign that will be set
 */
void LevelSetSignStorage::fill(Sign sign)
{
    m_signs->fill(sign);
}

/*!
 * Exchanges the content of the storage with the content the specified other
 * storage.
 *
 * \param other is another storage whose content is swapped with that of this
 * storage
 */
void LevelSetSignStorage::swap(LevelSetSignStorage &other) noexcept
{
    LevelSetExternalPiercedStorageManager::swap(other);

    std::swap(other.m_signs, m_signs);
}

/*!
 * \ingroup levelset
 * \class LevelSetSignedObjectInterface
 * \brief The class LevelSetSignedObjectInterface allows to define objects
 * that can store the signe on the whole domain.
 */

/*!
 * Initialize the storage.
 */
LevelSetSignStorage * LevelSetSignedObjectInterface::initializeSignStorage()
{
    m_signStorage = createSignStorage();

    return getSignStorage();
}

/*!
 * Get a pointer to the storage.
 *
 * \result A pointer to the storage.
 */
LevelSetSignStorage * LevelSetSignedObjectInterface::getSignStorage()
{
    return m_signStorage.get();
}

/*!
 * Get a constant pointer to the storage.
 *
 * \result A constant pointer to the storage.
 */
const LevelSetSignStorage * LevelSetSignedObjectInterface::getSignStorage() const
{
    return m_signStorage.get();
}

/*!
 * Check if the sign is dirty.
 *
 * The sign can be dirty if it was not yet initialized or if it is not
 * up-to-date.
 *
 * \result Return true if the sign is dirty, false otherwise.
 */
bool LevelSetSignedObjectInterface::isSignStorageDirty() const
{
    if (!m_signStorage) {
        return true;
    }

    return m_signStorage->isDirty();
}

/*!
 * Set the sign as dirty.
 *
 * The sign can be set as non-dirty only it it has been already initialized.
 *
 * \param dirty if set to true the sign will be set as dirty, otherwise, if
 * the sign has been initialized, it will be set as non-dirty
 */
void LevelSetSignedObjectInterface::setSignStorageDirty(bool dirty)
{
    if (!m_signStorage) {
        return;
    }

    m_signStorage->setDirty(dirty);
}

/*!
 * Clear the storage.
 */
void LevelSetSignedObjectInterface::clearSignStorage()
{
    // Stored sign
    if (m_signStorage) {
        m_signStorage->clear();
    }
}

/*!
 * Dump the storage.
 *
 * \param stream is the output stream
 */
void LevelSetSignedObjectInterface::dumpSignStorage(std::ostream &stream)
{
    // Stored sign
    if (m_signStorage) {
        m_signStorage->dump( stream );
    }
}

/*!
 * Restore the storage.
 *
 * \param stream is the input stream
 */
void LevelSetSignedObjectInterface::restoreSignStorage(std::istream &stream)
{
    // Stored sign
    if (m_signStorage) {
        m_signStorage->restore( stream );
    }
}

/*!
 * Exchanges the content of the object with the content the specified other
 * object.
 *
 * \param other is another object whose content is swapped with that of this
 * object
 */
void LevelSetSignedObjectInterface::swap(LevelSetSignedObjectInterface &other) noexcept
{
    m_signStorage.swap(other.m_signStorage);
}

}
