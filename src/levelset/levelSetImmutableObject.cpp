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

# include "levelSetImmutableObject.hpp"
# include "levelSetKernel.hpp"

namespace bitpit {

/*!
 * \ingroup levelset
 * \class LevelSetImmutableObject
 * \brief The class LevelSetImmutableObject allows to define immutable
 * objects.
 *
 * Immutable objects initialize levelset information in the constructor and
 * those information cannot be updated later. Morover, immutable objects
 * cannot evaluate the levelset on arbitray points.
 */

/*!
 * Constructor.
 *
 * \param id is the id that will be assigned to object
 */
LevelSetImmutableObject::LevelSetImmutableObject(int id)
    : LevelSetObject(id)
{
}

/*!
 * Constructor.
 *
 * \param[in,out] source is an object whose contents will be acquired by the
 * immutable object. Since source contents will be moved into the object, on
 * output the source will be in an undefined state
 */
LevelSetImmutableObject::LevelSetImmutableObject(LevelSetCachedObject *source)
    : LevelSetObject(*source),
      LevelSetCachedObjectInterface(std::move(*source)),
      LevelSetSignedObjectInterface(std::move(*source))
{
}

/*!
 * Constructor.
 *
    * \param[in,out] source is an object whose contents will be acquired by the
 * immutable object. Since source contents will be moved into the object, on
 * output the source will be in an undefined state
 */
LevelSetImmutableObject::LevelSetImmutableObject(LevelSetBooleanObject *source)
    : LevelSetObject(*source)
{
    // Mesh information
    const VolumeKernel *mesh = (const_cast<const LevelSetBooleanObject *>(source))->getKernel()->getMesh();

    VolumeKernel::CellConstIterator cellBegin = mesh->cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh->cellConstEnd();

    // Fill narrow band information
    LevelSetNarrowBandCache *narrowBandCache = initializeNarrowBandCache();

    bool isSignStorageNeeded = false;
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        long cellId = cellItr.getId();
        bool cellInsideNarrowBand = source->isInNarrowBand(cellId);

        if (cellInsideNarrowBand) {
            double value = source->getValue(cellId);
            std::array<double,3> gradient = source->getGradient(cellId);

            LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->insert(cellId);
            narrowBandCache->set(narrowBandCacheItr, value, gradient);
        } else if (!isSignStorageNeeded) {
            short cellSign = source->getSign(cellId);
            isSignStorageNeeded = (cellSign != levelSetDefaults::SIGN);
        }
    }

    // Fill sign information
    if (isSignStorageNeeded) {
        LevelSetSignStorage *signStorage = initializeSignStorage();

        for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
            long cellId = cellItr.getId();
            std::size_t cellRawId = cellItr.getRawIndex();
            LevelSetSignStorage::KernelIterator signStorageItr = signStorage->rawFind(cellRawId);
            signStorage->at(signStorageItr) = source->getSign(cellId);
        }

        signStorage->setDirty(false);
    }
}

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetImmutableObject::getLevelSetInfo(long id)const{

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id);
    if(narrowBandCacheItr != narrowBandCache->end()){
        double value = narrowBandCache->getValue(narrowBandCacheItr);
        const std::array<double, 3> &gradient = narrowBandCache->getGradient(narrowBandCacheItr);

        return LevelSetInfo(value, gradient);
    }

    return LevelSetInfo();

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetImmutableObject::getLS(long id)const {

    return getValue(id);

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetImmutableObject::getValue(long id)const {

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id);
    if(narrowBandCacheItr != narrowBandCache->end()){
        return narrowBandCache->getValue(narrowBandCacheItr);
    }

    return getSign(id) * levelSetDefaults::VALUE;

}

/*!
 * Get the sign of the levelset function
 * @param[in] id cell id
 * @return sign of levelset
 */
short LevelSetImmutableObject::getSign(long id) const {

    // Check if the sign can be evaluated from narrowband value
    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id);
    if(narrowBandCacheItr != narrowBandCache->end()){
        double value = narrowBandCache->getValue(narrowBandCacheItr);

        return evalValueSign(value);
    }

    // Check if the sign can be evaluated from the storage
    if (!isSignStorageDirty()) {
        const LevelSetSignStorage *signStorage = getSignStorage();
        LevelSetSignStorage::KernelIterator signStorageItr = signStorage->find(id);
        LevelSetSignStorage::Sign storedSign = signStorage->at(signStorageItr);
        if (storedSign != LevelSetSignStorage::SIGN_UNDEFINED) {
            return static_cast<short>(storedSign);
        }
    }

    // Unable to evaluate the sign
    //
    // The sign cannot be evaluated, let's return the defualt sign.
    return levelSetDefaults::SIGN;

}

/*!
 * Get the levelset gradient of cell
 * @param[in] id cell id
 * @return levelset gradient in cell
 */
std::array<double,3> LevelSetImmutableObject::getGradient(long id) const {

    const LevelSetNarrowBandCache *narrowBandCache = getNarrowBandCache();
    LevelSetNarrowBandCache::KernelIterator narrowBandCacheItr = narrowBandCache->find(id);
    if(narrowBandCacheItr != narrowBandCache->end()){
        return narrowBandCache->getGradient(narrowBandCacheItr);
    }

    return levelSetDefaults::GRADIENT;

}

/*!
 * Clones the object
 *
 * \result A pointer to cloned object.
 */
LevelSetImmutableObject * LevelSetImmutableObject::clone() const
{
    return new LevelSetImmutableObject(*this);
}

/*!
 * Computes levelset information at the specified point.
 *
 * Immutable objects cannot evaluate the levelset on arbitray points.
 *
 * \param coords are the coordinates of the point
 * \return The levelset information at the specified point.
 */
LevelSetInfo LevelSetImmutableObject::computeLevelSetInfo(const std::array<double,3> &coords) const
{
    BITPIT_UNUSED(coords);

    throw std::runtime_error("Immutable objects cannot evaluate the levelset on arbitray points.");
}

/*!
 * Create the storage for the narrow band data.
 */
std::shared_ptr<LevelSetNarrowBandCache> LevelSetImmutableObject::createNarrowBandCache()
{
    return std::shared_ptr<LevelSetNarrowBandCache>(new LevelSetNarrowBandCache());
}

/*!
 * Create the storage for sign.
 */
std::shared_ptr<LevelSetSignStorage> LevelSetImmutableObject::createSignStorage()
{
    VolumeKernel *mesh = m_kernel->getMesh();
    assert(mesh);
    return std::shared_ptr<LevelSetSignStorage>(new LevelSetSignStorage(&(mesh->getCells())));
}

/*!
 * Clears all levelset information.
 */
void LevelSetImmutableObject::_clear()
{
    // Clear narrow band entries
    clearNarrowBandCache();

    // Clear sign propgation storage
    clearSignStorage();
}

/*!
 * Write the object to the specified stream.
 *
 * \param stream output stream
 */
void LevelSetImmutableObject::_dump(std::ostream &stream)
{
    // Narrow band storage
    dumpNarrowBandCache(stream);

    // Stored sign
    dumpSignStorage(stream);
}

/*!
 * Read the object from the specified stream.
 *
 * \param stream output stream
 */
void LevelSetImmutableObject::_restore(std::istream &stream)
{
    // Narrow band storage
    restoreNarrowBandCache(stream);

    // Stored sign
    restoreSignStorage(stream);
}

}
