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

# ifndef __BITPIT_LEVELSET_OBJECT_FACTORY_TPP__
# define __BITPIT_LEVELSET_OBJECT_FACTORY_TPP__

#include "levelSetUnstructuredKernel.hpp"
namespace bitpit {

/*!
	@class      LevelSetObjectFactory
	@ingroup    levelset
	@brief      Allows to create levelset ojbects.
*/

/*!
 * Create a new boolean object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createBooleanObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args)
{
    BITPIT_UNUSED(kernel);
    BITPIT_UNUSED(storageType);

    return std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject(std::forward<Args>(args)...));
}

/*!
 * Create a new immutable object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param source is the source for the immutable object
 */
template<template<typename> class narrow_band_cache_t>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createImmutableObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, LevelSetObject *source)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createImmutableObject<LevelSetCartesianKernel, narrow_band_cache_t>(cartesianKernel, storageType, source);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createImmutableObject<LevelSetOctreeKernel, narrow_band_cache_t>(octreeKernel, storageType, source);
    } else if ( const LevelSetUnstructuredKernel *unstructuredKernel = dynamic_cast<const LevelSetUnstructuredKernel *>(kernel) ){
        return _createImmutableObject<LevelSetUnstructuredKernel, narrow_band_cache_t>(unstructuredKernel, storageType, source);
    }


    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new immutable object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param source is the source for the immutable object
 */
template<typename kernel_t, template<typename> class narrow_band_cache_t>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createImmutableObject(const kernel_t *kernel, LevelSetStorageType storageType, LevelSetObject *source)
{
    BITPIT_UNUSED(kernel);

    switch (storageType) {

    case LevelSetStorageType::SPARSE:
        if (LevelSetBooleanObject *booleanSource = dynamic_cast<LevelSetBooleanObject *>(source)) {
            return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>>(booleanSource));
        } else if (LevelSetCachedObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>> *cachedSource = dynamic_cast<LevelSetCachedObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>> *>(source)) {
            return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>>(cachedSource));
        }
        break;

    case LevelSetStorageType::DENSE:
        if (LevelSetBooleanObject *booleanSource = dynamic_cast<LevelSetBooleanObject *>(source)) {
            return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject<narrow_band_cache_t<typename kernel_t::DenseStorageManager>>(booleanSource));
        } else if (LevelSetCachedObject<narrow_band_cache_t<typename kernel_t::DenseStorageManager>> *cachedSource = dynamic_cast<LevelSetCachedObject<narrow_band_cache_t<typename kernel_t::DenseStorageManager>> *>(source)) {
            return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject<narrow_band_cache_t<typename kernel_t::DenseStorageManager>>(cachedSource));
        }
        break;

    default:
        BITPIT_UNREACHABLE("Storage type not supported");

    }

    BITPIT_UNREACHABLE("Unable to create an immutable object from the specified object.");
}

/*!
 * Create a new mask object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param args are the arguments that will be forwarded to the constructor
 */
template<template<typename> class narrow_band_cache_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createMaskObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createMaskObject<LevelSetCartesianKernel, narrow_band_cache_t, Args...>(cartesianKernel, storageType, std::forward<Args>(args)...);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createMaskObject<LevelSetOctreeKernel, narrow_band_cache_t, Args...>(octreeKernel, storageType, std::forward<Args>(args)...);
    } else if ( const LevelSetUnstructuredKernel *unstructuredKernel = dynamic_cast<const LevelSetUnstructuredKernel *>(kernel) ){
        return _createMaskObject<LevelSetUnstructuredKernel, narrow_band_cache_t, Args...>(unstructuredKernel, storageType, std::forward<Args>(args)...);
    }

    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new mask object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename kernel_t, template<typename> class narrow_band_cache_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createMaskObject(const kernel_t *kernel, LevelSetStorageType storageType, Args&&... args)
{
    BITPIT_UNUSED(kernel);

    switch (storageType) {

    case LevelSetStorageType::SPARSE:
        return std::unique_ptr<LevelSetObject>(new LevelSetMaskObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>>(std::forward<Args>(args)...));

    case LevelSetStorageType::DENSE:
        return std::unique_ptr<LevelSetObject>(new LevelSetMaskObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>>(std::forward<Args>(args)...));

    default:
        BITPIT_UNREACHABLE("Storage type not supported");

    }
}

/*!
 * Create a new segmentation object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param args are the arguments that will be forwarded to the constructor
 */
template<template<typename> class narrow_band_cache_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createSegmentationObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createSegmentationObject<LevelSetCartesianKernel, narrow_band_cache_t, Args...>(cartesianKernel, storageType, std::forward<Args>(args)...);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createSegmentationObject<LevelSetOctreeKernel, narrow_band_cache_t, Args...>(octreeKernel, storageType, std::forward<Args>(args)...);
    } else if ( const LevelSetUnstructuredKernel *unstructuredKernel = dynamic_cast<const LevelSetUnstructuredKernel *>(kernel) ){
        return _createSegmentationObject<LevelSetUnstructuredKernel, narrow_band_cache_t, Args...>(unstructuredKernel, storageType, std::forward<Args>(args)...);
    }

    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new segmentation object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param storageType is the storage type
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename kernel_t, template<typename> class narrow_band_cache_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createSegmentationObject(const kernel_t *kernel, LevelSetStorageType storageType, Args&&... args)
{
    BITPIT_UNUSED(kernel);

    switch (storageType) {

    case LevelSetStorageType::SPARSE:
        return std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject<narrow_band_cache_t<typename kernel_t::SparseStorageManager>>(std::forward<Args>(args)...));

    case LevelSetStorageType::DENSE:
        return std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject<narrow_band_cache_t<typename kernel_t::DenseStorageManager>>(std::forward<Args>(args)...));

    default:
        BITPIT_UNREACHABLE("Storage type not supported");

    }
}

}

# endif
