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

# ifndef __BITPIT_LEVELSET_OBJECT_FACTORY_HPP__
# define __BITPIT_LEVELSET_OBJECT_FACTORY_HPP__

#include "levelSetBooleanObject.hpp"
#include "levelSetCartesianKernel.hpp"
#include "levelSetImmutableObject.hpp"
#include "levelSetMaskObject.hpp"
#include "levelSetOctreeKernel.hpp"
#include "levelSetSegmentationObject.hpp"

namespace bitpit {

class LevelSetObjectFactory {

public:
    template<typename... Args>
    static std::unique_ptr<LevelSetObject> createBooleanObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args);

    template<template<typename> class narrow_band_cache_t>
    static std::unique_ptr<LevelSetObject> createImmutableObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, LevelSetObject *source);

    template<template<typename> class narrow_band_cache_t, typename... Args>
    static std::unique_ptr<LevelSetObject> createMaskObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args);

    template<template<typename> class narrow_band_cache_t, typename... Args>
    static std::unique_ptr<LevelSetObject> createSegmentationObject(const LevelSetKernel *kernel, LevelSetStorageType storageType, Args&&... args);

private:
    template<typename kernel_t, template<typename> class narrow_band_cache_t>
    static std::unique_ptr<LevelSetObject> _createImmutableObject(const kernel_t *kernel, LevelSetStorageType storageType, LevelSetObject *source);

    template<typename kernel_t, template<typename> class narrow_band_cache_t, typename... Args>
    static std::unique_ptr<LevelSetObject> _createMaskObject(const kernel_t *kernel, LevelSetStorageType storageType, Args&&... args);

    template<typename kernel_t, template<typename> class narrow_band_cache_t, typename... Args>
    static std::unique_ptr<LevelSetObject> _createSegmentationObject(const kernel_t *kernel, LevelSetStorageType storageType, Args&&... args);

};

}

// Include template implementations
#include "levelSetObjectFactory.tpp"

#endif
