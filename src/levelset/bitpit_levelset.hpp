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

#ifndef __BITPIT_MODULE_LEVELSET_HPP__
#define __BITPIT_MODULE_LEVELSET_HPP__
#include "moduleBegin.hpp"

/*!
 * @defgroup levelset LevelSet
 * @{
 *     @defgroup levelsetEnums LevelSet Enumerations
 * @}
 */

#include "levelSetCommon.hpp"

#include "levelSetKernel.hpp"
#include "levelSetCartesianKernel.hpp"
#include "levelSetOctreeKernel.hpp"
#include "levelSetUnstructuredKernel.hpp"

#include "levelSetObject.hpp"
#include "levelSetProxyObject.hpp"
#include "levelSetCachedObject.hpp"
#include "levelSetSegmentationObject.hpp"
#include "levelSetSignedObject.hpp"
#include "levelSetBooleanObject.hpp"
#include "levelSetMaskObject.hpp"

#include "levelSet.hpp"

#include "moduleEnd.hpp"
#endif
