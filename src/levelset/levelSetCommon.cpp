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


# include "levelSetCommon.hpp"

namespace bitpit {

/*!
 * @ingroup levelset
 * @class  LevelSetInfo
 *
 * @brief  A public container which includes all information provided by LevelSet
 *
 * LevelSetInfo conatins the following information
 * - distance to closest object
 * - gradient of level set function
 * - the id of the closest object
 * - the patch of closest object which contains the projection point
 * - the id of the surface element closest to the projection point 
 *
 * If a grid point lies within the narrow band of an object, all of these information are available, provided that the clsest object may provide them.
 * On the contrary, if the grid point lies were the levelset has been propagated not all of these information may be available; 
 * in particular the object, patch and segment may be set to the default values.
 *
*/

/*!
 * Default constructor
 */
LevelSetInfo::LevelSetInfo() :value(levelSetDefaults::VALUE), gradient(levelSetDefaults::GRADIENT) {
}

/*!
 * Complete constructor
 * @param[in] v value of levelset function
 * @param[in] g gradient of levelset function
 */
LevelSetInfo::LevelSetInfo(const double &v, const std::array<double,3> &g) :value(v), gradient(g) {
}

}

