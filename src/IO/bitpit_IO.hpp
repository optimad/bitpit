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

/*!
 * @defgroup IO Input/Output
 * @{
 * @defgroup GenericIO Generic
 * @defgroup DuneGridFormat DuneGridFormat (DGF)
 * @defgroup STereoLithography STereoLithography (STL)
 * @defgroup VisualizationToolKit VisualizationToolKit (VTK)
 * @{
 * @defgroup VTKEnums
 * @}
 * @}
 *
 */

#ifndef __BITPIT_IO__
#define __BITPIT_IO__

#include <bitpit_version.hpp>

#include "FileHandler.hpp"
#include "VTK.hpp"
#ifdef IO_ENABLE_VTK_WRAPPERS
#include "VTKWrappers.hpp"
#endif
#include "DGF.hpp"
#include "GenericIO.hpp"
#include "STL.hpp"

#endif
