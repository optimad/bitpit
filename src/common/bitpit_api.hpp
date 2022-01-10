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

#ifndef __BITPIT_API_HPP__
#define __BITPIT_API_HPP__

/* API keyword definition for bitpit lib import/export */

#if defined(_MSC_VER)
    #if defined(BITPIT_DLLGLOBALDATA_EXPORT)
        #define BITPIT_API  __declspec(dllexport)
    #else 
        #define BITPIT_API  __declspec(dllimport)
    #endif 
#else 
    #define BITPIT_API  
#endif      


#endif
