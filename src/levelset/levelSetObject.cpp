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

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# endif

# include "levelSet.hpp"

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@interface LevelSetObject
	@brief Interface class for all objects with respect to whom the levelset function may be computed.
*/

/*!
 * Destructor
 */
LevelSetObject::~LevelSetObject( ){
};

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetObject::LevelSetObject( int id) : m_id(id){
};

/*!
 * Get the id 
 * @return id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
};

/*!
 * Clears data structure after mesh modification
 * @param[in] mapper mapper describing mesh modifications
 */
void LevelSetObject::clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){
    BITPIT_UNUSED(mapper) ;
    return;
}

/*!
 * Clears data structure outside narrow band
 * @param[in] visitee LevelSetKernel with narrow band information
 */
void LevelSetObject::filterOutsideNarrowBand( LevelSetKernel *visitee ){
    BITPIT_UNUSED(visitee) ;
    return;
}

/*!
 * Gets the number of support items within the narrow band of cell
 * @param[in] id index of cell
 */
int LevelSetObject::getSupportCount( const long &id )const{
    BITPIT_UNUSED( id) ;
    return 0;
}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
long LevelSetObject::getClosestSupport( const long &id )const{
    BITPIT_UNUSED( id) ;
    return levelSetDefaults::SUPPORT;
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::fstream &stream ){
    bitpit::genericIO::flushBINARY(stream, m_id) ;
    dumpDerived(stream) ;
    return;
};

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::fstream &stream ){
    bitpit::genericIO::absorbBINARY(stream, m_id) ;
    restoreDerived(stream) ;
    return;
};

}
