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
	@ingroup    levelset
	@class      LevelSetObject
	@brief      LevelSet object.
*/

/*!
 * Destructor
 */
LevelSetObject::~LevelSetObject( ){
};

/*!
 * Constructor
 * @param[in] id Id assigned to object
 */
LevelSetObject::LevelSetObject( int id) : m_id(id){
};

/*!
 * Get the id 
 * @return Id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
};

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

# if BITPIT_ENABLE_MPI
void LevelSetObject::finalizeMPI(){
}
# endif

}
