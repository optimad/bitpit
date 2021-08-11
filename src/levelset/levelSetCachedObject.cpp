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

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"
# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSignPropagator.hpp"
# include "levelSet.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@interface LevelSetCachedObject
	@brief Interface class for all objects which need to store the discrete values of levelset function.
*/

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetCachedObject::LevelSetCachedObject(int id)
    : LevelSetObject(id), LevelSetSignStorage()
{
}

/*!
 * Get LevelSetInfo of cell
 * @param[in] i cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetCachedObject::getLevelSetInfo( long i)const{

    auto itr = m_ls.find(i);
    if ( itr == m_ls.end() ){
        return LevelSetInfo();
    }

    return *itr;

} 

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetCachedObject::getLS( long id)const {

    return getValue(id);

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetCachedObject::getValue( long id)const {

    auto itr = m_ls.find(id);
    if ( itr == m_ls.end() ){
        return getSign(id) * levelSetDefaults::VALUE;
    }

    return itr->value;

}

/*!
 * Get the sign of the levelset function
 * @param[in] id cell id
 * @return sign of levelset
 */
short LevelSetCachedObject::getSign( long id ) const {

    // Check if the sign can be evaluated from narrowband distance
    auto itr = m_ls.find(id);
    if ( itr != m_ls.end() ){
        return evalValueSign(itr->value);
    }

    // Check if the sign can be evaluated from the propagation
    if (!isStoredSignDirty()) {
        LevelSetSignStorage::Sign propagatedSign = getStoredSign(id);
        if (propagatedSign != LevelSetSignStorage::SIGN_UNDEFINED) {
            return static_cast<short>(propagatedSign);
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
std::array<double,3> LevelSetCachedObject::getGradient(long id) const {

    auto itr = m_ls.find(id);
    if ( itr == m_ls.end() ){
        return levelSetDefaults::GRADIENT;
    }

    return itr->gradient;

}

/*! 
 * Deletes non-existing items after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetCachedObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){

    for ( auto & map : mapper ){
        if( map.entity == adaption::Entity::ENTITY_CELL ){
            if( map.type == adaption::Type::TYPE_DELETION || 
                map.type == adaption::Type::TYPE_PARTITION_SEND  ||
                map.type == adaption::Type::TYPE_REFINEMENT  ||
                map.type == adaption::Type::TYPE_COARSENING  ){

                for ( auto & parent : map.previous){
                    if( isInNarrowBand(parent) )
                        m_ls.erase(parent,true) ;
                }
            }
        }
    }

    m_ls.flush() ;

    // Sign propgation needs to be updated
    setStoredSignDirty(true);

}

/*! 
 * Clears all levelset information
 */
void LevelSetCachedObject::_clear( ){
    m_ls.clear() ;

    clearSignStorage();
}

/*!
 * If cell lies within the narrow band and hence its levelset is computed
 * exactly.
 * @param[in] id cell id
 * @return true/false if the centroid is in narrow band
 */
bool LevelSetCachedObject::isInNarrowBand(long id)const{
    return m_ls.exists(id);
}

/*!
 * Writes LevelSetCachedObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_dump( std::ostream &stream ){

    // Levelset information
    utils::binary::write(stream, m_ls.size() ) ;
    bitpit::PiercedVector<LevelSetInfo>::iterator   infoItr, infoEnd = m_ls.end() ;

    for( infoItr=m_ls.begin(); infoItr!=infoEnd; ++infoItr){
        utils::binary::write(stream, infoItr.getId()) ;
        utils::binary::write(stream, infoItr->value) ;
        utils::binary::write(stream, infoItr->gradient) ;
    }

    // Stored sign
    LevelSetSignStorage::dumpStoredSign( stream );
}

/*!
 * Reads LevelSetCachedObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetCachedObject::_restore( std::istream &stream ){

    // Levelset information
    std::size_t nInfoItems;
    utils::binary::read(stream, nInfoItems);
    m_ls.reserve(nInfoItems);

    for( std::size_t i=0; i<nInfoItems; ++i){
        long id;
        utils::binary::read(stream, id) ;

        double value;
        utils::binary::read(stream, value) ;

        std::array<double,3> gradient;
        utils::binary::read(stream, gradient) ;

        m_ls.insert(id, LevelSetInfo(value, gradient)) ;
    }

    // Stored sign
    LevelSetSignStorage::restoreStoredSign( stream );
}

#if BITPIT_ENABLE_MPI

/*!
 * Flushing of data to communication buffers for partitioning
 * Sign data is not written into the buffer, because sign storage is kept in
 * sync with the mesh, hence, when this function is called, entries associated
 * with the cells to send as already been deleted.
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetCachedObject::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    // Evaluate the size of the buffer
    std::size_t nInfoItems = 0;
    for( long id : sendList){
        if( m_ls.exists(id)){
            ++nInfoItems ;
        }
    }

    dataBuffer.setSize(dataBuffer.getSize() + sizeof(std::size_t) + nInfoItems* (sizeof(std::size_t) +4*sizeof(double)) ) ;

    // Fill the buffer
    dataBuffer << nInfoItems ;

    for( std::size_t k = 0; k < sendList.size(); ++k){
        long id = sendList[k];
        auto infoItr = m_ls.find(id) ;
        if( infoItr != m_ls.end() ){
            dataBuffer << k ;
            dataBuffer << infoItr->value ;
            dataBuffer << infoItr->gradient ;
        }
    }
}

/*!
 * Processing of communication buffer into data structure
 * Sign data is not read from the buffer, because sign storage is kept in
 * sync with the mesh, hence, when the buffer is written, entries associated
 * with the cells to send as already been deleted.
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetCachedObject::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    std::size_t nInfoItems ;
    dataBuffer >> nInfoItems ;

    for( std::size_t i=0; i<nInfoItems; ++i){
        // Determine the id of the element
        std::size_t k ;
        dataBuffer >> k ;
        long id = recvList[k] ;

        // Assign the data of the element
        auto infoItr = m_ls.find(id) ;
        if( infoItr == m_ls.end() ){
            infoItr = m_ls.emplace(id) ;
        }

        dataBuffer >> infoItr->value ;
        dataBuffer >> infoItr->gradient ;
    }
}

#endif 

}
