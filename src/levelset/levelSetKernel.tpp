/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2022 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_KERNEL_TPP__
# define __BITPIT_LEVELSET_KERNEL_TPP__

namespace bitpit {

/*!
    @interface  LevelSetCachedKernel
    @ingroup levelset
    @brief  Mesh specific implementation to calculate the levelset function.
*/

/*!
 * Default constructor.
 */
template<typename CellCacheEntry>
LevelSetCachedKernel<CellCacheEntry>::LevelSetCachedKernel() : LevelSetKernel() {
}

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
template<typename CellCacheEntry>
LevelSetCachedKernel<CellCacheEntry>::LevelSetCachedKernel( VolumeKernel *patch): LevelSetKernel(patch) {
}

/*!
 * Updates the kernel after an adaption.
 *
 * @param[in] adaptionData are the information about the adaption
 */
template<typename CellCacheEntry>
void LevelSetCachedKernel<CellCacheEntry>::update( const std::vector<adaption::Info> &adaptionData ) {

    // Update base class
    LevelSetKernel::update( adaptionData );

    // Prune cell cache
    pruneCellCache( adaptionData );

}

/*!
 * Remove stale cell entries from the cache after an adaption.
 *
 * @param[in] adaptionData are the information about the adaption
 */
template<typename CellCacheEntry>
void LevelSetCachedKernel<CellCacheEntry>::pruneCellCache( const std::vector<adaption::Info> &adaptionData ) {

    // If there are no cells in the mesh we can just delete all the cache
    if ( m_mesh->getCellCount() == 0) {
        std::unordered_map<long, CellCacheEntry>().swap( m_cellCache ) ;
        return;
    }

    // Remove the previous cells from the cache
    for ( const adaption::Info &adaptionInfo : adaptionData ){
        if( adaptionInfo.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for ( long previousId : adaptionInfo.previous){
            auto cellCacheItr = m_cellCache.find( previousId ) ;
            if ( cellCacheItr == m_cellCache.end() ) {
                continue ;
            }

            m_cellCache.erase( cellCacheItr ) ;
        }
    }

}

/*!
 * Compute the geometrical information associated with the specified cell.
 *
 * If cell information have been already evaluated for the specified cell, the cached
 * value is returned. Otherwise cell information are evaluated and stored in the cache.
 *
 * @param[in] id is the id of cell
 * @return The geometrical information associated with the specified cell.
 */
template<typename CellCacheEntry>
const CellCacheEntry & LevelSetCachedKernel<CellCacheEntry>::computeCellCacheEntry( long id ) const {

    auto cellInfoItr = m_cellCache.find( id ) ;
    if ( cellInfoItr == m_cellCache.end() ) {
        cellInfoItr =  m_cellCache.insert( { id, CellCacheEntry( *m_mesh, id ) } ).first ;
    }

    return cellInfoItr->second ;

}

}

# endif
