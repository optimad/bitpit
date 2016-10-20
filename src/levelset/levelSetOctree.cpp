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

# include "bitpit_voloctree.hpp"
# include "levelSetKernel.hpp"
# include "levelSetOctree.hpp"

namespace bitpit {

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctree
 *  @brief      Implements LevelSetKernel for octree meshes
 */

/*!
 * Destructor
 */
LevelSetOctree::~LevelSetOctree( ){
    m_octree = NULL ;
};

/*!
 * Constructor
 */
LevelSetOctree::LevelSetOctree(VolOctree & patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_octree = &patch ;
};

VolOctree* LevelSetOctree::getOctreeMesh() const{
    return m_octree ;
}
/*!
 * Compute size of narrow band given a cell.
 * This function assumes that the octree is balanced 2:1
 * @param[in] id is the id of the cell
 */
double LevelSetOctree::computeRSearchFromCell( long id ){
    int level = m_octree->getCellLevel(id) ;
    return  computeRSearchFromLevel( level ) ;
};

/*!
 * Compute size of narrow band given the coarsest element level which is crossed by geometry.
 * This function assumes that the octree is balanced 2:1
 * @param[in] level the level of the coarsest octree
 */
double LevelSetOctree::computeRSearchFromLevel( uint8_t level){
    return  (m_octree->getTree()).levelToSize(level) *sqrt(11.) /2. ;
};

/*!
 * Compute size of smallest octants greater than a given size (typically the size of the narrow band)
 * @param[in] r limit size
 */
double LevelSetOctree::computeSizeFromRSearch( double r){

    PabloUniform &tree = m_octree->getTree() ;

    uint8_t     level ( tree.getLocalMaxDepth() ) ;
    double      size ;

    size = tree.levelToSize(level) ;

    while( size <= r ) {
        level-- ;
        size = tree.levelToSize(level) ;
    };

    return size ;

};

}
