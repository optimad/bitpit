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
# include "communications.hpp"
# endif

# include <unordered_set>

# include "bitpit_SA.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
    @ingroup levelset
    @interface  LevelSetKernel
    @brief  Mesh specific implementation to calculate the levelset function

*/

/*!
 * Default constructor.
 */
LevelSetKernel::LevelSetKernel() {

    m_mesh = NULL ;

#if BITPIT_ENABLE_MPI
    m_commMPI = MPI_COMM_NULL;
# endif

};

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
LevelSetKernel::LevelSetKernel( VolumeKernel *patch): LevelSetKernel() {
    m_mesh = patch ;
};

/*!
 * Destructor of LevelSetKernel
*/
LevelSetKernel::~LevelSetKernel(){
    m_mesh = NULL ;

# if BITPIT_ENABLE_MPI
    freeCommunicator();
# endif

};

/*!
 * Returns pointer to underlying mesh.
 * @return pointer to mesh
*/
VolumeKernel* LevelSetKernel::getMesh() const{
    return m_mesh ;
} 

/*!
 * Compute the size of the narrow band using the levelset value
 * @param[in] visitor reference object
 * @param[in] signd indicates if signed or unsigned distances are calculated
 */
double LevelSetKernel::computeSizeNarrowBandFromLS( LevelSetObject *visitor, const bool &signd ){

    // We need to consider only the cells with a levelset value less than
    // local narrow band (ie. size of the narrowband evalauted using the
    // cell).
    double newRSearch = 0.;
    int factor  ;

    for (auto &cell : getMesh()->getCells() ) {
        // Discard cells outside the narrow band
        long id = cell.getId() ;
        if( !visitor->isInNarrowBand(id) ) {
            continue;
        }

        // Evaluate local search radius
        std::array<double,3> myCenter = computeCellCentroid(id) ;

        const long* neighbours = cell.getAdjacencies() ;
        int N = cell.getAdjacencyCount() ;

        double localRSearch = 0. ;
        for(int n=0; n<N; ++n){
            long neighId = neighbours[n] ;
            if (neighId < 0) {
                continue;
            }

            if( visitor->isInNarrowBand(neighId)){
                factor = (int) signd * visitor->getSign(id) + (int) (!signd) ;
                std::array<double,3> diff = computeCellCentroid(neighId) - myCenter ;
                if( factor *dotProduct(diff, visitor->getGradient(id)) < 0){
                    localRSearch = std::max( localRSearch, norm2(diff) );
                }
            }
        } 

        // Discard cells with a levelset greater than the local narrow band
        if ( std::abs( visitor->getLS(id) ) > (localRSearch + m_mesh->getTol()) ) {
            continue;
        }

        // Update the levelset
        newRSearch = std::max( localRSearch, newRSearch );
    }

# if BITPIT_ENABLE_MPI
    if( assureMPI() ) {
        double reducedRSearch ;
        MPI_Allreduce( &newRSearch, &reducedRSearch, 1, MPI_DOUBLE, MPI_MAX, m_commMPI );
        newRSearch = reducedRSearch ;
    }
#endif

    return newRSearch;

};

/*!
 * Clears the geometry cache.
 */
void LevelSetKernel::clearGeometryCache(  ) {

    std::unordered_map<long, std::array<double,3>>().swap( m_cellCentroids ) ;

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetKernel::updateGeometryCache( const std::vector<adaption::Info> &mapper ) {

    // If there are no cells in the mesh we can just delete all the cache
    if ( m_mesh->getCellCount() == 0) {
        clearGeometryCache();
        return;
    }

    // Remove the previous cells from the cache
    for ( auto & map : mapper ){
        if( map.entity != adaption::Entity::ENTITY_CELL ){
            continue;
        }

        for ( auto & previousId : map.previous){
            auto centroidItr = m_cellCentroids.find( previousId ) ;
            if ( centroidItr == m_cellCentroids.end() ) {
                continue ;
            }

            m_cellCentroids.erase( previousId ) ;
        }
    }

}

/*!
 * Computes the centroid of the specfified cell.
 *
 * If the centroid of the cell has been already evaluated, the cached value is
 * returned. Otherwise the cell centroid is evaluated and stored in the cache.
 *
 * @param[in] id is the index of cell
 * @return The centroid of the cell.
 */
const std::array<double,3> & LevelSetKernel::computeCellCentroid( long id ) {

    auto centroidItr = m_cellCentroids.find( id ) ;
    if ( centroidItr == m_cellCentroids.end() ) {
        centroidItr = m_cellCentroids.insert( { id, m_mesh->evalCellCentroid( id ) } ).first ;
    }

    return centroidItr->second;

}


/*!
 * Checks if the specified cell is inside the given bounding box
 * @param[in] id is the id of the cell
 * @param[in] minPoint is the lower left point of the boungind box
 * @param[in] maxPoint is the upper right point of the boungind box
 * @result True if the specified cell is inside the given bounding box, false
 * otherwise.
 */
double LevelSetKernel::isCellInsideBoundingBox( long id, std::array<double, 3> minPoint, std::array<double, 3> maxPoint ){

    const Cell &cell = m_mesh->getCell(id);
    double tolerance = m_mesh->getTol();

    std::array<double, 3> cellMinPoint;
    std::array<double, 3> cellMaxPoint;

    cellMinPoint.fill(   std::numeric_limits<double>::max() ) ;
    cellMaxPoint.fill( - std::numeric_limits<double>::max() ) ;

    int nVertices = cell.getVertexCount();
    for (int i = 0; i < nVertices; ++i) {
        long vertexId = cell.getVertex(i);
        std::array<double, 3> vertexCoords = m_mesh->getVertexCoords(vertexId);
        for (int d = 0; d < 3; ++d) {
            cellMinPoint[d] = std::min( vertexCoords[d] - tolerance, cellMinPoint[d]) ;
            cellMaxPoint[d] = std::max( vertexCoords[d] + tolerance, cellMaxPoint[d]) ;
        }
    }

    return CGElem::intersectBoxBox(minPoint, maxPoint, cellMinPoint, cellMaxPoint);
};

# if BITPIT_ENABLE_MPI

/*!
 * Returns the MPI communicator stored within LevelSetKernel
 * @return MPI communicator
 */
MPI_Comm LevelSetKernel::getCommunicator() const {
    return m_commMPI;
}

/*!
    Checks if the communicator to be used for parallel communications has
    already been set.

    \result Returns true if the communicator has been set, false otherwise.
*/
bool LevelSetKernel::isCommunicatorSet() const {

    return (getCommunicator() != MPI_COMM_NULL);
}

/*!
    Frees the MPI communicator associated to the patch
*/
void LevelSetKernel::freeCommunicator() {

    if (!isCommunicatorSet()) {
        return;
    }

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled) {
        return;
    }

    MPI_Comm_free(&m_commMPI);
}

/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available MPI communicator is retreived from mesh and duplicated if necessary and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSetKernel::assureMPI( ){

    if (isCommunicatorSet()) {
        return true ;
    }

    MPI_Comm meshComm = m_mesh->getCommunicator() ;
    if( meshComm == MPI_COMM_NULL){
        return false;
    } else {
        MPI_Comm_dup(m_mesh->getCommunicator(), &m_commMPI);
        return true;
    }
}

#endif

}
