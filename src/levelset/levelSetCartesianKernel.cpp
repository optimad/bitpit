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

# include "bitpit_volcartesian.hpp"
# include "bitpit_CG.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesianKernel.hpp"

namespace bitpit {

/*!
	@ingroup    levelset
	@class      LevelSetCartesianKernel
	@brief      Implements LevelSetKernel for cartesian meshes
*/

/*!
 * Constructor
 */
LevelSetCartesianKernel::LevelSetCartesianKernel(VolCartesian &patch ): LevelSetKernel( &patch ){

    // Get mesh information
    const VolCartesian *mesh = getMesh();
    int dimension = mesh->getDimension();
    std::array<double, 3> spacing = mesh->getSpacing();

    // Initialize bounding and tangent radii
    m_cellTangentRadius  = std::numeric_limits<double>::max();
    m_cellBoundingRadius = 0.;
    for(int i=0; i<dimension; ++i){
        m_cellTangentRadius   = std::min(m_cellTangentRadius, spacing[i]);
        m_cellBoundingRadius += spacing[i] * spacing[i];
    }
    m_cellTangentRadius  = 0.5 * m_cellTangentRadius;
    m_cellBoundingRadius = 0.5 * std::sqrt(m_cellBoundingRadius);

}

/*!
 * Returns a pointer to VolCartesian
 * @return pointer to VolCartesian
 */
VolCartesian * LevelSetCartesianKernel::getMesh() const{

    return static_cast<VolCartesian *>(LevelSetKernel::getMesh()) ;

}

/*!
 * Computes the radius of the tangent sphere associated with the cells of the mesh.
 *
 * The tangent sphere is a sphere having the center in the cell centroid and tangent
 * to the cell.
 *
 * @param[in] id is the id of cell
 * @return The radius of the tangent sphere.
 */
double LevelSetCartesianKernel::getCellTangentRadius( ) const {

    return m_cellTangentRadius;

}

/*!
 * Get the radius of the bounding sphere associated with the cells of the mesh.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * @param[in] id is the id of cell
 * @return The radius of the bounding sphere.
 */
double LevelSetCartesianKernel::getCellBoundingRadius( ) const {

    return m_cellBoundingRadius;

}

/*!
 * Computes the centroid of the specfified cell.
 *
 * @param[in] id is the id of cell
 * @return The centroid of the cell.
 */
std::array<double, 3> LevelSetCartesianKernel::computeCellCentroid( long id ) const {

    const VolCartesian *mesh = getMesh();

    return mesh->evalCellCentroid( id );
}

/*!
 * Computes the radius of the tangent sphere associated with the specified cell.
 *
 * The tangent sphere is a sphere having the center in the cell centroid and tangent
 * to the cell.
 *
 * @param[in] id is the id of cell
 * @return The radius of the tangent sphere.
 */
double LevelSetCartesianKernel::computeCellTangentRadius( long id ) const {

    BITPIT_UNUSED( id );

    return getCellTangentRadius();

}

/*!
 * Computes the radius of the bounding sphere associated with the specified cell.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * @param[in] id is the id of cell
 * @return The radius of the bounding sphere.
 */
double LevelSetCartesianKernel::computeCellBoundingRadius( long id ) const {

    BITPIT_UNUSED( id );

    return getCellBoundingRadius();

}

}
