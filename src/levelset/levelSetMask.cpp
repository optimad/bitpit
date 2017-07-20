/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# include <cassert>
# include <memory>

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_surfunstructured.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"

# include "levelSetObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetSegmentation.hpp"
# include "levelSetMask.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@class  LevelSetMask
	@brief Implements the levelset around a set of cells or interfaces of the kernel     
*/

/*!
 * Destructor
 */
LevelSetMask::~LevelSetMask() {
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] mask the list of inner cells
 * @param[in] mesh the mesh hosting the cells
 */
LevelSetMask::LevelSetMask(int id, const std::unordered_set<long> &mask, const VolumeKernel &mesh ) :LevelSetSegmentation(id) {

    std::unordered_map<long,long> meshToEnvelope ;

    std::unique_ptr<SurfUnstructured> segmentation = extractCellEnvelope(mask,mesh,meshToEnvelope) ;

    long intrIndex = meshToEnvelope.begin()->first;
    long enveIndex = meshToEnvelope.begin()->second;

    bool sameOrientation = sameInterfaceEnvelopeOrientation(mesh, intrIndex, *segmentation, enveIndex);

    auto const &interface = mesh.getInterface(intrIndex);
    long ownerId = interface.getOwner();
    bool invert = (mask.count(ownerId)==0) ;

    bool flip = (sameOrientation == invert);

    bool orientable = segmentation->adjustCellOrientation( enveIndex, flip);
    if( !orientable){
        throw std::runtime_error ("Error in LevelSetMask");
    }

    setSegmentation(std::move(segmentation));
}

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] list the list of interfaces
 * @param[in] refId id of reference interface
 * @param[in] invert if orientation should be inverted with respect to the reference interface
 * @param[in] mesh the mesh hosting the cells
 */
LevelSetMask::LevelSetMask(int id, const std::vector<long> &list, const long &intrIndex, const bool &invert, const VolumeKernel &mesh ) :LevelSetSegmentation(id) {

    std::unordered_map<long,long> meshToEnvelope ;

    std::unique_ptr<SurfUnstructured> segmentation = extractFaceEnvelope(list,mesh,meshToEnvelope) ;

    long enveIndex = meshToEnvelope.at(intrIndex);
    bool sameOrientation = sameInterfaceEnvelopeOrientation(mesh, intrIndex, *segmentation, enveIndex);

    bool flip = (sameOrientation == invert);

    bool orientable = segmentation->adjustCellOrientation( enveIndex, flip);
    if( !orientable){
        throw std::runtime_error ("Error in LevelSetMask");
    }

    setSegmentation(std::move(segmentation));
}

/*!
 * Extracts the external envelope and create a new patch from it
 * The external envelope is composed by all the outer faces of the masked cells
 * @param[in] mask the list of inner cells
 * @param[in] mesh the mesh hosting the cells
 * @param[in] meshToEnvelope map which hosts the ndex association between the cells of the envelope and the faces of the mesh.
 * If the nullptr is passed, a local map will be used.
 * @return surface mesh
*/
std::unique_ptr<SurfUnstructured> LevelSetMask::extractCellEnvelope(const std::unordered_set<long> &mask, const VolumeKernel &mesh, std::unordered_map<long,long> &meshToEnvelope){

    std::vector<long> list;

    for( const long &cellIndex : mask){
        const auto &cell = mesh.getCell(cellIndex);
        const long *adjacencies = cell.getAdjacencies();
        const long *interfaceIndex = cell.getInterfaces();

        int interfaceCount= cell.getInterfaceCount() ;

        for(int i=0; i<interfaceCount; i++){
            const long &neigh = adjacencies[i];

            if(mask.count(neigh)==0){
                list.push_back(interfaceIndex[i]);
            }
        }
    }

    return extractFaceEnvelope(list,mesh,meshToEnvelope);
}

/*!
 * Extracts the external envelope and create a new patch from it.
 * The external envelope is composed by all interfaces in the list.
 * @param[in] list the list of interfaces
 * @param[in] mesh the mesh hosting the cells
 * @param[in] meshToEnvelope map which hosts the index association between the cells of the envelope and the faces of the mesh.
 * @return surface mesh
*/
std::unique_ptr<SurfUnstructured> LevelSetMask::extractFaceEnvelope(const std::vector<long> &list, const VolumeKernel &mesh, std::unordered_map<long,long> &meshToEnvelope){

    std::unique_ptr<SurfUnstructured> envelope = std::unique_ptr<SurfUnstructured>(new SurfUnstructured(0,mesh.getDimension()-1,mesh.getDimension()));

	// ====================================================================== //
	// RESIZE DATA STRUCTURES                                                 //
	// ====================================================================== //
    long nVertices(0);
    long nCells(list.size());

    for( const long &faceIndex : list){
        auto const &interface = mesh.getInterface(faceIndex);
        nVertices += interface.getVertexCount() ;
    }

	envelope->reserveVertices(nVertices);
	envelope->reserveCells(nCells);

	// ====================================================================== //
	// LOOP OVER CELLS                                                        //
	// ====================================================================== //
	std::unordered_map<long,long> vertexMap;

    for( const long &faceIndex : list){
        auto const &interface = mesh.getInterface(faceIndex);
        const long *faceConnect = interface.getConnect();
        int nFaceVertices = interface.getVertexCount();

        // Add face vertices to the envelope and get face
        // connectivity in the envelope
        std::unique_ptr<long[]> faceEnvelopeConnect = std::unique_ptr<long[]>(new long[nFaceVertices]);
        for (int j = 0; j < nFaceVertices; ++j) {
        	long vertexId = faceConnect[j];
        
        	// If the vertex is not yet in the envelope
        	// add it.
        	if (vertexMap.count(vertexId) == 0) {
        		const Vertex &vertex = mesh.getVertex(vertexId);
        		auto envelopeVertex = envelope->addVertex(vertex);
        		vertexMap[vertexId] = envelopeVertex->getId();
        	}
        
        	// Update face connectivity in the envelope
        	faceEnvelopeConnect[j] = vertexMap.at(vertexId);
        }

        // Add face to envelope
        ElementInfo::Type faceType = interface.getType();
        PatchKernel::CellIterator cellItr = envelope->addCell(faceType, true, std::move(faceEnvelopeConnect));
        meshToEnvelope.insert({{faceIndex,cellItr->getId()}});
	}

    envelope->squeeze();
    envelope->buildAdjacencies();
    envelope->buildInterfaces();

    envelope->getVTK().setName("geometry_002") ;
    envelope->write() ;

    return envelope;

}

/*!
 * Checks if the the corresponding mesh interface and envelope cell have the same orientation
 * @param[in] mesh the mesh hosting the cells
 * @param[in] faceIndex index of the mesh interface
 * @param[in] envelope surface mesh
 * @param[in] enveIndex index of the envelope cell
 * @return true if interface and cell have the same orientation
*/
bool LevelSetMask::sameInterfaceEnvelopeOrientation(const VolumeKernel &mesh, const long &faceIndex, SurfUnstructured &envelope, const long &enveIndex){

    std::array<double,3> facetNormal = envelope.evalFacetNormal(enveIndex);
    std::array<double,3> interfaceNormal = mesh.evalInterfaceNormal(faceIndex);
    
    return (dotProduct(facetNormal,interfaceNormal)>0) ;

}

}

