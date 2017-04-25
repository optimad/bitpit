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

#include "element_reference.hpp"

namespace bitpit {

/*!
    \class ReferenceElementInfo
    \ingroup patchelements

    \brief The ReferenceElementInfo class allows to define information about
    reference elements.

    The local numbering scheme of element vertices is shown below.

    \image html common_elements.png
*/

const ReferenceElementInfo ReferenceElementInfo::undefinedInfo  = ReferenceElementInfo(ElementType::UNDEFINED);
const ReferenceElementInfo ReferenceElementInfo::vertexInfo     = ReferenceElementInfo(ElementType::VERTEX);
const ReferenceElementInfo ReferenceElementInfo::lineInfo       = ReferenceElementInfo(ElementType::LINE);
const ReferenceElementInfo ReferenceElementInfo::triangleInfo   = ReferenceElementInfo(ElementType::TRIANGLE);
const ReferenceElementInfo ReferenceElementInfo::pixelInfo      = ReferenceElementInfo(ElementType::PIXEL);
const ReferenceElementInfo ReferenceElementInfo::quadInfo       = ReferenceElementInfo(ElementType::QUAD);
const ReferenceElementInfo ReferenceElementInfo::tetraInfo      = ReferenceElementInfo(ElementType::TETRA);
const ReferenceElementInfo ReferenceElementInfo::voxelInfo      = ReferenceElementInfo(ElementType::VOXEL);
const ReferenceElementInfo ReferenceElementInfo::hexahedronInfo = ReferenceElementInfo(ElementType::HEXAHEDRON);
const ReferenceElementInfo ReferenceElementInfo::pyramidInfo    = ReferenceElementInfo(ElementType::PYRAMID);
const ReferenceElementInfo ReferenceElementInfo::wedgeInfo      = ReferenceElementInfo(ElementType::WEDGE);

/*!
    Default constructor
*/
ReferenceElementInfo::ReferenceElementInfo()
{
    initializeUndefinedInfo();
}

/*!
    Creates a new set of element information.

    \param type is the type of element
*/
ReferenceElementInfo::ReferenceElementInfo(ElementType type)
{
    switch (type) {

    case (ElementType::VERTEX):
        initializeVertexInfo();
        break;

    case (ElementType::LINE):
        initializeLineInfo();
        break;

    case (ElementType::TRIANGLE):
        initializeTriangleInfo();
        break;

    case (ElementType::PIXEL):
        initializePixelInfo();
        break;

    case (ElementType::QUAD):
        initializeQuadInfo();
        break;

    case (ElementType::TETRA):
        initializeTetraInfo();
        break;

    case (ElementType::VOXEL):
        initializeVoxelInfo();
        break;

    case (ElementType::HEXAHEDRON):
        initializeHexahedronInfo();
        break;

    case (ElementType::PYRAMID):
        initializePyramidInfo();
        break;

    case (ElementType::WEDGE):
        initializeWedgeInfo();
        break;

    default:
        initializeUndefinedInfo();
        break;

    }
}

/*!
    Gets the information for the specified element type.

    \param type is the type of element
    \result The information for the specified element type.
*/
const ReferenceElementInfo & ReferenceElementInfo::getInfo(ElementType type)
{
    switch (type) {

    case (ElementType::VERTEX):
        return vertexInfo;

    case (ElementType::LINE):
        return lineInfo;

    case (ElementType::TRIANGLE):
        return triangleInfo;

    case (ElementType::PIXEL):
        return pixelInfo;

    case (ElementType::QUAD):
        return quadInfo;

    case (ElementType::TETRA):
        return tetraInfo;

    case (ElementType::VOXEL):
        return voxelInfo;

    case (ElementType::HEXAHEDRON):
        return hexahedronInfo;

    case (ElementType::PYRAMID):
        return pyramidInfo;

    case (ElementType::WEDGE):
        return wedgeInfo;

    default:
        BITPIT_UNREACHABLE("Unsupported element");
        return undefinedInfo;

    }
}

/*!
    Initializes the information for the undefined element.
*/
void ReferenceElementInfo::initializeUndefinedInfo()
{
    type      = UNDEFINED;
    dimension = -1;

    nVertices = -1;
    nEdges    = -1;
    nFaces    = -1;
}

/*!
    Initializes the information for the point element.
*/
void ReferenceElementInfo::initializeVertexInfo()
{
    type      = ElementType::VERTEX;
    dimension = 0;

    // Vertices data
    nVertices = 1;

    // Edge data
    nEdges = 1;

    edge_type = std::vector<ElementType>(nEdges);
    edge_type[0] = ElementType::VERTEX;

    edgeConnect = std::vector<std::vector<int>>(nEdges);
    edgeConnect[0] = std::vector<int>(nVertices);
    edgeConnect[0][0] = 0;

    // Face data
    nFaces = 1;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    face_type[0] = ElementType::VERTEX;

    facesInfo[0] = this;

    faceConnect = std::vector<std::vector<int>>(nFaces);
    faceConnect[0] = std::vector<int>(nVertices);
    faceConnect[0][0] = 0;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the line element.
*/
void ReferenceElementInfo::initializeLineInfo()
{
    ReferenceElementInfo vertexInfo(VERTEX);

    type      = LINE;
    dimension = 1;

    // Vertices data
    nVertices = 2;

    // Edge data
    nEdges = 2;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    nFaces = 2;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]       = ElementType::VERTEX;
        facesInfo[k]      = &vertexInfo;
        faceConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        faceConnect[k][0] = k;
    }

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the triangle element.
*/
void ReferenceElementInfo::initializeTriangleInfo()
{
    ReferenceElementInfo vertexInfo(VERTEX);
    ReferenceElementInfo lineInfo(LINE);

    type      = TRIANGLE;
    dimension = 2;

    // Vertices data
    nVertices = 3;

    // Edge data
    nEdges = 3;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    nFaces = 3;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]       = LINE;
        facesInfo[k]      = &lineInfo;
        faceConnect[k]    = std::vector<int>(lineInfo.nVertices);
        faceConnect[k][0] = k;
        faceConnect[k][1] = (k + 1) % nVertices;
    }

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the quadrangle element.
*/
void ReferenceElementInfo::initializePixelInfo()
{
    ReferenceElementInfo vertexInfo(VERTEX);
    ReferenceElementInfo lineInfo(LINE);

    type      = PIXEL;
    dimension = 2;

    // Vertices data
    nVertices = 4;

    // Edge data
    nEdges = 4;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    nFaces = 4;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]    = LINE;
        facesInfo[k]   = &lineInfo;
        faceConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    faceConnect[0][0] = 2;
    faceConnect[0][1] = 0;

    faceConnect[1][0] = 1;
    faceConnect[1][1] = 3;

    faceConnect[2][0] = 0;
    faceConnect[2][1] = 1;

    faceConnect[3][0] = 3;
    faceConnect[3][1] = 2;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the quadrangle element.
*/
void ReferenceElementInfo::initializeQuadInfo()
{
    ReferenceElementInfo vertexInfo(VERTEX);
    ReferenceElementInfo lineInfo(LINE);

    type      = QUAD;
    dimension = 2;

    // Vertices data
    nVertices = 4;

    // Edge data
    nEdges = 4;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    nFaces = 4;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]       = LINE;
        facesInfo[k]      = &lineInfo;
        faceConnect[k]    = std::vector<int>(lineInfo.nVertices);
        faceConnect[k][0] = k;
        faceConnect[k][1] = (k + 1) % nVertices;
    }

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the tetrahedron element.
*/
void ReferenceElementInfo::initializeTetraInfo()
{
    ReferenceElementInfo lineInfo(LINE);
    ReferenceElementInfo triangleInfo(TRIANGLE);

    type      = TETRA;
    dimension = 3;

    // Vertices data
    nVertices = 4;

    // Edge data
    nEdges = 6;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]   = LINE;
        edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    edgeConnect[0][0] = 0;
    edgeConnect[0][1] = 1;

    edgeConnect[1][0] = 1;
    edgeConnect[1][1] = 2;

    edgeConnect[2][0] = 2;
    edgeConnect[2][1] = 0;

    edgeConnect[3][0] = 3;
    edgeConnect[3][1] = 0;

    edgeConnect[4][0] = 3;
    edgeConnect[4][1] = 1;

    edgeConnect[5][0] = 3;
    edgeConnect[5][1] = 2;

    // Face data
    nFaces = 4;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]    = TRIANGLE;
        facesInfo[k]   = &triangleInfo;
        faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
    }

    faceConnect[0][0] = 1;
    faceConnect[0][1] = 0;
    faceConnect[0][2] = 2;

    faceConnect[1][0] = 0;
    faceConnect[1][1] = 3;
    faceConnect[1][2] = 2;

    faceConnect[2][0] = 3;
    faceConnect[2][1] = 1;
    faceConnect[2][2] = 2;

    faceConnect[3][0] = 0;
    faceConnect[3][1] = 1;
    faceConnect[3][2] = 3;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the brick element.
*/
void ReferenceElementInfo::initializeVoxelInfo()
{
    ReferenceElementInfo lineInfo(LINE);
    ReferenceElementInfo pixelInfo(PIXEL);

    type      = VOXEL;
    dimension = 3;

    // Vertices data
    nVertices = 8;

    // Edge data
    nEdges = 12;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]   = LINE;
        edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    edgeConnect[0][0] = 0;
    edgeConnect[0][1] = 2;

    edgeConnect[1][0] = 1;
    edgeConnect[1][1] = 3;

    edgeConnect[2][0] = 0;
    edgeConnect[2][1] = 1;

    edgeConnect[3][0] = 2;
    edgeConnect[3][1] = 3;

    edgeConnect[4][0] = 0;
    edgeConnect[4][1] = 4;

    edgeConnect[5][0] = 1;
    edgeConnect[5][1] = 5;

    edgeConnect[6][0] = 2;
    edgeConnect[6][1] = 6;

    edgeConnect[7][0] = 3;
    edgeConnect[7][1] = 7;

    edgeConnect[8][0] = 4;
    edgeConnect[8][1] = 6;

    edgeConnect[9][0] = 5;
    edgeConnect[9][1] = 7;

    edgeConnect[10][0] = 4;
    edgeConnect[10][1] = 5;

    edgeConnect[11][0] = 6;
    edgeConnect[11][1] = 7;

    // Face data
    nFaces = 6;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]    = PIXEL;
        facesInfo[k]   = &pixelInfo;
        faceConnect[k] = std::vector<int>(pixelInfo.nVertices);
    }

    faceConnect[0][0] = 2;
    faceConnect[0][1] = 0;
    faceConnect[0][2] = 6;
    faceConnect[0][3] = 4;

    faceConnect[1][0] = 1;
    faceConnect[1][1] = 3;
    faceConnect[1][2] = 5;
    faceConnect[1][3] = 7;

    faceConnect[2][0] = 0;
    faceConnect[2][1] = 1;
    faceConnect[2][2] = 4;
    faceConnect[2][3] = 5;

    faceConnect[3][0] = 3;
    faceConnect[3][1] = 2;
    faceConnect[3][2] = 7;
    faceConnect[3][3] = 6;

    faceConnect[4][0] = 2;
    faceConnect[4][1] = 3;
    faceConnect[4][2] = 0;
    faceConnect[4][3] = 1;

    faceConnect[5][0] = 4;
    faceConnect[5][1] = 5;
    faceConnect[5][2] = 6;
    faceConnect[5][3] = 7;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the hexahedron element.
*/
void ReferenceElementInfo::initializeHexahedronInfo()
{
    ReferenceElementInfo lineInfo(LINE);
    ReferenceElementInfo quadInfo(QUAD);

    type      = HEXAHEDRON;
    dimension = 3;

    // Vertices data
    nVertices = 8;

    // Edge data
    nEdges = 12;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]   = LINE;
        edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    edgeConnect[0][0] = 1;
    edgeConnect[0][1] = 0;

    edgeConnect[1][0] = 1;
    edgeConnect[1][1] = 2;

    edgeConnect[2][0] = 2;
    edgeConnect[2][1] = 3;

    edgeConnect[3][0] = 3;
    edgeConnect[3][1] = 0;

    edgeConnect[4][0] = 4;
    edgeConnect[4][1] = 5;

    edgeConnect[5][0] = 5;
    edgeConnect[5][1] = 6;

    edgeConnect[6][0] = 6;
    edgeConnect[6][1] = 7;

    edgeConnect[7][0] = 7;
    edgeConnect[7][1] = 4;

    edgeConnect[8][0] = 0;
    edgeConnect[8][1] = 4;

    edgeConnect[9][0] = 1;
    edgeConnect[9][1] = 5;

    edgeConnect[10][0] = 2;
    edgeConnect[10][1] = 6;

    edgeConnect[11][0] = 3;
    edgeConnect[11][1] = 7;

    // Face data
    nFaces = 6;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        face_type[k]    = QUAD;
        facesInfo[k]   = &quadInfo;
        faceConnect[k] = std::vector<int>(quadInfo.nVertices);
    }

    faceConnect[0][0] = 1;
    faceConnect[0][1] = 0;
    faceConnect[0][2] = 3;
    faceConnect[0][3] = 2;

    faceConnect[1][0] = 4;
    faceConnect[1][1] = 5;
    faceConnect[1][2] = 6;
    faceConnect[1][3] = 7;

    faceConnect[2][0] = 7;
    faceConnect[2][1] = 3;
    faceConnect[2][2] = 0;
    faceConnect[2][3] = 4;

    faceConnect[3][0] = 5;
    faceConnect[3][1] = 1;
    faceConnect[3][2] = 2;
    faceConnect[3][3] = 6;

    faceConnect[4][0] = 4;
    faceConnect[4][1] = 0;
    faceConnect[4][2] = 1;
    faceConnect[4][3] = 5;

    faceConnect[5][0] = 6;
    faceConnect[5][1] = 2;
    faceConnect[5][2] = 3;
    faceConnect[5][3] = 7;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the pyramid element.
*/
void ReferenceElementInfo::initializePyramidInfo()
{
    ReferenceElementInfo lineInfo(LINE);
    ReferenceElementInfo triangleInfo(TRIANGLE);
    ReferenceElementInfo quadInfo(QUAD);

    type      = PYRAMID;
    dimension = 3;

    // Vertices data
    nVertices = 5;

    // Edge data
    nEdges = 8;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]   = LINE;
        edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    edgeConnect[0][0] = 0;
    edgeConnect[0][1] = 1;

    edgeConnect[1][0] = 1;
    edgeConnect[1][1] = 2;

    edgeConnect[2][0] = 2;
    edgeConnect[2][1] = 3;

    edgeConnect[3][0] = 3;
    edgeConnect[3][1] = 0;

    edgeConnect[4][0] = 4;
    edgeConnect[4][1] = 0;

    edgeConnect[5][0] = 4;
    edgeConnect[5][1] = 1;

    edgeConnect[6][0] = 4;
    edgeConnect[6][1] = 2;

    edgeConnect[7][0] = 4;
    edgeConnect[7][1] = 3;

    // Face data
    nFaces = 5;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        if (k == 0) {
            face_type[k]   = ElementType::QUAD;
            facesInfo[k]   = &quadInfo;
            faceConnect[k] = std::vector<int>(quadInfo.nVertices);
        } else {
            face_type[k]   = ElementType::TRIANGLE;
            facesInfo[k]   = &triangleInfo;
            faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
        }
    }

    faceConnect[0][0] = 0;
    faceConnect[0][1] = 3;
    faceConnect[0][2] = 2;
    faceConnect[0][3] = 1;

    faceConnect[1][0] = 3;
    faceConnect[1][1] = 0;
    faceConnect[1][2] = 4;

    faceConnect[2][0] = 0;
    faceConnect[2][1] = 1;
    faceConnect[2][2] = 4;

    faceConnect[3][0] = 1;
    faceConnect[3][1] = 2;
    faceConnect[3][2] = 4;

    faceConnect[4][0] = 2;
    faceConnect[4][1] = 3;
    faceConnect[4][2] = 4;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the information for the wedge (triangular prism) element.
*/
void ReferenceElementInfo::initializeWedgeInfo()
{
    ReferenceElementInfo lineInfo(LINE);
    ReferenceElementInfo triangleInfo(TRIANGLE);
    ReferenceElementInfo quadInfo(QUAD);

    type      = WEDGE;
    dimension = 3;

    // Vertices data
    nVertices = 6;

    // Edge data
    nEdges = 9;

    edge_type = std::vector<ElementType>(nEdges);
    edgeConnect = std::vector<std::vector<int>>(nEdges);
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]   = LINE;
        edgeConnect[k] = std::vector<int>(lineInfo.nVertices);
    }

    edgeConnect[0][0] = 1;
    edgeConnect[0][1] = 0;

    edgeConnect[1][0] = 1;
    edgeConnect[1][1] = 2;

    edgeConnect[2][0] = 2;
    edgeConnect[2][1] = 0;

    edgeConnect[3][0] = 3;
    edgeConnect[3][1] = 4;

    edgeConnect[4][0] = 4;
    edgeConnect[4][1] = 5;

    edgeConnect[5][0] = 5;
    edgeConnect[5][1] = 3;

    edgeConnect[6][0] = 3;
    edgeConnect[6][1] = 0;

    edgeConnect[7][0] = 4;
    edgeConnect[7][1] = 1;

    edgeConnect[8][0] = 5;
    edgeConnect[8][1] = 2;

    // Face data
    nFaces = 5;

    std::vector<ReferenceElementInfo *> facesInfo(nFaces);

    face_type = std::vector<ElementType>(nFaces);
    faceConnect = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        if (k == 0 || k == 1) {
            face_type[k]    = TRIANGLE;
            facesInfo[k]   = &triangleInfo;
            faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
        } else {
            face_type[k]    = QUAD;
            facesInfo[k]   = &quadInfo;
            faceConnect[k] = std::vector<int>(quadInfo.nVertices);
        }
    }

    faceConnect[0][0] = 1;
    faceConnect[0][1] = 0;
    faceConnect[0][2] = 2;

    faceConnect[1][0] = 3;
    faceConnect[1][1] = 4;
    faceConnect[1][2] = 5;

    faceConnect[2][0] = 3;
    faceConnect[2][1] = 0;
    faceConnect[2][2] = 1;
    faceConnect[2][3] = 4;

    faceConnect[3][0] = 4;
    faceConnect[3][1] = 1;
    faceConnect[3][2] = 2;
    faceConnect[3][3] = 5;

    faceConnect[4][0] = 5;
    faceConnect[4][1] = 2;
    faceConnect[4][2] = 0;
    faceConnect[4][3] = 3;

    initializeFaceEdges(facesInfo);
}

/*!
    Initializes the list of edges associated to the faces
*/
void ReferenceElementInfo::initializeFaceEdges(const std::vector<ReferenceElementInfo *> &facesInfo)
{
    faceEdges = std::vector<std::vector<int>>(nFaces);
    for (int k = 0; k < nFaces; ++k) {
        const ReferenceElementInfo &faceInfo = *(facesInfo[k]);

        int nFaceEdges = faceInfo.nFaces;
        for (int i = 0; i < nFaceEdges; ++i) {
            // Number of vertices of the edge associated to the face
            std::size_t nFaceEdgeVertices = faceInfo.faceConnect[i].size();

            // Connectivity of the edge associated to the face
            const std::vector<int> &localFaceEdgeConnect = faceInfo.faceConnect[i];

            std::vector<int> faceEdgeConnect(nFaceEdgeVertices);
            for (std::size_t n = 0; n < nFaceEdgeVertices; ++n) {
                int localVertexId = localFaceEdgeConnect[n];
                int vertexId      = faceConnect[k][localVertexId];

                faceEdgeConnect[n] = vertexId;
            }

            // Search the edge that has the same connectivity of the face edge
            for (int j = 0; j < nEdges; ++j) {
                // If face edge and the guess edge have a different number of
                // vertices, the two edge cannot be the same.
                std::size_t nGuessEdgeVertices = edgeConnect[j].size();
                if (nGuessEdgeVertices != nFaceEdgeVertices) {
                    continue;
                }

                // If the connecitivity of the face edge and the one of the
                // guess edge are the same, the two edges coincides.
                const std::vector<int> commonVertices = utils::intersectionVector(faceEdgeConnect, edgeConnect[j]);
                if (commonVertices.size() == nFaceEdgeVertices) {
                    faceEdges[k].push_back(j);
                }
            }
        }

        assert(faceEdges[k].size() == nFaceEdges);
    }
}

}
