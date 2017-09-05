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

/*!
    Constructor

    \param _dimension is the space dimension of the element
    \param _type is the type of element
    \param _nVertices is the number of vertices
    \param _nFaces is the number of faces
    \param _nEdges is the number of edges
*/
ReferenceElementInfo::ReferenceElementInfo(int _dimension, ElementType _type, int _nVertices, int _nFaces, int _nEdges)
    : dimension(_dimension), type(_type),
      nVertices(_nVertices),
      nFaces(_nFaces), face_type(nFaces), faceConnect(nFaces), faceEdges(nFaces),
      nEdges(_nEdges), edge_type(nEdges), edgeConnect(nEdges)
{
}

/*!
    Destructor
*/
ReferenceElementInfo::~ReferenceElementInfo()
{
}

/*!
    Check if the sepcified element type is associated to a reference element.

    \param type is the type of element
    \result Return tru if the sepcified element type is associated to a
    reference element, false otherwis.
*/
bool ReferenceElementInfo::hasInfo(ElementType type)
{

    switch (type) {

    case (ElementType::VERTEX):
    case (ElementType::LINE):
    case (ElementType::TRIANGLE):
    case (ElementType::PIXEL):
    case (ElementType::QUAD):
    case (ElementType::TETRA):
    case (ElementType::VOXEL):
    case (ElementType::HEXAHEDRON):
    case (ElementType::PYRAMID):
    case (ElementType::WEDGE):
        return true;

    default:
        return false;

    }
}

/*!
    Gets the information for the specified element type.

    \param type is the type of element
    \result The information for the specified element type.
*/
const ReferenceElementInfo & ReferenceElementInfo::getInfo(ElementType type)
{
    const static ReferenceVertexInfo vertexInfo;
    const static ReferenceLineInfo lineInfo;
    const static ReferenceTriangleInfo triangleInfo;
    const static ReferencePixelInfo pixelInfo;
    const static ReferenceQuadInfo quadInfo;
    const static ReferenceTetraInfo tetraInfo;
    const static ReferenceVoxelInfo voxelInfo;
    const static ReferenceHexahedronInfo hexahedronInfo;
    const static ReferencePyramidInfo pyramidInfo;
    const static ReferenceWedgeInfo wedgeInfo;

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

    }
}

/*!
    Initializes the list of edges associated to the faces
*/
void ReferenceElementInfo::initializeFaceEdges(const std::vector<const ReferenceElementInfo *> &facesInfo)
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

/*!
    \class Reference3DElementInfo
    \ingroup patchelements

    \brief The Reference3DElementInfo class allows to define information about
    reference three-dimensional elements.

    The local numbering scheme of element vertices is shown below.

    \image html common_elements.png
*/

/*!
    Constructor

    \param type is the type of element
    \param nVertices is the number of vertices
    \param nFaces is the number of faces
*/
Reference3DElementInfo::Reference3DElementInfo(ElementType type, int nVertices, int nFaces)
    : ReferenceElementInfo(3, type, nVertices, nFaces, nVertices + nFaces - 2)
{
}

/*!
    \ingroup patchelements

    \brief The ReferenceTetraInfo class defines the information about the
    reference tetrahedron.
*/

/*!
    Default constructor
*/
ReferenceTetraInfo::ReferenceTetraInfo()
    : Reference3DElementInfo(ElementType::TETRA, 4, 4)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceTriangleInfo triangleInfo;

    // Edge data
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
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        face_type[k]   = TRIANGLE;
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
    \class ReferenceVoxelInfo
    \ingroup patchelements

    \brief The ReferenceVoxelInfo class defines the information about the
    reference voxel.
*/

/*!
    Default constructor
*/
ReferenceVoxelInfo::ReferenceVoxelInfo()
    : Reference3DElementInfo(ElementType::VOXEL, 8, 6)
{
    const ReferenceLineInfo lineInfo;
    const ReferencePixelInfo pixelInfo;

    // Edge data
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
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        face_type[k]   = PIXEL;
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
    \class ReferenceHexahedronInfo
    \ingroup patchelements

    \brief The ReferenceHexahedronInfo class defines the information about the
    reference hexahedron.
*/

/*!
    Default constructor
*/
ReferenceHexahedronInfo::ReferenceHexahedronInfo()
    : Reference3DElementInfo(ElementType::HEXAHEDRON, 8, 6)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceQuadInfo quadInfo;

    // Edge data
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
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

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
    \class ReferencePyramidInfo
    \ingroup patchelements

    \brief The ReferencePyramidInfo class defines the information about the
    reference pyramid.
*/

/*!
    Default constructor
*/
ReferencePyramidInfo::ReferencePyramidInfo()
    : Reference3DElementInfo(ElementType::PYRAMID, 5, 5)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceTriangleInfo triangleInfo;
    const ReferenceQuadInfo quadInfo;

    // Edge data
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
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

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
    \class ReferenceWedgeInfo
    \ingroup patchelements

    \brief The ReferenceWedgeInfo class defines the information about the
    reference wedge.
*/

/*!
    Default constructor
*/
ReferenceWedgeInfo::ReferenceWedgeInfo()
    : Reference3DElementInfo(ElementType::WEDGE, 6, 5)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceTriangleInfo triangleInfo;
    const ReferenceQuadInfo quadInfo;

    // Edge data
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
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        if (k == 0 || k == 1) {
            face_type[k]   = TRIANGLE;
            facesInfo[k]   = &triangleInfo;
            faceConnect[k] = std::vector<int>(triangleInfo.nVertices);
        } else {
            face_type[k]   = QUAD;
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
    \class Reference2DElementInfo
    \ingroup patchelements

    \brief The Reference2DElementInfo class allows to define information about
    reference two-dimensional elements.

    The local numbering scheme of element vertices is shown below.

    \image html common_elements.png
*/

/*!
    Constructor

    \param type is the type of element
    \param nVertices is the number of vertices
    \param nFaces is the number of faces
*/
Reference2DElementInfo::Reference2DElementInfo(ElementType type, int nVertices)
    : ReferenceElementInfo(2, type, nVertices, nVertices, nVertices)
{
}

/*!
    \class ReferenceTriangleInfo
    \ingroup patchelements

    \brief The ReferenceTriangleInfo class defines the information about the
    reference triangle.
*/

/*!
    Default constructor
*/
ReferenceTriangleInfo::ReferenceTriangleInfo()
    : Reference2DElementInfo(ElementType::TRIANGLE, 3)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

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
    \class ReferencePixelInfo
    \ingroup patchelements

    \brief The ReferencePixelInfo class defines the information about the
    reference pixel.
*/

/*!
    Default constructor
*/
ReferencePixelInfo::ReferencePixelInfo()
    : Reference2DElementInfo(ElementType::PIXEL, 4)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

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
    \class ReferenceQuadInfo
    \ingroup patchelements

    \brief The ReferenceQuadInfo class defines the information about the
    reference quadrangle.
*/

/*!
    Default constructor
*/
ReferenceQuadInfo::ReferenceQuadInfo()
    : Reference2DElementInfo(ElementType::QUAD, 4)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

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
    \class Reference1DElementInfo
    \ingroup patchelements

    \brief The Reference1DElementInfo class allows to define information about
    reference one-dimensional elements.

    The local numbering scheme of element vertices is shown below.

    \image html common_elements.png
*/

/*!
    Constructor

    \param type is the type of element
*/
Reference1DElementInfo::Reference1DElementInfo(ElementType type)
    : ReferenceElementInfo(1, type, 2, 2, 2)
{
}

/*!
    \class ReferenceLineInfo
    \ingroup patchelements

    \brief The ReferenceLineInfo class defines the information about the
    reference line.
*/

/*!
    Default constructor
*/
ReferenceLineInfo::ReferenceLineInfo()
    : Reference1DElementInfo(ElementType::LINE)
{
    const ReferenceVertexInfo vertexInfo;

    // Edge data
    for (int k = 0; k < nEdges; ++k) {
        edge_type[k]      = ElementType::VERTEX;
        edgeConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        edgeConnect[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        face_type[k]       = ElementType::VERTEX;
        facesInfo[k]      = &vertexInfo;
        faceConnect[k]    = std::vector<int>(vertexInfo.nVertices);
        faceConnect[k][0] = k;
    }

    initializeFaceEdges(facesInfo);
}

/*!
    \class Reference0DElementInfo
    \ingroup patchelements

    \brief The Reference0DElementInfo class allows to define information about
    reference zero-dimensional elements.

    The local numbering scheme of element vertices is shown below.

    \image html common_elements.png
*/

/*!
    Constructor

    \param type is the type of element
*/
Reference0DElementInfo::Reference0DElementInfo(ElementType type)
    : ReferenceElementInfo(0, type, 1, 1, 1)
{
}

/*!
    \class ReferenceVertexInfo
    \ingroup patchelements

    \brief The ReferenceVertexInfo class defines the information about the
    reference vertex.
*/

/*!
    Default constructor
*/
ReferenceVertexInfo::ReferenceVertexInfo()
    : Reference0DElementInfo(ElementType::VERTEX)
{
    // Edge data
    edge_type[0] = ElementType::VERTEX;

    edgeConnect[0] = std::vector<int>(nVertices);
    edgeConnect[0][0] = 0;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    facesInfo[0] = this;

    face_type[0] = ElementType::VERTEX;

    faceConnect[0] = std::vector<int>(nVertices);
    faceConnect[0][0] = 0;

    initializeFaceEdges(facesInfo);
}

}
