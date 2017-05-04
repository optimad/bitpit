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

#include "bitpit_operators.hpp"

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
        throw std::runtime_error("Unsupported element");

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
    Evaluates the surface area of a three-dimensional element with the
    specified vertex coordinates. The surface area is a measure of the
    total area that the surface of the object occupies.

    \param vertexCoords are the coordinate of the vertices
    \result The area of the element.
*/
double Reference3DElementInfo::evalSurfaceArea(const std::array<double, 3> *vertexCoords) const
{
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> faceVertexCoords;

    double area = 0;
    for (int i = 0; i < nFaces; ++i) {
        for (int n = 0; n < nFaces; ++n) {
            faceVertexCoords[i] = vertexCoords[faceConnect[i][n]];
        }

        const Reference2DElementInfo &faceInfo = static_cast<const Reference2DElementInfo &>(getInfo(face_type[i]));
        area += faceInfo.evalArea(faceVertexCoords.data());
    }

    return area;
}

/*!
    Evaluates the perimeter of a three-dimensional element with the specified
    vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The perimeter of the element.
*/
double Reference3DElementInfo::evalEdgePerimeter(const std::array<double, 3> *vertexCoords) const
{
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> edgeVertexCoords;
    const ReferenceLineInfo &edgeInfo = static_cast<const ReferenceLineInfo &>(getInfo(ElementType::LINE));

    double perimeter = 0;
    for (int i = 0; i < nEdges; ++i) {
        for (int n = 0; n < edgeInfo.nVertices; ++n) {
            edgeVertexCoords[n] = vertexCoords[edgeConnect[i][n]];
        }

        perimeter += edgeInfo.evalLength(edgeVertexCoords.data());
    }

    return perimeter;
}

/*!
    \class ReferenceTetraInfo
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
    Evaluates the volume of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The volume of the element.
*/
double ReferenceTetraInfo::evalVolume(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];
    const std::array<double, 3> &V_D = vertexCoords[3];

    double volume = std::abs(dotProduct(V_A - V_D, crossProduct(V_B - V_D, V_C - V_D))) / 6.;

    return volume;
}

/*!
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the tetrahedron is the height of the regular
    tetrahedron having the same inscribed sphere.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceTetraInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double volume = evalVolume(vertexCoords);
    double area   = evalSurfaceArea(vertexCoords);

    double inscribedRadius = 3 * volume / area;

    double length = 4 * inscribedRadius * sqrt(3. / 2.);

    return length;
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
    Evaluates the volume of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The volume of the element.
*/
double ReferenceVoxelInfo::evalVolume(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];
    const std::array<double, 3> &V_D = vertexCoords[3];

    double volume = std::abs(dotProduct(V_A - V_D, crossProduct(V_B - V_D, V_C - V_D))) / 6.;

    return volume;
}

/*!
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the voxel is evaluated as the volume divied by
    the mean side area.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceVoxelInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double volume = evalVolume(vertexCoords);
    double area   = evalSurfaceArea(vertexCoords);

    double meanFaceArea = area / nFaces;

    double length = volume / meanFaceArea;

    return length;
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
    Evaluates the volume of an element with the specified vertex coordinates.

    The hexahedron is divided into six pyramids with a common vertex (the
    centroid of the hexahedron) and with one of the faces as the base of each.

    \param vertexCoords are the coordinate of the vertices
    \result The volume of the element.
*/
double ReferenceHexahedronInfo::evalVolume(const std::array<double, 3> *vertexCoords) const
{
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> pyramidVertexCoords;

    const ReferenceQuadInfo    &quadInfo    = static_cast<const ReferenceQuadInfo &>(getInfo(ElementType::QUAD));
    const ReferencePyramidInfo &pyramidInfo = static_cast<const ReferencePyramidInfo &>(getInfo(ElementType::PYRAMID));

    // The centroid is the apex of all the pyramids
    pyramidVertexCoords[4] = vertexCoords[0];
    for (int i = 1; i < nVertices; ++i) {
        pyramidVertexCoords[4] += vertexCoords[i];
    }
    pyramidVertexCoords[4] = pyramidVertexCoords[4] / double(nVertices);

    // Sum the volume of the pyramids having a face as the base and the
    // centroid as the apex.
    double volume = 0.;
    for (int i = 0; i < nFaces; ++i) {
        for (int n = 0; n < quadInfo.nVertices; ++n) {
            pyramidVertexCoords[quadInfo.nVertices - n - 1] = vertexCoords[faceConnect[i][n]];
        }

        volume += pyramidInfo.evalVolume(pyramidVertexCoords.data());
    }

    return volume;
}

/*!
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the hexahedron is evaluated as the volume
    divied by the mean side area.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceHexahedronInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double volume = evalVolume(vertexCoords);
    double area   = evalSurfaceArea(vertexCoords);

    double meanFaceArea = area / nFaces;

    double length = volume / meanFaceArea;

    return length;
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
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the pyramid is the height of the regular
    pyramid with a square base having the same inscribed sphere.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferencePyramidInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double volume = evalVolume(vertexCoords);
    double area   = evalSurfaceArea(vertexCoords);

    double inscribedRadius = 3 * volume / area;

    double length = inscribedRadius / (sqrt(5.) + 1.);

    return length;
}

/*!
    Evaluates the volume of an element with the specified vertex coordinates.

    The formula used to evaluate the volume is taken from "Calculation of the
    volume of a general hexahedron for flow predictions", D. E. Davies;
    D. J. Dalmond, AIAA Journal, June, Vol. 23, No. 6 : pp. 954-956.

    \param vertexCoords are the coordinate of the vertices
    \result The volume of the element.
*/
double ReferencePyramidInfo::evalVolume(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> RA = vertexCoords[0] - vertexCoords[4];
    const std::array<double, 3> DB = vertexCoords[3] - vertexCoords[1];
    const std::array<double, 3> AC = vertexCoords[2] - vertexCoords[0];
    const std::array<double, 3> AD = vertexCoords[1] - vertexCoords[0];
    const std::array<double, 3> AB = vertexCoords[3] - vertexCoords[0];

    double volume = dotProduct(RA, crossProduct(DB, AC)) / 6. + dotProduct(AC, crossProduct(AD, AB)) / 12.;

    return volume;
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

    faceConnect[0][0] = 0;
    faceConnect[0][1] = 1;
    faceConnect[0][2] = 2;

    faceConnect[1][0] = 4;
    faceConnect[1][1] = 3;
    faceConnect[1][2] = 5;

    faceConnect[2][0] = 4;
    faceConnect[2][1] = 1;
    faceConnect[2][2] = 0;
    faceConnect[2][3] = 3;

    faceConnect[3][0] = 1;
    faceConnect[3][1] = 4;
    faceConnect[3][2] = 5;
    faceConnect[3][3] = 2;

    faceConnect[4][0] = 3;
    faceConnect[4][1] = 0;
    faceConnect[4][2] = 2;
    faceConnect[4][3] = 5;

    initializeFaceEdges(facesInfo);
}

/*!
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the wedge is evaluated as the volume divied by
    the mean side area.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceWedgeInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double volume = evalVolume(vertexCoords);
    double area   = evalSurfaceArea(vertexCoords);

    double meanFaceArea = area / nFaces;

    double length = volume / meanFaceArea;

    return length;
}

/*!
    Evaluates the volume of an element with the specified vertex coordinates.

    The wedge is divided into three pyramids and two tetrahedron with a common
    vertex (the centroid of the hexahedron) and with one of the faces as the
    base of each.

    \param vertexCoords are the coordinate of the vertices
    \result The volume of the element.
*/
double ReferenceWedgeInfo::evalVolume(const std::array<double, 3> *vertexCoords) const
{
    std::array<int, 2> triaFaces = {{0, 1}};
    std::array<int, 3> quadFaces = {{2, 3, 4}};

    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> pyramidVertexCoords;
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> tetraVertexCoords;

    const ReferenceTriangleInfo &triaInfo    = static_cast<const ReferenceTriangleInfo &>(getInfo(ElementType::TRIANGLE));
    const ReferenceQuadInfo     &quadInfo    = static_cast<const ReferenceQuadInfo &>(getInfo(ElementType::QUAD));
    const ReferenceTetraInfo    &tetraInfo   = static_cast<const ReferenceTetraInfo &>(getInfo(ElementType::TETRA));
    const ReferencePyramidInfo  &pyramidInfo = static_cast<const ReferencePyramidInfo &>(getInfo(ElementType::PYRAMID));

    // The centroid is the apex of all the pyramids
    pyramidVertexCoords[4] = vertexCoords[0];
    for (int i = 1; i < nVertices; ++i) {
        pyramidVertexCoords[4] += vertexCoords[i];
    }
    pyramidVertexCoords[4] = pyramidVertexCoords[4] / double(nVertices);

    // The centroid is the apex of all the tetras
    tetraVertexCoords[3] = pyramidVertexCoords[4];

    // Sum the volume of the pyramids/tetra having a face as the base and the
    // centroid as the apex.
    double volume = 0.;
    for (int i : triaFaces) {
        for (int n = 0; n < triaInfo.nVertices; ++n) {
            tetraVertexCoords[triaInfo.nVertices - n - 1] = vertexCoords[faceConnect[i][n]];
        }

        volume += tetraInfo.evalVolume(tetraVertexCoords.data());
    }

    for (int i : quadFaces) {
        for (int n = 0; n < quadInfo.nVertices; ++n) {
            pyramidVertexCoords[quadInfo.nVertices - n - 1] = vertexCoords[faceConnect[i][n]];
        }

        volume += pyramidInfo.evalVolume(pyramidVertexCoords.data());
    }

    return volume;
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
    Evaluates the perimeter of a two-dimensional element with the specified
    vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The perimeter of the element.
*/
double Reference2DElementInfo::evalPerimeter(const std::array<double, 3> *vertexCoords) const
{
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> sideVertexCoords;
    const ReferenceLineInfo &sideInfo = static_cast<const ReferenceLineInfo &>(getInfo(ElementType::LINE));

    double perimeter = 0;
    for (int i = 0; i < nFaces; ++i) {
        for (int n = 0; n < sideInfo.nVertices; ++n) {
            sideVertexCoords[n] = vertexCoords[faceConnect[i][n]];
        }

        perimeter += sideInfo.evalLength(sideVertexCoords.data());
    }

    return perimeter;
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
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the triangle is the height of the regular
    triangle having the same inscribed circle.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceTriangleInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double area      = evalArea(vertexCoords);
    double perimeter = evalPerimeter(vertexCoords);

    double inscribedRadius = 2 * area / perimeter;

    double length = 3. * inscribedRadius;

    return length;
}

/*!
    Evaluates the area of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The area of the element.
*/
double ReferenceTriangleInfo::evalArea(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];

    double area = 0.5 * norm2(crossProduct(V_B - V_A, V_C - V_A));

    return area;
}

/*!
    Evaluates the normal of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \param point are the element reference coordinates of the point where the
    normal should be evaluated
    \result The normal of the element.
*/
std::array<double, 3> ReferenceTriangleInfo::evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];

    std::array<double, 3> normal = crossProduct(V_B - V_A, V_C - V_A);
    normal = normal / norm2(normal);

    return normal;
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
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the pixel is evaluated as the area divied by
    the mean side length.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferencePixelInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double area      = evalArea(vertexCoords);
    double perimeter = evalPerimeter(vertexCoords);

    double meanSideLength = perimeter / nFaces;

    double length = area / meanSideLength;

    return length;
}

/*!
    Evaluates the area of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The area of the element.
*/
double ReferencePixelInfo::evalArea(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];

    double area = norm2(crossProduct(V_B - V_A, V_C - V_A));

    return area;
}

/*!
    Evaluates the normal of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \param point are the element reference coordinates of the point where the
    normal should be evaluated
    \result The normal of the element.
*/
std::array<double, 3> ReferencePixelInfo::evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];

    std::array<double, 3> normal = crossProduct(V_B - V_A, V_C - V_A);
    normal = normal / norm2(normal);

    return normal;
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
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    The characteristics size of the quadrangle is evaluated as the area divied
    by the mean side length.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceQuadInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    double area      = evalArea(vertexCoords);
    double perimeter = evalPerimeter(vertexCoords);

    double meanSideLength = perimeter / nFaces;

    double length = area / meanSideLength;

    return length;
}

/*!
    Evaluates the area of an element with the specified vertex coordinates.

    NOTE: the formula used to evaluate the area assumes that the quadrilateral
    lies on a plane.

    \param vertexCoords are the coordinate of the vertices
    \result The area of the element.
*/
double ReferenceQuadInfo::evalArea(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];
    const std::array<double, 3> &V_C = vertexCoords[2];
    const std::array<double, 3> &V_D = vertexCoords[3];

    std::array<double, 3> mapping = crossProduct(V_B - V_A, V_D - V_A);
    mapping += 0.5 * crossProduct(V_B - V_A, V_C - V_D);
    mapping += 0.5 * crossProduct(V_C - V_B, V_D - V_A);

    double area = norm2(mapping);

    return area;
}

/*!
    Evaluates the normal of an element with the specified vertex coordinates.

    The formula used to evaluate the normal is taken from "Calculation of the
    volume of a general hexahedron for flow predictions", D. E. Davies;
    D. J. Dalmond, AIAA Journal, June, Vol. 23, No. 6 : pp. 954-956.

    \param vertexCoords are the coordinate of the vertices
    \param point are the element reference coordinates of the point where the
    normal should be evaluated
    \result The normal of the element.
*/
std::array<double, 3> ReferenceQuadInfo::evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point) const
{
    const std::array<double, 3> AB = vertexCoords[3] - vertexCoords[0];
    const std::array<double, 3> AD = vertexCoords[1] - vertexCoords[0];
    const std::array<double, 3> BC = vertexCoords[2] - vertexCoords[3];
    const std::array<double, 3> DC = vertexCoords[2] - vertexCoords[1];

    const double csi = point[0];
    const double eta = point[1];

    std::array<double, 3> normal = crossProduct(AB, AD);
    normal += csi * crossProduct(AB, DC);
    normal += eta * crossProduct(BC, AD);
    normal = normal / (- norm2(normal));

    return normal;
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
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the element.
*/
double ReferenceLineInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    return evalLength(vertexCoords);
}

/*!
    Evaluates the length of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the line.
*/
double ReferenceLineInfo::evalLength(const std::array<double, 3> *vertexCoords) const
{
    const std::array<double, 3> &V_A = vertexCoords[0];
    const std::array<double, 3> &V_B = vertexCoords[1];

    double length = norm2(V_B - V_A);

    return length;
}

/*!
    Evaluates the normal of an element with the specified vertex coordinates.

    \param vertexCoords are the coordinate of the vertices
    \param orientation is a vector perpendicular to the plane where the normal
    should lie
    \param point are the element reference coordinates of the point where the
    normal should be evaluated
    \result The normal of the element.
*/
std::array<double, 3> ReferenceLineInfo::evalNormal(const std::array<double, 3> *vertexCoords,
                                                    const std::array<double, 3> &orientation,
                                                    const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

    const std::array<double, 3> tangent = vertexCoords[1] - vertexCoords[0];

    std::array<double, 3> normal = crossProduct(tangent, orientation);
    normal = normal / norm2(normal);

    return normal;
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

/*!
    Evaluates the characteristics size of an element with the specified vertex
    coordinates.

    \param vertexCoords are the coordinate of the vertices
    \result The length of the element.
*/
double ReferenceVertexInfo::evalSize(const std::array<double, 3> *vertexCoords) const
{
    BITPIT_UNUSED(vertexCoords);

    return 0.;
}

/*!
    Evaluates the normal of an element with the specified vertex coordinates.

    For a vertex it makes no sense to evaluate a normal, this is only a dummy
    function that returns the direction that receives in input.

    \param vertexCoords are the coordinate of the vertices
    \param orientation is a vector that specifies the direction of the normal
    \result The normal of the element.
*/
std::array<double, 3> ReferenceVertexInfo::evalNormal(const std::array<double, 3> *vertexCoords,
                                                      const std::array<double, 3> &orientation) const
{
    BITPIT_UNUSED(vertexCoords);

    std::array<double, 3> normal = orientation;
    normal = normal / norm2(normal);

    return normal;
}

}
