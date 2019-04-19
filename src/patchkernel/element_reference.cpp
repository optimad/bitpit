/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include <set>

#include "bitpit_CG.hpp"
#include "bitpit_containers.hpp"
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
      nVertices(_nVertices), nFaces(_nFaces), nEdges(_nEdges)
{
    faceTypeStorage.fill(ElementType::UNDEFINED);
    for (int i = 0; i < MAX_ELEM_FACES; ++i) {
        faceConnectStorage[i].fill(-1);
        faceEdgeStorage[i].fill(-1);
    }

    edgeTypeStorage.fill(ElementType::UNDEFINED);
    for (int i = 0; i < MAX_ELEM_FACES; ++i) {
        edgeConnectStorage[i].fill(-1);
    }
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
    switch (type) {

    case (ElementType::VERTEX):
        return ReferenceVertexInfo::info;

    case (ElementType::LINE):
        return ReferenceLineInfo::info;

    case (ElementType::TRIANGLE):
        return ReferenceTriangleInfo::info;

    case (ElementType::PIXEL):
        return ReferencePixelInfo::info;

    case (ElementType::QUAD):
        return ReferenceQuadInfo::info;

    case (ElementType::TETRA):
        return ReferenceTetraInfo::info;

    case (ElementType::VOXEL):
        return ReferenceVoxelInfo::info;

    case (ElementType::HEXAHEDRON):
        return ReferenceHexahedronInfo::info;

    case (ElementType::PYRAMID):
        return ReferencePyramidInfo::info;

    case (ElementType::WEDGE):
        return ReferenceWedgeInfo::info;

    default:
        BITPIT_UNREACHABLE("Unsupported element");
        throw std::runtime_error("Unsupported element");

    }
}

/*!
    Initializes the list of edges associated to the faces
*/
void ReferenceElementInfo::initializeFaceEdges(const std::vector<const ReferenceElementInfo *> &facesInfo,
                                               const std::vector<const ReferenceElementInfo *> &edgesInfo)
{
    for (int k = 0; k < nFaces; ++k) {
        const ReferenceElementInfo &faceInfo = *(facesInfo[k]);

        int nFaceEdges = faceInfo.nFaces;
        int faceEdgeCounter = 0;
        for (int i = 0; i < nFaceEdges; ++i) {
            const ReferenceElementInfo &faceEdgeInfo = *(edgesInfo[i]);

            // Connectivity of the edge associated to the face
            const int *localFaceEdgeConnect = faceInfo.faceConnectStorage[i].data();

            std::set<int> faceEdgeConnect;
            for (int n = 0; n < faceEdgeInfo.nVertices; ++n) {
                int localVertexId = localFaceEdgeConnect[n];
                int vertexId      = faceConnectStorage[k][localVertexId];

                faceEdgeConnect.insert(vertexId);
            }

            // Search the edge that has the same connectivity of the face edge
            for (int j = 0; j < nEdges; ++j) {
                const ReferenceElementInfo &guessEdgeInfo = *(edgesInfo[j]);

                // If face edge and the guess edge have a different type, the
                // two edge cannot be the same.
                if (guessEdgeInfo.type != faceEdgeInfo.type) {
                    continue;
                }

                // If the connecitivity of the face edge and the one of the
                // guess edge are the same, the two edges coincides.
                const std::set<int> guessEdgeConnect = std::set<int>(edgeConnectStorage[j].begin(), edgeConnectStorage[j].begin() + guessEdgeInfo.nVertices);
                if (faceEdgeConnect == guessEdgeConnect) {
                    faceEdgeStorage[k][faceEdgeCounter] = j;
                    ++faceEdgeCounter;
                }
            }
        }

        assert(faceEdgeCounter == nFaceEdges);
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
        ElementType faceType = faceTypeStorage[i];
        const Reference2DElementInfo &faceInfo = static_cast<const Reference2DElementInfo &>(getInfo(faceType));
        for (int n = 0; n < faceInfo.nVertices; ++n) {
            faceVertexCoords[n] = vertexCoords[faceConnectStorage[i][n]];
        }

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
            edgeVertexCoords[n] = vertexCoords[edgeConnectStorage[i][n]];
        }

        perimeter += edgeInfo.evalLength(edgeVertexCoords.data());
    }

    return perimeter;
}

/*!
    Evaluates the distance between the element and the specified point.

    \param[in] point is the point
    \param vertexCoords are the coordinate of the vertices
    \result The distance between the element and the specified point.
*/
double Reference3DElementInfo::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const
{
    std::array<std::array<double, 3>, MAX_ELEM_VERTICES> faceVertexCoords;

    double distance = std::numeric_limits<double>::max();
    for (int i = 0; i < nFaces; ++i) {
        ElementType faceType = faceTypeStorage[i];
        const Reference2DElementInfo &faceInfo = static_cast<const Reference2DElementInfo &>(getInfo(faceType));
        for (int n = 0; n < faceInfo.nVertices; ++n) {
            faceVertexCoords[n] = vertexCoords[faceConnectStorage[i][n]];
        }

        distance = std::min(faceInfo.evalPointDistance(point, faceVertexCoords.data()), distance);
    }

    return distance;
}

/*!
    \class ReferenceTetraInfo
    \ingroup patchelements

    \brief The ReferenceTetraInfo class defines the information about the
    reference tetrahedron.
*/

const ReferenceTetraInfo ReferenceTetraInfo::info;

/*!
    Default constructor
*/
ReferenceTetraInfo::ReferenceTetraInfo()
    : Reference3DElementInfo(ElementType::TETRA, 4, 4)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceTriangleInfo triangleInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &lineInfo;

        edgeTypeStorage[k] = LINE;
    }

    edgeConnectStorage[0][0] = 0;
    edgeConnectStorage[0][1] = 1;

    edgeConnectStorage[1][0] = 1;
    edgeConnectStorage[1][1] = 2;

    edgeConnectStorage[2][0] = 2;
    edgeConnectStorage[2][1] = 0;

    edgeConnectStorage[3][0] = 3;
    edgeConnectStorage[3][1] = 0;

    edgeConnectStorage[4][0] = 3;
    edgeConnectStorage[4][1] = 1;

    edgeConnectStorage[5][0] = 3;
    edgeConnectStorage[5][1] = 2;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &triangleInfo;

        faceTypeStorage[k] = TRIANGLE;
    }

    faceConnectStorage[0][0] = 1;
    faceConnectStorage[0][1] = 0;
    faceConnectStorage[0][2] = 2;

    faceConnectStorage[1][0] = 0;
    faceConnectStorage[1][1] = 3;
    faceConnectStorage[1][2] = 2;

    faceConnectStorage[2][0] = 3;
    faceConnectStorage[2][1] = 1;
    faceConnectStorage[2][2] = 2;

    faceConnectStorage[3][0] = 0;
    faceConnectStorage[3][1] = 1;
    faceConnectStorage[3][2] = 3;

    initializeFaceEdges(facesInfo, edgesInfo);
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

const ReferenceVoxelInfo ReferenceVoxelInfo::info;

/*!
    Default constructor
*/
ReferenceVoxelInfo::ReferenceVoxelInfo()
    : Reference3DElementInfo(ElementType::VOXEL, 8, 6)
{
    const ReferenceLineInfo lineInfo;
    const ReferencePixelInfo pixelInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &lineInfo;

        edgeTypeStorage[k] = LINE;
    }

    edgeConnectStorage[0][0] = 0;
    edgeConnectStorage[0][1] = 2;

    edgeConnectStorage[1][0] = 1;
    edgeConnectStorage[1][1] = 3;

    edgeConnectStorage[2][0] = 0;
    edgeConnectStorage[2][1] = 1;

    edgeConnectStorage[3][0] = 2;
    edgeConnectStorage[3][1] = 3;

    edgeConnectStorage[4][0] = 0;
    edgeConnectStorage[4][1] = 4;

    edgeConnectStorage[5][0] = 1;
    edgeConnectStorage[5][1] = 5;

    edgeConnectStorage[6][0] = 2;
    edgeConnectStorage[6][1] = 6;

    edgeConnectStorage[7][0] = 3;
    edgeConnectStorage[7][1] = 7;

    edgeConnectStorage[8][0] = 4;
    edgeConnectStorage[8][1] = 6;

    edgeConnectStorage[9][0] = 5;
    edgeConnectStorage[9][1] = 7;

    edgeConnectStorage[10][0] = 4;
    edgeConnectStorage[10][1] = 5;

    edgeConnectStorage[11][0] = 6;
    edgeConnectStorage[11][1] = 7;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &pixelInfo;

        faceTypeStorage[k] = PIXEL;
    }

    faceConnectStorage[0][0] = 2;
    faceConnectStorage[0][1] = 0;
    faceConnectStorage[0][2] = 6;
    faceConnectStorage[0][3] = 4;

    faceConnectStorage[1][0] = 1;
    faceConnectStorage[1][1] = 3;
    faceConnectStorage[1][2] = 5;
    faceConnectStorage[1][3] = 7;

    faceConnectStorage[2][0] = 0;
    faceConnectStorage[2][1] = 1;
    faceConnectStorage[2][2] = 4;
    faceConnectStorage[2][3] = 5;

    faceConnectStorage[3][0] = 3;
    faceConnectStorage[3][1] = 2;
    faceConnectStorage[3][2] = 7;
    faceConnectStorage[3][3] = 6;

    faceConnectStorage[4][0] = 2;
    faceConnectStorage[4][1] = 3;
    faceConnectStorage[4][2] = 0;
    faceConnectStorage[4][3] = 1;

    faceConnectStorage[5][0] = 4;
    faceConnectStorage[5][1] = 5;
    faceConnectStorage[5][2] = 6;
    faceConnectStorage[5][3] = 7;

    initializeFaceEdges(facesInfo, edgesInfo);
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

const ReferenceHexahedronInfo ReferenceHexahedronInfo::info;

/*!
    Default constructor
*/
ReferenceHexahedronInfo::ReferenceHexahedronInfo()
    : Reference3DElementInfo(ElementType::HEXAHEDRON, 8, 6)
{
    const ReferenceLineInfo lineInfo;
    const ReferenceQuadInfo quadInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &lineInfo;

        edgeTypeStorage[k] = LINE;
    }

    edgeConnectStorage[0][0] = 1;
    edgeConnectStorage[0][1] = 0;

    edgeConnectStorage[1][0] = 1;
    edgeConnectStorage[1][1] = 2;

    edgeConnectStorage[2][0] = 2;
    edgeConnectStorage[2][1] = 3;

    edgeConnectStorage[3][0] = 3;
    edgeConnectStorage[3][1] = 0;

    edgeConnectStorage[4][0] = 4;
    edgeConnectStorage[4][1] = 5;

    edgeConnectStorage[5][0] = 5;
    edgeConnectStorage[5][1] = 6;

    edgeConnectStorage[6][0] = 6;
    edgeConnectStorage[6][1] = 7;

    edgeConnectStorage[7][0] = 7;
    edgeConnectStorage[7][1] = 4;

    edgeConnectStorage[8][0] = 0;
    edgeConnectStorage[8][1] = 4;

    edgeConnectStorage[9][0] = 1;
    edgeConnectStorage[9][1] = 5;

    edgeConnectStorage[10][0] = 2;
    edgeConnectStorage[10][1] = 6;

    edgeConnectStorage[11][0] = 3;
    edgeConnectStorage[11][1] = 7;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &quadInfo;

        faceTypeStorage[k] = QUAD;
    }

    faceConnectStorage[0][0] = 1;
    faceConnectStorage[0][1] = 0;
    faceConnectStorage[0][2] = 3;
    faceConnectStorage[0][3] = 2;

    faceConnectStorage[1][0] = 4;
    faceConnectStorage[1][1] = 5;
    faceConnectStorage[1][2] = 6;
    faceConnectStorage[1][3] = 7;

    faceConnectStorage[2][0] = 7;
    faceConnectStorage[2][1] = 3;
    faceConnectStorage[2][2] = 0;
    faceConnectStorage[2][3] = 4;

    faceConnectStorage[3][0] = 5;
    faceConnectStorage[3][1] = 1;
    faceConnectStorage[3][2] = 2;
    faceConnectStorage[3][3] = 6;

    faceConnectStorage[4][0] = 4;
    faceConnectStorage[4][1] = 0;
    faceConnectStorage[4][2] = 1;
    faceConnectStorage[4][3] = 5;

    faceConnectStorage[5][0] = 6;
    faceConnectStorage[5][1] = 2;
    faceConnectStorage[5][2] = 3;
    faceConnectStorage[5][3] = 7;

    initializeFaceEdges(facesInfo, edgesInfo);
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
            pyramidVertexCoords[quadInfo.nVertices - n - 1] = vertexCoords[faceConnectStorage[i][n]];
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

const ReferencePyramidInfo ReferencePyramidInfo::info;

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
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &lineInfo;

        edgeTypeStorage[k] = LINE;
    }

    edgeConnectStorage[0][0] = 0;
    edgeConnectStorage[0][1] = 1;

    edgeConnectStorage[1][0] = 1;
    edgeConnectStorage[1][1] = 2;

    edgeConnectStorage[2][0] = 2;
    edgeConnectStorage[2][1] = 3;

    edgeConnectStorage[3][0] = 3;
    edgeConnectStorage[3][1] = 0;

    edgeConnectStorage[4][0] = 4;
    edgeConnectStorage[4][1] = 0;

    edgeConnectStorage[5][0] = 4;
    edgeConnectStorage[5][1] = 1;

    edgeConnectStorage[6][0] = 4;
    edgeConnectStorage[6][1] = 2;

    edgeConnectStorage[7][0] = 4;
    edgeConnectStorage[7][1] = 3;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        if (k == 0) {
            facesInfo[k] = &quadInfo;

            faceTypeStorage[k] = ElementType::QUAD;
        } else {
            facesInfo[k] = &triangleInfo;

            faceTypeStorage[k] = ElementType::TRIANGLE;
        }
    }

    faceConnectStorage[0][0] = 0;
    faceConnectStorage[0][1] = 3;
    faceConnectStorage[0][2] = 2;
    faceConnectStorage[0][3] = 1;

    faceConnectStorage[1][0] = 3;
    faceConnectStorage[1][1] = 0;
    faceConnectStorage[1][2] = 4;

    faceConnectStorage[2][0] = 0;
    faceConnectStorage[2][1] = 1;
    faceConnectStorage[2][2] = 4;

    faceConnectStorage[3][0] = 1;
    faceConnectStorage[3][1] = 2;
    faceConnectStorage[3][2] = 4;

    faceConnectStorage[4][0] = 2;
    faceConnectStorage[4][1] = 3;
    faceConnectStorage[4][2] = 4;

    initializeFaceEdges(facesInfo, edgesInfo);
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

const ReferenceWedgeInfo ReferenceWedgeInfo::info;

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
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgeTypeStorage[k]   = LINE;
        edgesInfo[k]   = &lineInfo;
    }

    edgeConnectStorage[0][0] = 1;
    edgeConnectStorage[0][1] = 0;

    edgeConnectStorage[1][0] = 1;
    edgeConnectStorage[1][1] = 2;

    edgeConnectStorage[2][0] = 2;
    edgeConnectStorage[2][1] = 0;

    edgeConnectStorage[3][0] = 3;
    edgeConnectStorage[3][1] = 4;

    edgeConnectStorage[4][0] = 4;
    edgeConnectStorage[4][1] = 5;

    edgeConnectStorage[5][0] = 5;
    edgeConnectStorage[5][1] = 3;

    edgeConnectStorage[6][0] = 3;
    edgeConnectStorage[6][1] = 0;

    edgeConnectStorage[7][0] = 4;
    edgeConnectStorage[7][1] = 1;

    edgeConnectStorage[8][0] = 5;
    edgeConnectStorage[8][1] = 2;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        if (k == 0 || k == 1) {
            facesInfo[k] = &triangleInfo;

            faceTypeStorage[k] = TRIANGLE;
        } else {
            facesInfo[k] = &quadInfo;

            faceTypeStorage[k] = QUAD;
        }
    }

    faceConnectStorage[0][0] = 0;
    faceConnectStorage[0][1] = 1;
    faceConnectStorage[0][2] = 2;

    faceConnectStorage[1][0] = 4;
    faceConnectStorage[1][1] = 3;
    faceConnectStorage[1][2] = 5;

    faceConnectStorage[2][0] = 4;
    faceConnectStorage[2][1] = 1;
    faceConnectStorage[2][2] = 0;
    faceConnectStorage[2][3] = 3;

    faceConnectStorage[3][0] = 1;
    faceConnectStorage[3][1] = 4;
    faceConnectStorage[3][2] = 5;
    faceConnectStorage[3][3] = 2;

    faceConnectStorage[4][0] = 3;
    faceConnectStorage[4][1] = 0;
    faceConnectStorage[4][2] = 2;
    faceConnectStorage[4][3] = 5;

    initializeFaceEdges(facesInfo, edgesInfo);
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
            tetraVertexCoords[triaInfo.nVertices - n - 1] = vertexCoords[faceConnectStorage[i][n]];
        }

        volume += tetraInfo.evalVolume(tetraVertexCoords.data());
    }

    for (int i : quadFaces) {
        for (int n = 0; n < quadInfo.nVertices; ++n) {
            pyramidVertexCoords[quadInfo.nVertices - n - 1] = vertexCoords[faceConnectStorage[i][n]];
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
            sideVertexCoords[n] = vertexCoords[faceConnectStorage[i][n]];
        }

        perimeter += sideInfo.evalLength(sideVertexCoords.data());
    }

    return perimeter;
}

/*!
    Evaluates the distance between the element and the specified point.

    \param[in] point is the point
    \param vertexCoords are the coordinate of the vertices
    \result The distance between the element and the specified point.
*/
double Reference2DElementInfo::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const
{
    return CGElem::distancePointPolygon(point, nVertices, vertexCoords);
}

/*!
    \class ReferenceTriangleInfo
    \ingroup patchelements

    \brief The ReferenceTriangleInfo class defines the information about the
    reference triangle.
*/

const ReferenceTriangleInfo ReferenceTriangleInfo::info;

/*!
    Default constructor
*/
ReferenceTriangleInfo::ReferenceTriangleInfo()
    : Reference2DElementInfo(ElementType::TRIANGLE, 3)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &vertexInfo;

        edgeTypeStorage[k]       = ElementType::VERTEX;
        edgeConnectStorage[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &lineInfo;

        faceTypeStorage[k]       = LINE;
        faceConnectStorage[k][0] = k;
        faceConnectStorage[k][1] = (k + 1) % nVertices;
    }

    initializeFaceEdges(facesInfo, edgesInfo);
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
    Evaluates the distance between the element and the specified point.

    \param[in] point is the point
    \param vertexCoords are the coordinate of the vertices
    \result The distance between the element and the specified point.
*/
double ReferenceTriangleInfo::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const
{
    return CGElem::distancePointTriangle(point, vertexCoords[0], vertexCoords[1], vertexCoords[2]);
}

/*!
    \class ReferencePixelInfo
    \ingroup patchelements

    \brief The ReferencePixelInfo class defines the information about the
    reference pixel.
*/

const ReferencePixelInfo ReferencePixelInfo::info;

/*!
    Default constructor
*/
ReferencePixelInfo::ReferencePixelInfo()
    : Reference2DElementInfo(ElementType::PIXEL, 4)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &vertexInfo;

        edgeTypeStorage[k]       = ElementType::VERTEX;
        edgeConnectStorage[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &lineInfo;

        faceTypeStorage[k] = LINE;
    }

    faceConnectStorage[0][0] = 2;
    faceConnectStorage[0][1] = 0;

    faceConnectStorage[1][0] = 1;
    faceConnectStorage[1][1] = 3;

    faceConnectStorage[2][0] = 0;
    faceConnectStorage[2][1] = 1;

    faceConnectStorage[3][0] = 3;
    faceConnectStorage[3][1] = 2;

    initializeFaceEdges(facesInfo, edgesInfo);
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

const ReferenceQuadInfo ReferenceQuadInfo::info;

/*!
    Default constructor
*/
ReferenceQuadInfo::ReferenceQuadInfo()
    : Reference2DElementInfo(ElementType::QUAD, 4)
{
    const ReferenceVertexInfo vertexInfo;
    const ReferenceLineInfo lineInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &vertexInfo;

        edgeTypeStorage[k]       = ElementType::VERTEX;
        edgeConnectStorage[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &lineInfo;

        faceTypeStorage[k]       = LINE;
        faceConnectStorage[k][0] = k;
        faceConnectStorage[k][1] = (k + 1) % nVertices;
    }

    initializeFaceEdges(facesInfo, edgesInfo);
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

const ReferenceLineInfo ReferenceLineInfo::info;

/*!
    Default constructor
*/
ReferenceLineInfo::ReferenceLineInfo()
    : Reference1DElementInfo(ElementType::LINE)
{
    const ReferenceVertexInfo vertexInfo;

    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    for (int k = 0; k < nEdges; ++k) {
        edgesInfo[k] = &vertexInfo;

        edgeTypeStorage[k]       = ElementType::VERTEX;
        edgeConnectStorage[k][0] = k;
    }

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    for (int k = 0; k < nFaces; ++k) {
        facesInfo[k] = &vertexInfo;

        faceTypeStorage[k]       = ElementType::VERTEX;
        faceConnectStorage[k][0] = k;
    }

    initializeFaceEdges(facesInfo, edgesInfo);
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
    Evaluates the distance between the element and the specified point.

    \param[in] point is the point
    \param vertexCoords are the coordinate of the vertices
    \result The distance between the element and the specified point.
*/
double ReferenceLineInfo::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const
{
    return CGElem::distancePointSegment(point, vertexCoords[0], vertexCoords[1]);
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

const ReferenceVertexInfo ReferenceVertexInfo::info;

/*!
    Default constructor
*/
ReferenceVertexInfo::ReferenceVertexInfo()
    : Reference0DElementInfo(ElementType::VERTEX)
{
    // Edge data
    std::vector<const ReferenceElementInfo *> edgesInfo(nEdges);

    edgesInfo[0] = this;

    edgeTypeStorage[0] = ElementType::VERTEX;

    edgeConnectStorage[0][0] = 0;

    // Face data
    std::vector<const ReferenceElementInfo *> facesInfo(nFaces);

    facesInfo[0] = this;

    faceTypeStorage[0] = ElementType::VERTEX;

    faceConnectStorage[0][0] = 0;

    initializeFaceEdges(facesInfo, edgesInfo);
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

/*!
    Evaluates the distance between the element and the specified point.

    \param[in] point is the point
    \param vertexCoords are the coordinate of the vertices
    \result The distance between the element and the specified point.
*/
double ReferenceVertexInfo::evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const
{
    return norm2(point - vertexCoords[0]);
}

}
