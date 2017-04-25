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

#ifndef __BITPIT_ELEMENT_REFERENCE_HPP__
#define __BITPIT_ELEMENT_REFERENCE_HPP__

#include <assert.h>
#include <vector>

#include "bitpit_common.hpp"

#include "element_type.hpp"

namespace bitpit {

class ReferenceElementInfo {

public:
    ElementType type;
    int dimension;

    int nVertices;
    int nEdges;
    int nFaces;

    static const ReferenceElementInfo undefinedInfo;
    static const ReferenceElementInfo vertexInfo;
    static const ReferenceElementInfo lineInfo;
    static const ReferenceElementInfo triangleInfo;
    static const ReferenceElementInfo pixelInfo;
    static const ReferenceElementInfo quadInfo;
    static const ReferenceElementInfo tetraInfo;
    static const ReferenceElementInfo voxelInfo;
    static const ReferenceElementInfo hexahedronInfo;
    static const ReferenceElementInfo pyramidInfo;
    static const ReferenceElementInfo wedgeInfo;

    std::vector<ElementType> face_type;
    std::vector<std::vector<int>> faceConnect;
    std::vector<std::vector<int>> faceEdges;

    std::vector<ElementType> edge_type;
    std::vector<std::vector<int>> edgeConnect;

    ReferenceElementInfo();
    ReferenceElementInfo(ElementType type);

    static const ReferenceElementInfo & getInfo(ElementType type);

private:
    void initializeUndefinedInfo();
    void initializeVertexInfo();
    void initializeLineInfo();
    void initializeTriangleInfo();
    void initializePixelInfo();
    void initializeQuadInfo();
    void initializeTetraInfo();
    void initializeVoxelInfo();
    void initializeHexahedronInfo();
    void initializePyramidInfo();
    void initializeWedgeInfo();

    void initializeFaceEdges(const std::vector<ReferenceElementInfo *> &facesInfo = std::vector<ReferenceElementInfo *>());
};

}

#endif
