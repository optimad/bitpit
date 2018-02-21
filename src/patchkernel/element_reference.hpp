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
    static const int MAX_ELEM_VERTICES = 8;
    static const int MAX_ELEM_FACES    = 6;
    static const int MAX_ELEM_EDGES    = 12;

    static bool hasInfo(ElementType type);
    static const ReferenceElementInfo & getInfo(ElementType type);

    int dimension;
    ElementType type;

    int nVertices;

    int nFaces;
    std::vector<ElementType> face_type;
    std::vector<std::vector<int>> faceConnect;
    std::vector<std::vector<int>> faceEdges;

    int nEdges;
    std::vector<ElementType> edge_type;
    std::vector<std::vector<int>> edgeConnect;

    virtual ~ReferenceElementInfo();

    virtual double evalSize(const std::array<double, 3> *vertexCoords) const = 0;

    virtual double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const = 0;

protected:
    ReferenceElementInfo(int _dimension, ElementType _type, int _nVertices, int _nFaces, int _nEdges);

    ReferenceElementInfo(ReferenceElementInfo const &) = delete;
    ReferenceElementInfo & operator=(ReferenceElementInfo const &) = delete;

    void initializeFaceEdges(const std::vector<const ReferenceElementInfo *> &facesInfo);

};

class Reference3DElementInfo : public ReferenceElementInfo {

public:
    virtual double evalVolume(const std::array<double, 3> *vertexCoords) const = 0;
    double evalSurfaceArea(const std::array<double, 3> *vertexCoords) const;
    double evalEdgePerimeter(const std::array<double, 3> *vertexCoords) const;

    double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const override;

protected:
    Reference3DElementInfo(ElementType type, int nVertices, int nFaces);

};

class ReferenceTetraInfo : public Reference3DElementInfo {

friend class ReferenceElementInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalVolume(const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceTetraInfo info;

protected:
    ReferenceTetraInfo();

    ReferenceTetraInfo(ReferenceTetraInfo const &) = delete;
    ReferenceTetraInfo & operator=(ReferenceTetraInfo const &) = delete;

};

class ReferenceVoxelInfo : public Reference3DElementInfo {

friend class ReferenceElementInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalVolume(const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceVoxelInfo info;

protected:
    ReferenceVoxelInfo();

    ReferenceVoxelInfo(ReferenceVoxelInfo const &) = delete;
    ReferenceVoxelInfo & operator=(ReferenceVoxelInfo const &) = delete;

};

class ReferenceHexahedronInfo : public Reference3DElementInfo {

friend class ReferenceElementInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalVolume(const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceHexahedronInfo info;

protected:
    ReferenceHexahedronInfo();

    ReferenceHexahedronInfo(ReferenceHexahedronInfo const &) = delete;
    ReferenceHexahedronInfo & operator=(ReferenceHexahedronInfo const &) = delete;

};

class ReferencePyramidInfo : public Reference3DElementInfo {

friend class ReferenceElementInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalVolume(const std::array<double, 3> *vertexCoords) const override;

    const static ReferencePyramidInfo info;

protected:
    ReferencePyramidInfo();

    ReferencePyramidInfo(ReferencePyramidInfo const &) = delete;
    ReferencePyramidInfo & operator=(ReferencePyramidInfo const &) = delete;

};

class ReferenceWedgeInfo : public Reference3DElementInfo {

friend class ReferenceElementInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalVolume(const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceWedgeInfo info;

protected:
    ReferenceWedgeInfo();

    ReferenceWedgeInfo(ReferenceWedgeInfo const &) = delete;
    ReferenceWedgeInfo & operator=(ReferenceWedgeInfo const &) = delete;

};

class Reference2DElementInfo : public ReferenceElementInfo {

public:
    virtual double evalArea(const std::array<double, 3> *vertexCoords) const = 0;
    double evalPerimeter(const std::array<double, 3> *vertexCoords) const;

    virtual std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const = 0;

    double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const override;

protected:
    Reference2DElementInfo(ElementType type, int nVertices);

};

class ReferenceTriangleInfo : public Reference2DElementInfo {

friend class ReferenceElementInfo;
friend class ReferenceTetraInfo;
friend class ReferencePyramidInfo;
friend class ReferenceWedgeInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalArea(const std::array<double, 3> *vertexCoords) const override;

    std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const override;

    double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceTriangleInfo info;

protected:
    ReferenceTriangleInfo();

    ReferenceTriangleInfo(ReferenceTriangleInfo const &) = delete;
    ReferenceTriangleInfo & operator=(ReferenceTriangleInfo const &) = delete;

};

class ReferencePixelInfo : public Reference2DElementInfo {

friend class ReferenceElementInfo;
friend class ReferenceVoxelInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalArea(const std::array<double, 3> *vertexCoords) const override;

    std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const override;

    const static ReferencePixelInfo info;

protected:
    ReferencePixelInfo();

    ReferencePixelInfo(ReferencePixelInfo const &) = delete;
    ReferencePixelInfo & operator=(ReferencePixelInfo const &) = delete;

};

class ReferenceQuadInfo : public Reference2DElementInfo {

friend class ReferenceElementInfo;
friend class ReferenceHexahedronInfo;
friend class ReferencePyramidInfo;
friend class ReferenceWedgeInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalArea(const std::array<double, 3> *vertexCoords) const override;

    std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const override;

    const static ReferenceQuadInfo info;

protected:
    ReferenceQuadInfo();

    ReferenceQuadInfo(ReferenceQuadInfo const &) = delete;
    ReferenceQuadInfo & operator=(ReferenceQuadInfo const &) = delete;

};

class Reference1DElementInfo : public ReferenceElementInfo {

public:
    virtual double evalLength(const std::array<double, 3> *vertexCoords) const = 0;

    virtual std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &orientation = {{0., 0., 1.}}, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const = 0;

protected:
    Reference1DElementInfo(ElementType type);

};

class ReferenceLineInfo : public Reference1DElementInfo {

friend class ReferenceElementInfo;
friend class ReferenceTriangleInfo;
friend class ReferencePixelInfo;
friend class ReferenceQuadInfo;
friend class ReferenceTetraInfo;
friend class ReferenceVoxelInfo;
friend class ReferenceHexahedronInfo;
friend class ReferencePyramidInfo;
friend class ReferenceWedgeInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    double evalLength(const std::array<double, 3> *vertexCoords) const override;

    std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &orientation = {{0., 0., 1.}}, const std::array<double, 3> &point = {{0.5, 0.5, 0.5}}) const override;

    double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceLineInfo info;

protected:
    ReferenceLineInfo();

    ReferenceLineInfo(ReferenceLineInfo const &) = delete;
    ReferenceLineInfo & operator=(ReferenceLineInfo const &) = delete;

};

class Reference0DElementInfo : public ReferenceElementInfo {

public:
    virtual std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &orientation = {{1., 0., 0.}}) const = 0;

protected:
    Reference0DElementInfo(ElementType type);

};

class ReferenceVertexInfo : public Reference0DElementInfo {

friend class ReferenceElementInfo;
friend class ReferenceLineInfo;
friend class ReferenceTriangleInfo;
friend class ReferencePixelInfo;
friend class ReferenceQuadInfo;

public:
    double evalSize(const std::array<double, 3> *vertexCoords) const override;

    std::array<double, 3> evalNormal(const std::array<double, 3> *vertexCoords, const std::array<double, 3> &orientation = {{1., 0., 0.}}) const override;

    double evalPointDistance(const std::array<double, 3> &point, const std::array<double, 3> *vertexCoords) const override;

    const static ReferenceVertexInfo info;

protected:
    ReferenceVertexInfo();

    ReferenceVertexInfo(ReferenceVertexInfo const &) = delete;
    ReferenceVertexInfo & operator=(ReferenceVertexInfo const &) = delete;

};

}

#endif
