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

# ifndef __BITPIT_LEVELSET_SEGMENTATION_OBJECT_HPP__
# define __BITPIT_LEVELSET_SEGMENTATION_OBJECT_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_set>
# include <unordered_map>

#include "bitpit_CG.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_volcartesian.hpp"
#include "bitpit_voloctree.hpp"

#include "levelSetCartesianKernel.hpp"
#include "levelSetCommon.hpp"
#include "levelSetKernel.hpp"
#include "levelSetObject.hpp"

namespace bitpit{

class SurfaceSkdTree;

class LevelSetSegmentationKernel {

public:
    LevelSetSegmentationKernel();
    LevelSetSegmentationKernel(const LevelSetSegmentationKernel &other);
    LevelSetSegmentationKernel(LevelSetSegmentationKernel &&other);
    LevelSetSegmentationKernel(const SurfUnstructured *surface, double featureAngle);
    LevelSetSegmentationKernel(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle);

    virtual ~LevelSetSegmentationKernel() = default;

    const SurfUnstructured & getSurface() const;
    void setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle = 2. * BITPIT_PI);
    void setSurface(const SurfUnstructured *surface, double featureAngle = 2. * BITPIT_PI);

    double getFeatureAngle() const;

    const SurfaceSkdTree & getSearchTree() const;

    int evalLevelsetInfo(const std::array<double,3> &point, long segmentId, bool signedLevelSet,
                         double *value, std::array<double,3> *gradient) const;

    double getSegmentSize(long segmentId) const;
    double getMinSegmentSize() const;
    double getMaxSegmentSize() const;

private:
    typedef std::pair<long, int> SegmentVertexKey;

    const SurfUnstructured *m_surface;
    std::unique_ptr<const SurfUnstructured> m_ownedSurface;
    double m_featureAngle;

    std::unique_ptr<SurfaceSkdTree> m_searchTree;

    PiercedStorage<std::size_t> m_segmentVertexOffset;

    mutable PiercedStorage<bool> m_segmentNormalsValid;
    mutable PiercedStorage<std::array<double,3>> m_segmentNormalsStorage;
    mutable PiercedStorage<bool> m_unlimitedVertexNormalsValid;
    mutable PiercedStorage<std::array<double,3>> m_unlimitedVertexNormalsStorage;
    mutable std::vector<bool> m_limitedSegmentVertexNormalValid;
    mutable std::unordered_map<SegmentVertexKey, std::array<double,3>, utils::hashing::hash<SegmentVertexKey>> m_limitedSegmentVertexNormalStorage;

    std::array<double,3> computePseudoNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const;
    std::array<double,3> computeSurfaceNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const;

    std::array<double,3> computeSegmentNormal( const SurfUnstructured::CellConstIterator &segmentIterator ) const;
    std::array<double,3> computeSegmentEdgeNormal( const SurfUnstructured::CellConstIterator &segmentIterator, int edge ) const;
    std::array<double,3> computeSegmentVertexNormal( const SurfUnstructured::CellConstIterator &segmentIterator, int vertex, bool limited ) const;

};

class LevelSetSegmentationObject : public LevelSetObject, public LevelSetSegmentationKernel {

public:
    LevelSetSegmentationObject(int);
    LevelSetSegmentationObject(int, std::unique_ptr<const SurfUnstructured> &&, double featureAngle = 2. * BITPIT_PI);
    LevelSetSegmentationObject(int, const SurfUnstructured*, double featureAngle = 2. * BITPIT_PI);

    LevelSetSegmentationObject * clone() const override;

    void setFieldCache(LevelSetField field, LevelSetCacheMode cacheMode) override;

    LevelSetFieldset getSupportedFields() const override;

    double getMinSurfaceFeatureSize() const;
    double getMaxSurfaceFeatureSize() const;

    int evalCellPart(long id) const;
    std::array<double,3> evalCellNormal(long id, bool signedLevelSet) const;
    double evalCellSurfaceFeatureSize(long id) const;
    long evalCellSupport(long id) const;

    int evalPart(const std::array<double,3> &point) const;
    std::array<double,3> evalNormal(const std::array<double,3> &point, bool signedLevelSet) const;
    double evalSurfaceFeatureSize(const std::array<double,3> &point) const;
    long evalSupport(const std::array<double,3> &point) const;

    BITPIT_DEPRECATED(int getPart(long cellId) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3> getNormal(long cellId) const);
    BITPIT_DEPRECATED(long getSupport(long cellId) const);
    BITPIT_DEPRECATED(double getSurfaceFeatureSize(long cellId) const);

    BITPIT_DEPRECATED(int getPart(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3> getNormal(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(long getSupport(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(double getSurfaceFeatureSize(const std::array<double,3> &point) const);

protected:
    void fillNarrowBandCellCaches() override;
    void fillNarrowBandCellCaches(const std::vector<adaption::Info> &adaptionData) override;

    void fillFieldCellCache(long id, LevelSetField field, double searchRadius = std::numeric_limits<double>::max()) override;

    bool _isCellInNarrowBand(long id, bool checkNeighbours, double *maximumDistance = nullptr) const override;

    short _evalCellSign(long id) const override;
    double _evalCellValue(long id, bool signedLevelSet) const override;
    std::array<double,3> _evalCellGradient(long id, bool signedLevelSet) const override;
    virtual long _evalCellSupport(long id) const;

    short _evalSign(const std::array<double,3> &point) const override;
    double _evalValue(const std::array<double,3> &point, bool signedLevelSet) const override;
    std::array<double,3> _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const override;
    virtual long _evalSupport(const std::array<double,3> &point) const;

    void flushField(LevelSetField field, std::fstream &stream, VTKFormat format) const override;

private:
    void fillCartesianNarrowBandCellCaches();

    long evalCellSupport(long id, double searchRadius) const;
    long evalSupport(const std::array<double,3> &point, double searchRadius) const;

    long _evalCellSupport(long id, double searchRadius) const;
    long _evalSupport(const std::array<double,3> &point, double searchRadius) const;

};

}

#endif
