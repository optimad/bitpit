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

#include "levelSetCartesianKernel.hpp"
#include "levelSetBooleanObject.hpp"
#include "levelSetBooleanObject.tpp"
#include "levelSetCommon.hpp"
#include "levelSetComplementObject.hpp"
#include "levelSetComplementObject.tpp"
#include "levelSetKernel.hpp"
#include "levelSetObject.hpp"

#include "bitpit_common.hpp"
#include "bitpit_surfunstructured.hpp"
#include "bitpit_volcartesian.hpp"
#include "bitpit_voloctree.hpp"

# include <array>
# include <vector>
# include <unordered_set>
# include <unordered_map>

namespace bitpit{

class SurfaceSkdTree;

class LevelSetSegmentationSurfaceInfo {

public:
    typedef SurfUnstructured::CellIterator SegmentIterator;
    typedef SurfUnstructured::CellConstIterator SegmentConstIterator;

    BITPIT_PUBLIC_API static const double DEFAULT_FEATURE_ANGLE;

    LevelSetSegmentationSurfaceInfo();
    LevelSetSegmentationSurfaceInfo(LevelSetSurfaceSmoothing surfaceSmoothing);
    LevelSetSegmentationSurfaceInfo(const LevelSetSegmentationSurfaceInfo &other);
    LevelSetSegmentationSurfaceInfo(LevelSetSegmentationSurfaceInfo &&other) = default;
    LevelSetSegmentationSurfaceInfo(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing = LevelSetSurfaceSmoothing::LOW_ORDER);
    LevelSetSegmentationSurfaceInfo(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing = LevelSetSurfaceSmoothing::LOW_ORDER);

    const SurfUnstructured & getSurface() const;
    void setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle = DEFAULT_FEATURE_ANGLE);
    void setSurface(const SurfUnstructured *surface, double featureAngle = DEFAULT_FEATURE_ANGLE);
    void setSurface(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing);

    const SurfaceSkdTree & getSearchTree() const;

    double getFeatureAngle() const;
    LevelSetSurfaceSmoothing getSurfaceSmoothing() const;

    double evalDistance(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, bool signedDistance, std::array<double, 3> *distanceVector) const;
    double evalDistance(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, bool signedDistance) const;
    std::array<double, 3> evalDistanceVector(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr) const;

    std::array<double, 3> evalNormal(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr) const;
    std::array<double,3> evalPseudoNormal(const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda) const;

    void evalProjection(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalProjection(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint) const;
    std::array<double, 3> evalProjection(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, double *lambda) const;

private:
    typedef std::pair<long, int> SegmentVertexKey;

    const SurfUnstructured *m_surface;
    std::unique_ptr<const SurfUnstructured> m_ownedSurface;
    double m_featureAngle;
    LevelSetSurfaceSmoothing m_surfaceSmoothing;

    std::unique_ptr<SurfaceSkdTree> m_searchTree;

    PiercedStorage<std::size_t> m_segmentVertexOffset;

    mutable PiercedStorage<bool> m_segmentNormalsValid;
    mutable PiercedStorage<std::array<double,3>> m_segmentNormalsStorage;
    mutable PiercedStorage<bool> m_unlimitedVertexNormalsValid;
    mutable PiercedStorage<std::array<double,3>> m_unlimitedVertexNormalsStorage;
    mutable std::vector<bool> m_limitedSegmentVertexNormalValid;
    mutable std::unordered_map<SegmentVertexKey, std::array<double,3>, utils::hashing::hash<SegmentVertexKey>> m_limitedSegmentVertexNormalStorage;

    void evalProjectionOnVertex(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalHighOrderProjectionOnLine(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalHighOrderProjectionOnTriangle(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalHighOrderProjectionOnPolygon(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalHighOrderProjection(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;

    void evalLowOrderProjectionOnLine(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalLowOrderProjectionOnTriangle(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalLowOrderProjectionOnPolygon(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;
    void evalLowOrderProjection(const std::array<double, 3> &point, const SegmentConstIterator &segmentItr, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;

    std::array<double,3> computePseudoNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const;
    std::array<double,3> computeSurfaceNormal( const SurfUnstructured::CellConstIterator &segmentIterator, const double *lambda ) const;

    std::array<double,3> computeSegmentNormal( const SurfUnstructured::CellConstIterator &segmentIterator ) const;
    std::array<double,3> computeSegmentEdgeNormal(const SegmentConstIterator &segmentItr, int edge, bool limited ) const;
    std::array<double,3> computeSegmentVertexNormal( const SurfUnstructured::CellConstIterator &segmentIterator, int vertex, bool limited ) const;

};

class LevelSetSegmentationBaseObject : public LevelSetObject {

public:
    BITPIT_PUBLIC_API static const double AUTOMATIC_SEARCH_RADIUS;

    friend class LevelSetBooleanObject<LevelSetSegmentationBaseObject>;

    using LevelSetObject::LevelSetObject;

    LevelSetFieldset getSupportedFields() const override;

    const SurfUnstructured & evalCellSurface(long id) const;
    long evalCellSupport(long id, double searchRadius = AUTOMATIC_SEARCH_RADIUS) const;
    int evalCellPart(long id) const;
    std::array<double,3> evalCellNormal(long id, bool signedLevelSet) const;

    const SurfUnstructured & evalSurface(const std::array<double,3> &point) const;
    long evalSupport(const std::array<double,3> &point) const;
    long evalSupport(const std::array<double,3> &point, double searchRadius) const;
    int evalPart(const std::array<double,3> &point) const;
    std::array<double,3> evalNormal(const std::array<double,3> &point, bool signedLevelSet) const;
    void evalProjection(const std::array<double,3> &point, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;

    BITPIT_DEPRECATED(int getPart(long cellId) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3> getNormal(long cellId) const);
    BITPIT_DEPRECATED(long getSupport(long cellId) const);
    BITPIT_DEPRECATED(double getSurfaceFeatureSize(long cellId) const);

    BITPIT_DEPRECATED(int getPart(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3> getNormal(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(long getSupport(const std::array<double,3> &point) const);
    BITPIT_DEPRECATED(double getSurfaceFeatureSize(const std::array<double,3> &point) const);

protected:
    LevelSetSegmentationBaseObject(int, const LevelSetSegmentationSurfaceInfo *surfaceInfo);

    using LevelSetObject::createFieldCellCache;
    std::size_t createFieldCellCache(LevelSetField field, std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID) override;
    void fillFieldCellCache(LevelSetField field, long id) override;
    using LevelSetObject::fillFieldCellCache;

    LevelSetIntersectionStatus _intersectSurface(long, double distance, LevelSetIntersectionMode=LevelSetIntersectionMode::FAST_FUZZY) const override;

    virtual const SurfUnstructured & _evalCellSurface(long id) const = 0;
    virtual int _evalCellPart(long id) const;
    virtual std::array<double,3> _evalCellNormal(long id, bool signedLevelSet) const = 0;
    virtual long _evalCellSupport(long id, double searchRadius = AUTOMATIC_SEARCH_RADIUS) const = 0;

    virtual const SurfUnstructured & _evalSurface(const std::array<double,3> &point) const = 0;
    virtual int _evalPart(const std::array<double,3> &point) const;
    virtual long _evalSupport(const std::array<double,3> &point) const = 0;
    virtual long _evalSupport(const std::array<double,3> &point, double searchRadius) const = 0;
    virtual void _evalProjection(const std::array<double,3> &point, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const = 0;

    void addVTKOutputData(LevelSetField field, const std::string &objectName) override;
    std::string getVTKOutputFieldName(LevelSetField field) const override;
    void flushVTKOutputData(std::fstream &stream, VTKFormat format, LevelSetField field) const override;
    using LevelSetObject::flushVTKOutputData;
};

class LevelSetSegmentationObject : public LevelSetSegmentationBaseObject {

public:
    LevelSetSegmentationObject(int);
    LevelSetSegmentationObject(int, std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle = 2. * BITPIT_PI);
    LevelSetSegmentationObject(int, const SurfUnstructured *surface, double featureAngle = 2. * BITPIT_PI);
    LevelSetSegmentationObject(int, const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing);
    LevelSetSegmentationObject(const LevelSetSegmentationObject &other);
    LevelSetSegmentationObject(LevelSetSegmentationObject &&other) = default;

    bool empty() const override;

    LevelSetSegmentationObject * clone() const override;

    const SurfUnstructured & getSurface() const;
    void setSurface(std::unique_ptr<const SurfUnstructured> &&surface, bool force = false);
    void setSurface(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle, bool force = false);
    void setSurface(const SurfUnstructured *surface, bool force = false);
    void setSurface(const SurfUnstructured *surface, double featureAngle, bool force = false);
    void setSurface(const SurfUnstructured *surface, double featureAngle, LevelSetSurfaceSmoothing surfaceSmoothing, bool force = false);

    const SurfaceSkdTree & getSearchTree() const;

    double getFeatureAngle() const;
    LevelSetSurfaceSmoothing getSurfaceSmoothing() const;

    BITPIT_DEPRECATED(double getMinSurfaceFeatureSize() const);
    BITPIT_DEPRECATED(double getMaxSurfaceFeatureSize() const);

protected:
    void fillCellLocationCache() override;
    void fillCellLocationCache(const std::vector<adaption::Info> &adaptionData) override;
    LevelSetCellLocation fillCellGeometricNarrowBandLocationCache(long id) override;

    short _evalCellSign(long id) const override;
    double _evalCellValue(long id, bool signedLevelSet) const override;
    std::array<double,3> _evalCellGradient(long id, bool signedLevelSet) const override;
    const SurfUnstructured & _evalCellSurface(long id) const override;
    long _evalCellSupport(long id, double searchRadius = AUTOMATIC_SEARCH_RADIUS) const override;
    std::array<double,3> _evalCellNormal(long id, bool signedLevelSet) const override;

    short _evalSign(const std::array<double,3> &point) const override;
    double _evalValue(const std::array<double,3> &point, bool signedLevelSet) const override;
    std::array<double,3> _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const override;
    const SurfUnstructured & _evalSurface(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point, double searchRadius) const override;
    void _evalProjection(const std::array<double,3> &point, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const override;

private:
    std::unique_ptr<LevelSetSegmentationSurfaceInfo> m_surfaceInfo;

    void fillCartesianCellZoneCache();

    short _evalSign(const std::array<double,3> &point, long support) const;
    double _evalValue(const std::array<double,3> &point, long support, bool signedLevelSet) const;
    std::array<double,3> _evalGradient(const std::array<double,3> &point, long support, bool signedLevelSet) const;
    void _evalProjection(const std::array<double,3> &point, long support, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const;

};

template<>
class LevelSetBooleanObject<LevelSetSegmentationBaseObject>: public LevelSetBooleanBaseObject<LevelSetSegmentationBaseObject> {

public:
    LevelSetBooleanObject(int, LevelSetBooleanOperation, const LevelSetSegmentationBaseObject *, const LevelSetSegmentationBaseObject *);
    LevelSetBooleanObject(int, LevelSetBooleanOperation, const std::vector<const LevelSetSegmentationBaseObject *> &);

    LevelSetBooleanObject<LevelSetSegmentationBaseObject> * clone() const override;

protected:
    const SurfUnstructured & _evalCellSurface(long id) const override;
    long _evalCellSupport(long id, double searchRadius = AUTOMATIC_SEARCH_RADIUS) const override;
    int _evalCellPart(long id) const override;
    std::array<double,3> _evalCellNormal(long id, bool signedLevelSet) const override;

    const SurfUnstructured & _evalSurface(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point, double searchRadius) const override;
    int _evalPart(const std::array<double,3> &point) const override;
    void _evalProjection(const std::array<double,3> &point, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const override;

};

template<>
class LevelSetComplementObject<LevelSetSegmentationBaseObject>: public LevelSetComplementBaseObject<LevelSetSegmentationBaseObject> {

public:
    LevelSetComplementObject(int id, const LevelSetSegmentationBaseObject *source);

    LevelSetComplementObject<LevelSetSegmentationBaseObject> * clone() const override;

protected:
    const SurfUnstructured & _evalCellSurface(long id) const override;
    long _evalCellSupport(long id, double searchRadius = AUTOMATIC_SEARCH_RADIUS) const override;
    int _evalCellPart(long id) const override;
    std::array<double,3> _evalCellNormal(long id, bool signedLevelSet) const override;

    const SurfUnstructured & _evalSurface(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point) const override;
    long _evalSupport(const std::array<double,3> &point, double searchRadius) const override;
    int _evalPart(const std::array<double,3> &point) const override;
    void _evalProjection(const std::array<double,3> &point, bool signedLevelSet, std::array<double, 3> *projectionPoint, std::array<double, 3> *projectionNormal) const override;

};

// Typedefs for compatibility with older versions
typedef LevelSetSegmentationObject LevelSetSegmentation;

}

#endif
