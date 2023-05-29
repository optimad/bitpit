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

#include "levelSetCommon.hpp"
#include "levelSetBoundedObject.hpp"
#include "levelSetCachedObject.hpp"
#include "levelSetCartesianKernel.hpp"
#include "levelSetKernel.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}
class SurfaceSkdTree;

class SendBuffer;
class RecvBuffer;

template<typename storage_manager_t>
class LevelSetSegmentationNarrowBandCacheBase : public virtual LevelSetNarrowBandCacheBase<storage_manager_t>
{
    public:
    using typename LevelSetNarrowBandCacheBase<storage_manager_t>::Kernel;
    using typename LevelSetNarrowBandCacheBase<storage_manager_t>::KernelIterator;

    template<typename T>
    using Storage = typename LevelSetNarrowBandCache<storage_manager_t>::template Storage<T>;

    LevelSetSegmentationNarrowBandCacheBase();

    virtual long &                                      getSupportId(const KernelIterator &itr) = 0;
    virtual long                                        getSupportId(const KernelIterator &itr) const = 0;

    virtual std::array<double, 3> &                     getSurfaceNormal(const KernelIterator &itr) = 0;
    virtual const std::array<double, 3> &               getSurfaceNormal(const KernelIterator &itr) const = 0;

    void                                        set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient) = delete ;
    void                                        set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient, long semgnetId, const std::array<double, 3> &surfaceNormal) ;

    void                                        swap(LevelSetSegmentationNarrowBandCacheBase &other) noexcept;

    protected:
    Storage<long>                              *m_supportIds;      /** Support ids of the cells inside the narrow band */
    Storage<std::array<double,3>>              *m_surfaceNormals;  /** Surface normal associated with the cells inside the narrow band */

};

template<typename storage_manager_t>
class LevelSetSegmentationNarrowBandCache : public virtual storage_manager_t, public virtual LevelSetNarrowBandCache<storage_manager_t>, public virtual LevelSetSegmentationNarrowBandCacheBase<storage_manager_t>
{

};

template<>
class LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager> : public virtual LevelSetExternalPiercedStorageManager, public virtual LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>, public virtual LevelSetSegmentationNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>
{

public:
    LevelSetSegmentationNarrowBandCache(Kernel *kernel);

    long &                                      getSupportId(const KernelIterator &itr) override;
    long                                        getSupportId(const KernelIterator &itr) const override;

    std::array<double, 3> &                     getSurfaceNormal(const KernelIterator &itr) override;
    const std::array<double, 3> &               getSurfaceNormal(const KernelIterator &itr) const override;

};

template<>
class LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager> : public virtual LevelSetInternalPiercedStorageManager, public virtual LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>, public virtual LevelSetSegmentationNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>
{

public:
    LevelSetSegmentationNarrowBandCache();

    long &                                      getSupportId(const KernelIterator &itr) override;
    long                                        getSupportId(const KernelIterator &itr) const override;

    std::array<double, 3> &                     getSurfaceNormal(const KernelIterator &itr) override;
    const std::array<double, 3> &               getSurfaceNormal(const KernelIterator &itr) const override;

};

template<>
class LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager> : public virtual LevelSetDirectStorageManager, public virtual LevelSetNarrowBandCache<LevelSetDirectStorageManager>, public virtual LevelSetSegmentationNarrowBandCacheBase<LevelSetDirectStorageManager>
{

public:
    LevelSetSegmentationNarrowBandCache(std::size_t nItems);

    long &                                      getSupportId(const KernelIterator &itr) override;
    long                                        getSupportId(const KernelIterator &itr) const override;

    std::array<double, 3> &                     getSurfaceNormal(const KernelIterator &itr) override;
    const std::array<double, 3> &               getSurfaceNormal(const KernelIterator &itr) const override;

};

template<>
class LevelSetNarrowBandCacheFactory<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>>
{

public:
    static std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>> create(LevelSetCachedObjectInterface<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>> *object);

};

template<>
class LevelSetNarrowBandCacheFactory<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>>
{

public:
    static std::shared_ptr<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>> create(LevelSetCachedObjectInterface<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>> *object);

};

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

    int getSegmentInfo( const std::array<double,3> &pointCoords, long segmentId, bool signd, double &distance, std::array<double,3> &gradient, std::array<double,3> &normal ) const;

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

class LevelSetSegmentationObjectInterface : public virtual LevelSetObjectInterface
{
public:
    virtual int getPart(long cellId) const = 0;
    virtual std::array<double BITPIT_COMMA 3> getNormal(long cellId) const = 0;
    virtual long getSupport(long id) const = 0;

    virtual double getSurfaceFeatureSize(long cellId) const = 0;
    virtual double getMinSurfaceFeatureSize() const = 0;
    virtual double getMaxSurfaceFeatureSize() const = 0;
};

template<typename narrow_band_cache_t>
class LevelSetSegmentationObject : public LevelSetSegmentationObjectInterface, public LevelSetSegmentationKernel, public LevelSetCachedObject<narrow_band_cache_t>, public LevelSetBoundedObject {

    protected:

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> &) const override;
# if BITPIT_ENABLE_MPI
    void                                        getGlobalBoundingBox( std::array<double,3> &, std::array<double,3> &) const override;
#endif

    void                                        computeNarrowBand(bool) override;
    void                                        computeNarrowBand( LevelSetCartesianKernel *, bool);
    void                                        computeNarrowBand( LevelSetKernel *, bool);
    void                                        updateNarrowBand(const std::vector<long> &cellIds, bool signd) override;
    void                                        updateNarrowBand(LevelSetKernel *, const std::vector<long> &cellIds, bool signd);

    void                                        addVTKOutputData(LevelSetField field, const std::string &objectName) override;
    std::string                                 getVTKOutputFieldName(LevelSetField field) const override;
    void                                        flushVTKOutputData(LevelSetField field, std::fstream &stream, VTKFormat format) const override;

    public:

    LevelSetSegmentationObject(int);
    LevelSetSegmentationObject(int, std::unique_ptr<const SurfUnstructured> &&, double featureAngle = 2. * BITPIT_PI);
    LevelSetSegmentationObject(int, const SurfUnstructured*, double featureAngle = 2. * BITPIT_PI);

    LevelSetSegmentationObject *                clone() const override ;

    LevelSetFieldset                            getSupportedFields() const override;

    int                                         getPart(long ) const override;
    std::array<double,3>                        getNormal(long ) const override;
    long                                        getSupport(long id) const override;

    double                                      getSurfaceFeatureSize(long ) const override;
    double                                      getMinSurfaceFeatureSize() const override;
    double                                      getMaxSurfaceFeatureSize() const override;

    LevelSetInfo                                computeLevelSetInfo(const std::array<double,3> &) const override;

};

// Typdefs for compatibility with older versions
typedef LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>> LevelSetSegmentation;

}

// Include template implementations
#include "levelSetSegmentationObject.tpp"


// Explicit instantization
#ifndef __BITPIT_LEVELSET_SEGMENTATION_OBJECT_SRC__
namespace bitpit {

extern template class LevelSetSegmentationNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>;
extern template class LevelSetSegmentationNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>;
extern template class LevelSetSegmentationNarrowBandCacheBase<LevelSetDirectStorageManager>;

extern template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
extern template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
extern template class LevelSetSegmentationObject<LevelSetSegmentationNarrowBandCache<LevelSetDirectStorageManager>>;

}
#endif

#endif
