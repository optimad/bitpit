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

# ifndef __BITPIT_LEVELSET_SEGMENTATION_HPP__
# define __BITPIT_LEVELSET_SEGMENTATION_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_set>
# include <unordered_map>

#include "levelSetCommon.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}
class SurfUnstructured;
class SurfaceSkdTree;

class SendBuffer;
class RecvBuffer;

class LevelSetKernel;
class LevelSetCartesian;
class LevelSetOctree;
class LevelSetCachedObject;

class SegmentationKernel {

public:
    SegmentationKernel();
    SegmentationKernel(const SurfUnstructured *surface, double featureAngle);
    SegmentationKernel(std::unique_ptr<const SurfUnstructured> &&surface, double featureAngle);

    const SurfUnstructured & getSurface() const;
    double getFeatureAngle() const;

    const SurfaceSkdTree & getSearchTree() const ;

    int getSegmentInfo( const std::array<double,3> &p, long i, bool signd, double &d, std::array<double,3> &x, std::array<double,3> &n ) const;

private:
    const SurfUnstructured *m_surface;
    std::unique_ptr<const SurfUnstructured> m_ownedSurface;
    double m_featureAngle;

    std::unique_ptr<SurfaceSkdTree> m_searchTree;

    PiercedStorage<std::vector< std::array<double,3>>> m_segmentVertexNormals;

    void setSurface( const SurfUnstructured *surface, double featureAngle);

};

class LevelSetSegmentation : public LevelSetCachedObject {

    private:
    struct DistanceComparator
    {
        const std::vector<double> & m_vector;

        DistanceComparator(const std::vector<double> & vector)
            : m_vector(vector)
        {

        }

        bool operator()(double i1, double i2)
        {
            return m_vector[i1] < m_vector[i2];
        }
    };

    struct SurfaceInfo
    {
        long support;
        std::array<double,3> normal;

        SurfaceInfo();
        SurfaceInfo(long, const std::array<double,3> &);
    };

    std::shared_ptr<const SegmentationKernel> m_segmentation;
    PiercedVector<SurfaceInfo>                         m_surfaceInfo;                      /**< cell support information  */

    double                                      getSegmentSize( long ) const;

    protected:

    void                                        __clear() override ;
    void                                        __dump( std::ostream &) override ;
    void                                        __restore( std::istream &) override ;

    void                                        __clearAfterMeshAdaption(const std::vector<adaption::Info> &) override ;

# if BITPIT_ENABLE_MPI
    void                                        __writeCommunicationBuffer(const std::vector<long> &, SendBuffer &) override ;
    void                                        __readCommunicationBuffer(const std::vector<long> &, RecvBuffer &)  override ;
# endif

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> &) const override;
    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, double, std::vector<long> &);

    void                                        computeLSInNarrowBand( LevelSetCartesian *, bool);
    void                                        computeLSInNarrowBand( LevelSetOctree *, bool);
    void                                        updateLSInNarrowBand(LevelSetOctree *, const std::vector<adaption::Info> &, bool);

    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int);
    LevelSetSegmentation(int, std::unique_ptr<const SurfUnstructured> &&, double featureAngle = 2. * BITPIT_PI);
    LevelSetSegmentation(int, const SurfUnstructured*, double featureAngle = 2. * BITPIT_PI);

    LevelSetSegmentation*                       clone() const override ;

    void                                        setSegmentation(std::unique_ptr<const SurfUnstructured> &&patch, double featureAngle = 2. * BITPIT_PI) ;
    void                                        setSegmentation(const SurfUnstructured *patch, double featureAngle = 2. * BITPIT_PI) ;
    const SegmentationKernel &                  getSegmentation() const ;

    virtual int                                 getPart(long ) const override;
    virtual std::array<double,3>                getNormal(long ) const override;
    long                                        getSupport(long i) const;

    double                                      getSurfaceFeatureSize(long ) const override;
    double                                      getMinSurfaceFeatureSize() const override;
    double                                      getMaxSurfaceFeatureSize() const override;

    void                                        computeLSInNarrowBand(bool) override;
    void                                        updateLSInNarrowBand(const std::vector<adaption::Info> &, bool) override;

    LevelSetInfo                                computeLevelSetInfo(const std::array<double,3> &) const override;

};

}

#endif
