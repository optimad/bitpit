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
    class Info;
}
class SurfUnstructured;
class SendBuffer;
class RecvBuffer;

class LevelSetKernel;
class LevelSetCartesian;
class LevelSetOctree;
class LevelSetCachedObject;

class LevelSetSegmentation : public LevelSetCachedObject {

    private:
    struct SegInfo{
        std::vector<long>                segments ;                /**< list of segments within narrow band */
        std::vector<double>              distances ;               /**< list of segments distances within narrow band */

        SegInfo( ) ;
        SegInfo( const std::vector<long> &_segments, const std::vector<double> &_distances ) ;
    };

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

    typedef std::unordered_map<long, std::vector<long>> SegmentToCellMap ;

    int                                         m_dimension ;               /**< number of space dimensions */

    SurfUnstructured*                           m_segmentation;             /**< surface segmentation */
    std::unique_ptr<SurfUnstructured>           m_own;                      /**< owner of surface segmentation */
    double                                      m_featureAngle;             /**< critical angle between facets */

    std::unordered_map< long, std::vector< std::array<double,3>> > m_vertexNormal;            /**< vertex normals */
    std::unordered_map< long, std::vector< std::array<double,3>> > m_vertexGradient;            /**< vertex gradient */
    PiercedVector<SegInfo>                      m_seg;                      /**< cell -> segment association information */

    double                                      getSegmentSize( long ) const;

    protected:

    void                                        __clear() ;
    void                                        __dump( std::ostream &) ;
    void                                        __restore( std::istream &) ;

    void                                        __clearAfterMeshAdaption(const std::vector<adaption::Info> &);
    void                                        __filterOutsideNarrowBand(double);

# if BITPIT_ENABLE_MPI
    void                                        __writeCommunicationBuffer(const std::vector<long> &, SendBuffer &);
    void                                        __readCommunicationBuffer(const std::vector<long> &, RecvBuffer &);
# endif

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> &) const;
    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &);

    std::vector<std::array<double,3>>           getSimplexVertices( const long &) const;

    std::unordered_set<long>                    createSegmentInfo( LevelSetKernel *, const double &, SegmentToCellMap &) ;
    void                                        updateSegmentList( const double &) ;

    void                                        createLevelsetInfo( LevelSetKernel *, const bool &, std::unordered_set<long> &) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    SegmentToCellMap                            extractSegmentToCellMap( LevelSetCartesian *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( LevelSetOctree *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( const std::vector<adaption::Info> & ) ;

    int                                         getNarrowBandResizeDirection( LevelSetOctree *, const double &) ;

    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int,double angle=2.*M_PI);
    LevelSetSegmentation(int, std::unique_ptr<SurfUnstructured> &&, double angle=2.*M_PI);
    LevelSetSegmentation(int, SurfUnstructured*, double angle=2.*M_PI);
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    void                                        setSegmentation(std::unique_ptr<SurfUnstructured> &&) ;
    void                                        setSegmentation(SurfUnstructured *) ;
    const SurfUnstructured &                    getSegmentation() const ;
    void                                        setFeatureAngle(double) ;

    virtual int                                 getPart(const long &) const ;
    long                                        getSupport(const long &i) const;
    int                                         getSupportCount(const long &) const ;
    const std::vector<long> &                   getSimplexList(const long &) const ;

    double                                      getSurfaceFeatureSize(const long &) const;
    double                                      getMinSurfaceFeatureSize() const;
    double                                      getMaxSurfaceFeatureSize() const;

    void                                        computeLSInNarrowBand(const double &, const bool &);
    void                                        updateLSInNarrowBand(const std::vector<adaption::Info> &, const double &, const bool &) ;
};

}

#endif
