/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

# ifndef __BITPIT_LEVELSET_HPP__
# define __BITPIT_LEVELSET_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_set>
# include <unordered_map>

# if BITPIT_ENABLE_MPI
# include "bitpit_communications.hpp"
# endif
# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"

namespace bitpit{

class LevelSetKernel ;
class LevelSetObject ;
class LevelSetSegmentation ;

namespace levelSetDefaults{
    const double                            VALUE = 1.e18 ;             /**< Default value for levelset function */
    const std::array<double,3>              GRADIENT = {{0.,0.,0.}};    /**< Default value for levelset gradient */
    const short                             SIGN = 1;                   /**< Default value for the sign */
    const int                               OBJECT = -1 ;               /**< Default value for closest object  */
    const int                               PART  = -1 ;               /**< Default value for closest patch  */
    const long                              SUPPORT = -1 ;              /**< Default value for closest support element */
    const std::vector<long>                 LIST;                       /**< Default value for list of elements in narrow band */
};

class LevelSetInfo{
    public:
    double                                  value ;                     /**< Levelset value */
    std::array<double,3>                    gradient ;                  /**< Levelset gradient */
    int                                     object ;                    /**< Id of closest object */
    int                                     part  ;                     /**< Id of closest patch */

    LevelSetInfo() ;
};

class LevelSet{

    private:
    std::unique_ptr<LevelSetKernel>                             m_kernel ;              /**< LevelSet computational kernel */
    std::unordered_map<int,std::unique_ptr<LevelSetObject>>     m_object ;              /**< Objects defining the boundaries */

    bool                                        m_userRSearch;          /**< Flag if user has set size of narrow band (default=false)  */
    bool                                        m_signedDF;             /**< Flag for sigend/unsigned distance function (default = true) */
    bool                                        m_propagateS;           /**< Flag for sign propagation from narrow band (default = false) */
    bool                                        m_propagateV;           /**< Flag for value propagation from narrow band (default = false) */

    public:
    ~LevelSet() ;
    LevelSet() ;

    LevelSet(LevelSet&& other) = default;

    void                                        setMesh( VolumeKernel* ) ;
    void                                        setMesh( VolCartesian* ) ;
    void                                        setMesh( VolOctree* ) ;

    int                                         addObject( std::unique_ptr<SurfaceKernel> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( SurfaceKernel *, double, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( std::unique_ptr<SurfUnstructured> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( SurfUnstructured *, double, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( std::unique_ptr<LevelSetObject> && ) ;
    int                                         addObject( const std::unique_ptr<LevelSetObject> & ) ;
    const LevelSetObject &                      getObject( int ) const ;
    int                                         getObjectCount( ) const ;
    std::vector<int>                            getObjectIds( ) const ;

    void                                        clear();
    void                                        clearObject();
    void                                        clearObject( int );

    LevelSetInfo                                getLevelSetInfo(const long &) const ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;
    int                                         getClosestObject(const long &) const ;
    std::pair<int,int>                          getClosestPart(const long &) const ;
    std::pair<int,long>                         getClosestSupport(const long &) const ;
    int                                         getSupportCount(const long &) const ;
    short                                       getSign(const long &) const;
    bool                                        isInNarrowBand(const long &) const;

    double                                      getSizeNarrowBand() const;

    void                                        setSizeNarrowBand(double) ;
    void                                        setSign(bool);
    void                                        setPropagateSign(bool) ;
    void                                        setPropagateValue(bool) ;

    void                                        dump( std::fstream &);
    void                                        restore( std::fstream &);

    void                                        compute( ) ;
    void                                        update( const std::vector<adaption::Info> & ) ;
# if BITPIT_ENABLE_MPI
    void                                        exchangeGhosts( ) ;
    void                                        communicate( std::unordered_map<int,std::vector<long>> &, std::unordered_map<int,std::vector<long>> &, std::vector<adaption::Info> const *mapper=NULL ) ;
# endif

    private:
# if BITPIT_ENABLE_MPI
    bool                                        assureMPI() ;
# endif

};

class LevelSetKernel{

    public:

    protected:
    PiercedVector<LevelSetInfo>                 m_ls ;          /**< Levelset information for each cell */
    VolumeKernel*                               m_mesh ;        /**< Pointer to underlying mesh*/

    double                                      m_RSearch;      /**< Size of narrow band */

# if BITPIT_ENABLE_MPI
    MPI_Comm                                    m_commMPI ;     /**< MPI communicator */
# endif

    private:
    std::unordered_map<long, std::array<double,3>> m_cellCentroids;

    public:
    virtual ~LevelSetKernel() ;
    LevelSetKernel() ;
    LevelSetKernel( VolumeKernel *) ;

    VolumeKernel*                               getMesh() const ;

    PiercedVector<LevelSetInfo>&                getLevelSetInfo() ;
    LevelSetInfo                                getLevelSetInfo(const long &) const ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;
    int                                         getClosestObject(const long &) const ;
    std::pair<int,int>                          getClosestPart(const long &) const ;
    int                                         getSupportCount(const long &) const ;
    short                                       getSign(const long &) const;
    double                                      getSizeNarrowBand() const;
    bool                                        isInNarrowBand(const long &) const;

    void                                        setSizeNarrowBand(double) ;

    virtual double                              computeSizeNarrowBand( LevelSetObject * )=0;
    virtual double                              computeSizeNarrowBandFromLS(  );
    virtual double                              updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, std::unique_ptr<LevelSetObject>> &objects )=0;
    virtual double                              computeRSearchFromCell( long id ) = 0;

    void                                        clear() ;
    void                                        clearAfterMeshAdaption( const std::vector<adaption::Info> & ) ;
    void                                        filterOutsideNarrowBand( double ) ;

    void                                        propagateSign( std::unordered_map<int, std::unique_ptr<LevelSetObject>> & ) ;
    void                                        propagateValue( LevelSetObject *) ;

    void                                        dump( std::fstream &);
    void                                        restore( std::fstream &);

    void                                        clearGeometryCache(  ) ;
    void                                        updateGeometryCache( const std::vector<adaption::Info> &mapper ) ;

    const std::array<double,3> &                computeCellCentroid( long id ) ;
    double                                      isCellInsideBoundingBox( long id, std::array<double, 3> minPoint, std::array<double, 3> maxPoint );

# if BITPIT_ENABLE_MPI

    MPI_Comm                                    getCommunicator() const ;
    void                                        freeCommunicator();
    bool                                        isCommunicatorSet() const;
    bool                                        assureMPI() ;
    void                                        writeCommunicationBuffer( const std::vector<long> &, SendBuffer &, SendBuffer & );
    void                                        readCommunicationBuffer(  const std::vector<long> &, const long &, RecvBuffer & ) ;
# endif

    protected:
    void                                        solveEikonal( double, double );
    virtual double                              updateEikonal( double, double, const long &, const std::unordered_map<long,short> & ) ; 

    std::array<double,3>                        computeGradientUpwind(const long &) ;
    std::array<double,3>                        computeGradientCentral(const long &) ;



};

class LevelSetCartesian : public LevelSetKernel{

    private:
    VolCartesian*                               m_cartesian ;       /**< Pointer to underlying cartesian mesh*/

    private:
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, std::unique_ptr<LevelSetObject>> &objects );
    double                                      computeRSearchFromCell( long id ) ;

    double                                      updateEikonal( double, double, const long &, const std::unordered_map<long,short> & ) ; 

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( VolCartesian & );

    double                                      computeSizeNarrowBand( LevelSetObject * );
};

class LevelSetOctree : public LevelSetKernel{

    private:
    VolOctree*                                  m_octree ;       /**< Pointer to underlying octree mesh*/

    public:
    double                                      computeRSearchFromLevel( uint8_t ) ;
    double                                      computeSizeFromRSearch( double ) ;

    virtual ~LevelSetOctree();
    LevelSetOctree( VolOctree & );

    double                                      computeSizeNarrowBand( LevelSetObject * );
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, std::unique_ptr<LevelSetObject>> &objects );
    double                                      computeRSearchFromCell( long id ) ;

};

class LevelSetObject{


    private:
    int                                         m_id;           /**< identifier of object */

    public:
    virtual ~LevelSetObject();
    LevelSetObject(int);

    virtual int                                 getId() const ;
    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;
    virtual LevelSetObject*                     clone() const = 0;

    virtual void                                clear( LevelSetKernel *visitee = nullptr )=0 ;

    virtual void                                computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &)=0 ;
    virtual void                                updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool &)=0 ;
    virtual void                                clearAfterMeshAdaption( const std::vector<adaption::Info> & ) ;
    virtual void                                filterOutsideNarrowBand( LevelSetKernel *) ;
    virtual int                                 getSupportCount(const long &) const ;
    virtual long                                getClosestSupport(const long &i) const;

    virtual void                                dumpDerived( std::fstream &) =0 ;
    virtual void                                restoreDerived( std::fstream &) =0 ;

    virtual double                              evaluateLS( LevelSetKernel *, long) const =0;

    void                                        dump( std::fstream &) ;
    void                                        restore( std::fstream &) ;

# if BITPIT_ENABLE_MPI
    virtual void                                writeCommunicationBuffer( const std::vector<long> &, SendBuffer &, SendBuffer & ) =0 ;
    virtual void                                readCommunicationBuffer(  const std::vector<long> &, const long &, RecvBuffer & ) =0 ;
# endif 

};

class LevelSetSegmentation : public LevelSetObject {

    private:
    struct SegInfo{
        std::vector<long>                segments ;                /**< list of segments within narrow band */
        std::vector<double>              distances ;               /**< list of segments distances within narrow band */

        SegInfo( ) ;
        SegInfo( const std::vector<long> &_segments, const std::vector<double> &_distances ) ;
    };

    struct DistanceComparator
    {
        const vector<double> & m_vector;

        DistanceComparator(const vector<double> & vector)
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


    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int,double angle=2.*M_PI);
    LevelSetSegmentation(int, std::unique_ptr<SurfUnstructured> &&, double angle=2.*M_PI);
    LevelSetSegmentation(int, SurfUnstructured*, double angle=2.*M_PI);
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    void                                        setSegmentation( std::unique_ptr<SurfUnstructured> && ) ;
    void                                        setSegmentation( SurfUnstructured * ) ;
    const SurfUnstructured &                    getSegmentation() const ;
    void                                        setFeatureAngle(double) ;

    const std::vector<long> &                   getSimplexList(const long &) const ;
    bool                                        isInNarrowBand( const long &) ;

    void                                        dumpDerived( std::fstream &) ;
    void                                        restoreDerived( std::fstream &) ;

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;

    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;
    double                                      evaluateLS( LevelSetKernel *, long) const ;

    void                                        clear( LevelSetKernel *visitee = nullptr ) ;

    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &);
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool & ) ;
    void                                        clearAfterMeshAdaption( const std::vector<adaption::Info> & ) ;
    void                                        filterOutsideNarrowBand( LevelSetKernel *) ;
    int                                         getSupportCount(const long &) const ;
    long                                        getClosestSupport(const long &i) const;

# if BITPIT_ENABLE_MPI
    void                                        writeCommunicationBuffer( const std::vector<long> &, SendBuffer &, SendBuffer & ) ;
    void                                        readCommunicationBuffer(  const std::vector<long> &, const long &, RecvBuffer & ) ;
# endif

    protected:
    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;

    std::unordered_set<long>                    createSegmentInfo( LevelSetKernel *visitee, const double &search, SegmentToCellMap &segmentToCellMap ) ;
    void                                        updateSegmentList( const double &search ) ;

    void                                        createLevelsetInfo( LevelSetKernel *visitee, const bool & signd, std::unordered_set<long> &cellList ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    SegmentToCellMap                            extractSegmentToCellMap( LevelSetCartesian *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( LevelSetOctree *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( const std::vector<adaption::Info> & ) ;

    int                                         getNarrowBandResizeDirection( LevelSetOctree *visitee, const double &newRSearch ) ;
};

}

#endif /* LEVELSET_HPP */
