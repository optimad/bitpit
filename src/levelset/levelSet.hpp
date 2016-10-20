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

enum class LevelSetBooleanOperation{
    UNION =0,
    INTERSECTION =1,
    SUBTRACTION =2
};

class LevelSetInfo{
    public:
    double                                  value ;                     /**< Levelset value */
    std::array<double,3>                    gradient ;                  /**< Levelset gradient */

    LevelSetInfo() ;
    LevelSetInfo( const double &, const std::array<double,3> &) ;
};

class LevelSet{

    private:
    std::unique_ptr<LevelSetKernel>                             m_kernel ;              /**< LevelSet computational kernel */
    std::unordered_map<int,std::unique_ptr<LevelSetObject>>     m_object ;              /**< Objects defining the boundaries */

    std::vector<int>                            m_order ;               /**< Processing order of objects */
    bool                                        m_userRSearch;          /**< Flag if user has set size of narrow band (default=false)  */
    bool                                        m_signedDF;             /**< Flag for sigend/unsigned distance function (default = true) */
    bool                                        m_propagateS;           /**< Flag for sign propagation from narrow band (default = false) */

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
    int                                         addObject( const LevelSetBooleanOperation &, const int &, const int &, int id=levelSetDefaults::OBJECT ) ;
    void                                        addProcessingOrder( int) ;
    const LevelSetObject &                      getObject( int ) const ;
    int                                         getObjectCount( ) const ;
    std::vector<int>                            getObjectIds( ) const ;

    void                                        clear();
    void                                        clearObject();
    void                                        clearObject( int );

    int                                         getClosestObject(const long &) const ;
    LevelSetInfo                                getLevelSetInfo(const long &, int objectId=levelSetDefaults::OBJECT) const ;
    double                                      getLS(const long &, int objectId=levelSetDefaults::OBJECT) const;
    std::array<double,3>                        getGradient(const long &, int objectId=levelSetDefaults::OBJECT) const ;
    std::pair<int,int>                          getClosestPart(const long &, int objectId=levelSetDefaults::OBJECT) const ;
    std::pair<int,long>                         getClosestSupport(const long &, int objectId=levelSetDefaults::OBJECT) const ;
    int                                         getSupportCount(const long &, int objectId=levelSetDefaults::OBJECT) const ;
    short                                       getSign(const long &, int objectId=levelSetDefaults::OBJECT) const;
    bool                                        isInNarrowBand(const long &, int ) const ;

    double                                      getSizeNarrowBand( const int &) const;
    void                                        setSizeNarrowBand(double) ;

    void                                        setSign(bool);
    void                                        setPropagateSign(bool) ;

    void                                        dump( std::ostream &);
    void                                        restore( std::istream &);

    void                                        compute( ) ;
    void                                        update( const std::vector<adaption::Info> & ) ;
    private:
# if BITPIT_ENABLE_MPI
    bool                                        assureMPI() ;
# endif

};

class LevelSetKernel{

    private:
    std::unordered_map<long, std::array<double,3>> m_cellCentroids; /**< Cached cell center coordinates*/

    protected:
    VolumeKernel*                               m_mesh ;        /**< Pointer to underlying mesh*/

# if BITPIT_ENABLE_MPI
    MPI_Comm                                    m_commMPI ;     /**< MPI communicator */
# endif

    public:
    virtual ~LevelSetKernel() ;
    LevelSetKernel() ;
    LevelSetKernel( VolumeKernel *) ;

    VolumeKernel*                               getMesh() const ;

    virtual double                              computeSizeNarrowBandFromLS( LevelSetObject*, const bool & );
    virtual double                              computeRSearchFromCell( long id ) = 0;

    void                                        clearGeometryCache(  ) ;
    void                                        updateGeometryCache( const std::vector<adaption::Info> &mapper ) ;

    const std::array<double,3> &                computeCellCentroid( long id ) ;
    double                                      isCellInsideBoundingBox( long id, std::array<double, 3> minPoint, std::array<double, 3> maxPoint );

# if BITPIT_ENABLE_MPI
    MPI_Comm                                    getCommunicator() const ;
    void                                        freeCommunicator();
    bool                                        isCommunicatorSet() const;
    bool                                        assureMPI() ;
# endif

};

class LevelSetCartesian : public LevelSetKernel{

    private:
    VolCartesian*                               m_cartesian ;       /**< Pointer to underlying cartesian mesh*/

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( VolCartesian & );

    VolCartesian *                              getCartesianMesh() const;
    double                                      computeRSearchFromCell( long id ) ;

};

class LevelSetOctree : public LevelSetKernel{

    private:
    VolOctree*                                  m_octree ;       /**< Pointer to underlying octree mesh*/

    public:
    virtual ~LevelSetOctree();
    LevelSetOctree( VolOctree & );

    VolOctree *                                 getOctreeMesh() const;
    double                                      computeRSearchFromCell( long id ) ;

    double                                      computeRSearchFromLevel( uint8_t ) ;
    double                                      computeSizeFromRSearch( double ) ;
};

class LevelSetObject{

    private:
    int                                         m_id;           /**< identifier of object */
    bool                                        m_primary;      /**< identifier of object */


    protected:
    double                                      m_RSearch;      /**< Size of narrow band */

    virtual void                                _clear();
    virtual void                                _clearAfterMeshAdaption(const std::vector<adaption::Info>&) ;
    virtual void                                _filterOutsideNarrowBand(double) ;
    virtual void                                _dump(std::ostream &);
    virtual void                                _restore( std::istream &);

# if BITPIT_ENABLE_MPI
    virtual void                                writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ; 
    virtual void                                _writeCommunicationBuffer(const std::vector<long>&, SendBuffer&)  ;
    virtual void                                readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ; 
    virtual void                                _readCommunicationBuffer(const std::vector<long>&, RecvBuffer&)  ;
# endif 

    public:
    virtual ~LevelSetObject();
    LevelSetObject(int,bool);

    virtual LevelSetObject*                     clone() const =0;
    void                                        clear();

    int                                         getId() const ;
    bool                                        isPrimary() const ;

    virtual LevelSetInfo                        getLevelSetInfo(const long &) const =0;
    virtual double                              getLS(const long &) const =0; 
    virtual std::array<double,3>                getGradient(const long &) const =0 ; 

    virtual int                                 getPart(const long &) const ;
    virtual long                                getSupport(const long &) const;
    virtual int                                 getSupportCount(const long &) const ;

    short                                       getSign(const long &) const;
    virtual void                                propagateSign(LevelSetKernel*) ; 

    bool                                        isInNarrowBand(const long &) const;
    double                                      getSizeNarrowBand() const;
    virtual void                                setSizeNarrowBand(double) ;

    virtual double                              computeSizeNarrowBand(LevelSetKernel*)=0;
    virtual void                                computeLSInNarrowBand(LevelSetKernel*, const double &, const bool &) ;

    virtual double                              updateSizeNarrowBand(LevelSetKernel*, const std::vector<adaption::Info> &);
    virtual void                                updateLSInNarrowBand(LevelSetKernel*, const std::vector<adaption::Info> &, const double &, const bool &) ;
    void                                        clearAfterMeshAdaption(const std::vector<adaption::Info>&);
    void                                        filterOutsideNarrowBand(double);  ;

    void                                        dump(std::ostream &); 
    void                                        restore(std::istream &); 

# if BITPIT_ENABLE_MPI
    bool                                        assureMPI(LevelSetKernel * ) ;
    void                                        exchangeGhosts(LevelSetKernel * ) ;
    void                                        communicate( LevelSetKernel *, std::unordered_map<int,std::vector<long>> &, std::unordered_map<int,std::vector<long>> &, std::vector<adaption::Info> const *mapper=NULL ) ;
# endif 

};

class LevelSetCachedObject : public LevelSetObject{

    private:
    void                                        assignSign( int sign, const std::unordered_set<long> &cells ) ;

    protected:
    PiercedVector<LevelSetInfo>                 m_ls ;          /**< Levelset information for each cell */
    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;

    void                                        _clear( ) ;
    virtual void                                __clear() ;

    double                                      _computeSizeNarrowBand(LevelSetCartesian*);
    double                                      _computeSizeNarrowBand(LevelSetOctree*);

    double                                      _updateSizeNarrowBand(LevelSetCartesian*, const std::vector<adaption::Info> &);
    double                                      _updateSizeNarrowBand(LevelSetOctree*, const std::vector<adaption::Info> &);
    void                                        _clearAfterMeshAdaption(const std::vector<adaption::Info> & ) ;
    virtual void                                __clearAfterMeshAdaption(const std::vector<adaption::Info> & ) ;
    void                                        _filterOutsideNarrowBand(double);
    virtual void                                __filterOutsideNarrowBand(double);

    void                                        _dump( std::ostream &) ; 
    virtual void                                __dump(std::ostream &) ; 
    void                                        _restore( std::istream &) ; 
    virtual void                                __restore(std::istream &) ; 

# if BITPIT_ENABLE_MPI
    void                                        _writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ; 
    virtual void                                __writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ; 
    void                                        _readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ; 
    virtual void                                __readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ; 
# endif 

    public:
    virtual ~LevelSetCachedObject();
    LevelSetCachedObject(int);

    LevelSetInfo                                getLevelSetInfo(const long &) const ;
    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;

    void                                        propagateSign(LevelSetKernel*);

    double                                      computeSizeNarrowBand(LevelSetKernel*);
    double                                      updateSizeNarrowBand(LevelSetKernel*, const std::vector<adaption::Info> &);
};

class LevelSetBoolean: public LevelSetObject{

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    LevelSetObject*                             m_objPtr1;              /**< pointer to first object */
    LevelSetObject*                             m_objPtr2;              /**< pointer to second object */
    int                                         m_objId1;               /**< identifier of first object */
    int                                         m_objId2;               /**< identifier of second object */

    protected:
    void                                        _dump( std::ostream &);
    void                                        _restore( std::istream &);

    public:
    ~LevelSetBoolean();
    LevelSetBoolean(int, LevelSetBooleanOperation, LevelSetObject*, LevelSetObject*);
    LevelSetBoolean(const LevelSetBoolean &);

    LevelSetBoolean*                            clone() const ;

    LevelSetInfo                                getLevelSetInfo(const long &) const;
    double                                      getLS(const long &) const; 
    std::array<double,3>                        getGradient(const long &) const; 

    int                                         getPart(const long &) const ;
    long                                        getSupport(const long &) const;
    int                                         getSupportCount(const long &) const ;

    void                                        setSizeNarrowBand(double) ;

    double                                      computeSizeNarrowBand(LevelSetKernel*);
    double                                      updateSizeNarrowBand(LevelSetKernel*,const std::vector<adaption::Info> &);
    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &) ;
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool &);

    LevelSetBooleanOperation                    getBooleanOperation() const;
    LevelSetObject*                             getClosestObject(const long &) const ;
    LevelSetInfo                                booleanOperation(const long &) const ;
};

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

    protected:

    void                                        __clear() ;
    void                                        __dump( std::ostream &) ;
    void                                        __restore( std::istream &) ;

    void                                        __clearAfterMeshAdaption( const std::vector<adaption::Info> & ) ;
    void                                        __filterOutsideNarrowBand( double) ;

# if BITPIT_ENABLE_MPI
    void                                        __writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ;
    void                                        __readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ;
# endif

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;
    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;

    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;

    std::unordered_set<long>                    createSegmentInfo( LevelSetKernel *visitee, const double &search, SegmentToCellMap &segmentToCellMap ) ;
    void                                        updateSegmentList( const double &search ) ;

    void                                        createLevelsetInfo( LevelSetKernel *visitee, const bool & signd, std::unordered_set<long> &cellList ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    SegmentToCellMap                            extractSegmentToCellMap( LevelSetCartesian *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( LevelSetOctree *, const double &);
    SegmentToCellMap                            extractSegmentToCellMap( const std::vector<adaption::Info> & ) ;

    int                                         getNarrowBandResizeDirection( LevelSetOctree *visitee, const double &newRSearch ) ;

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

    virtual int                                 getPart(const long &) const ;
    long                                        getSupport(const long &i) const;
    int                                         getSupportCount(const long &) const ;
    const std::vector<long> &                   getSimplexList(const long &) const ;

    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &);
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool & ) ;
};

}

#endif /* LEVELSET_HPP */
