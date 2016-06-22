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
    const std::unordered_set<long>          LIST;                       /**< Default value for list of elements in narrow band */
};

class LevelSetInfo{
    public:
    double                                  value ;                     /**< Levelset value */
    std::array<double,3>                    gradient ;                  /**< Levelset gradient */
    int                                     object ;                    /**< Id of closest object */
    int                                     part  ;                     /**< Id of closest patch */
    long                                    support ;                   /**< Id of closest support */

    LevelSetInfo() ;
};

class LevelSet{

    private:
    LevelSetKernel*                             m_kernel ;              /**< LevelSet computational kernel */
    std::unordered_map<int,LevelSetObject*>     m_object ;              /**< Objects defining the boundaries */

    bool                                        m_userRSearch;          /**< Flag if user has set size of narrow band (default=false)  */
    bool                                        m_signedDF;             /**< Flag for sigend/unsigned distance function (default = true) */
    bool                                        m_propagateS;           /**< Flag for sign propagation from narrow band (default = false) */
    bool                                        m_propagateV;           /**< Flag for value propagation from narrow band (default = false) */

    public:
    ~LevelSet() ;
    LevelSet() ;

    void                                        setMesh( VolumeKernel* ) ;
    void                                        setMesh( VolCartesian* ) ;
    void                                        setMesh( VolOctree* ) ;

    int                                         addObject( std::unique_ptr<SurfaceKernel> &&, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( SurfaceKernel *, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( std::unique_ptr<SurfUnstructured> &&, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( SurfUnstructured *, int id = levelSetDefaults::OBJECT ) ;
    int                                         addObject( LevelSetObject* ) ;
    const LevelSetObject &                      getObject( int ) const ;
    int                                         getObjectCount( ) const ;

    void                                        clear();
    void                                        clearObject();
    void                                        clearObject( int );

    LevelSetInfo                                getLevelSetInfo(const long &) const ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;
    int                                         getClosestObject(const long &) const ;
    std::pair<int,int>                          getClosestPart(const long &) const ;
    std::pair<int,long>                         getClosestSupport(const long &) const ;
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
    std::pair<int,long>                         getClosestSupport(const long &) const ;
    short                                       getSign(const long &) const;
    double                                      getSizeNarrowBand() const;
    bool                                        isInNarrowBand(const long &) const;

    void                                        setSizeNarrowBand(double) ;

    virtual double                              computeSizeNarrowBand( LevelSetObject * )=0;
    virtual double                              computeSizeNarrowBandFromLS(  );
    virtual double                              updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, LevelSetObject *> &objects )=0;
    virtual double                              computeRSearchFromCell( long id ) = 0;

    void                                        clear() ;
    void                                        clearAfterMeshMovement( const std::vector<adaption::Info> & ) ;
    void                                        filterOutsideNarrowBand( double ) ;

    void                                        propagateSign( std::unordered_map<int,LevelSetObject*> ) ;
    void                                        propagateValue( LevelSetObject *) ;

    void                                        dump( std::fstream &);
    void                                        restore( std::fstream &);

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
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, LevelSetObject *> &objects );
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
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> &, std::unordered_map<int, LevelSetObject *> &objects );
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

    virtual void                                computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &)=0 ;
    virtual void                                updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool &)=0 ;
    virtual void                                clearAfterMeshMovement( const std::vector<adaption::Info> & ) ;
    virtual void                                filterOutsideNarrowBand( LevelSetKernel *) ;

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
        std::unordered_set<long>                m_segments ;                /**< list of segments within narrow band */
        bool                                    m_checked ;                 /**< true if list of segments has been processed */

        SegInfo( ) ;
        SegInfo( const std::unordered_set<long> & ) ;
    };

    int                                         m_dimension ;               /**< number of space dimensions */

    SurfUnstructured*                           m_segmentation;             /**< surface segmentation */
    std::unique_ptr<SurfUnstructured>           m_own;                      /**< owner of surface segmentation */

    std::unordered_map< long, std::vector< std::array<double,3>> > m_vertexNormal;            /**< vertex normals */
    PiercedVector<SegInfo>                      m_seg;                      /**< cell -> segment association information */


    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int);
    LevelSetSegmentation(int, std::unique_ptr<SurfUnstructured> && );
    LevelSetSegmentation(int, SurfUnstructured* );
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    void                                        setSegmentation( std::unique_ptr<SurfUnstructured> && ) ;
    void                                        setSegmentation( SurfUnstructured * ) ;
    const SurfUnstructured &                    getSegmentation() const ;

    const std::unordered_set<long> &            getSimplexList(const long &) const ;
    bool                                        isInNarrowBand( const long &) ;

    void                                        dumpDerived( std::fstream &) ;
    void                                        restoreDerived( std::fstream &) ;

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;

    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;
    double                                      evaluateLS( LevelSetKernel *, long) const ;

    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &);
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool & ) ;
    void                                        clearAfterMeshMovement( const std::vector<adaption::Info> & ) ;
    void                                        filterOutsideNarrowBand( LevelSetKernel *) ;

# if BITPIT_ENABLE_MPI
    void                                        writeCommunicationBuffer( const std::vector<long> &, SendBuffer &, SendBuffer & ) ;
    void                                        readCommunicationBuffer(  const std::vector<long> &, const long &, RecvBuffer & ) ;
# endif

    protected:
    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;
    void                                        lsFromSimplex( LevelSetKernel *, const double &, const bool &, bool filter = false ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;
    void                                        associateSimplexToCell( LevelSetCartesian *, const double &);
    void                                        associateSimplexToCell( LevelSetOctree *, const double &);

    void                                        updateSimplexToCell( LevelSetOctree *, const std::vector<adaption::Info> &, const double & ) ;
};

}

#endif /* LEVELSET_HPP */
