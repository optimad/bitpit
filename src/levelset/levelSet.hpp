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
# include <cmath>
# include <array>
# include <vector>
# include <deque>
# include <string>
# include <iostream>
# include <set>

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
    const int                               OBJECT = -1 ;               /**< Default value for closest object id */
    const std::unordered_set<long>          LIST = { } ;                /**< Default value for closest segment */
    const long                              ELEMENT = -1 ;              /**< Default value for segmments in narrow band */
};

class LevelSet{

    private:
    LevelSetKernel*                             m_kernel ;              /**< LevelSet computational kernel */
    std::vector<LevelSetObject*>                m_object ;              /**< Objects defining the boundaries */

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

    void                                        addObject( SurfaceKernel* ) ;
    void                                        addObject( SurfUnstructured * ) ;
    void                                        addObject( LevelSetObject* ) ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;
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
    void                                        loadBalance( const std::vector<adaption::Info> & ) ;
    void                                        finalizeMPI() ;
# endif

    private:
# if BITPIT_ENABLE_MPI
    bool                                        assureMPI() ;
# endif

};

class LevelSetKernel{

    friend LevelSet;

    protected:
    struct LSInfo{
        double                                  value ;         /**< Levelset value */
        std::array<double,3>                    gradient ;      /**< Levelset gradient */
        int                                     object ;        /**< Id of closest object */
        short                                   active ;        /**< Flag for nodes for fast marching */
    };

    PiercedVector<LSInfo>                       m_ls ;          /**< Levelset information for each cell */
    VolumeKernel*                               m_mesh ;        /**< Pointer to underlying mesh*/

    double                                      m_RSearch;      /**< Size of narrow band */

# if BITPIT_ENABLE_MPI
    MPI_Comm                                    m_commMPI ;     /**< MPI communicator */
# endif

    public:
    virtual ~LevelSetKernel() ;
    LevelSetKernel() ;
    LevelSetKernel( VolumeKernel *) ;

    PiercedVector<LSInfo>&                      getLSInfo() ;
    VolumeKernel*                               getMesh() const ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;
    short                                       getSign(const long &) const;
    double                                      getSizeNarrowBand() const;
    bool                                        isInNarrowBand(const long &) const;

    void                                        setSizeNarrowBand(double) ;

    virtual double                              computeSizeNarrowBand( LevelSetObject * )=0;
    virtual double                              updateSizeNarrowBand( const std::vector<adaption::Info> & )=0;

    void                                        clearAfterAdaption( const std::vector<adaption::Info> &, double & ) ;

    void                                        propagateSign( LevelSetObject *) ;
    void                                        propagateValue( LevelSetObject *) ;

    void                                        dump( std::fstream &);
    void                                        restore( std::fstream &);

    protected:
    void                                        solveEikonal( double, double );
    virtual double                              updateEikonal( double, double, const long & ) ;

    std::array<double,3>                        computeGradientUpwind(const long &) ;
    std::array<double,3>                        computeGradientCentral(const long &) ;


# if BITPIT_ENABLE_MPI
    MPI_Comm                                    getCommunicator() ;
    void                                        finalizeMPI() ;
    bool                                        assureMPI() ;
    void                                        writeCommunicationBuffer( const std::vector<long> &, OBinaryStream &, OBinaryStream & );
    void                                        readCommunicationBuffer( const long &, IBinaryStream & ) ;
# endif

};

class LevelSetCartesian : public LevelSetKernel{

    private:
    VolCartesian*                               m_cartesian ;       /**< Pointer to underlying cartesian mesh*/

    private:
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> & );
    double                                      updateEikonal( double, double, const long & ) ; 

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( VolCartesian & );

    double                                      computeSizeNarrowBand( LevelSetObject * );
    void                                        update( LevelSetObject *, const std::vector<adaption::Info> & ) ;
};

class LevelSetOctree : public LevelSetKernel{

    friend  LevelSetObject ;
    friend  LevelSetSegmentation ;

    private:
    VolOctree*                                  m_octree ;       /**< Pointer to underlying octree mesh*/

    private:
    double                                      computeRSearchFromLevel( uint8_t ) ;
    double                                      computeSizeFromRSearch( double ) ;

    public:
    virtual ~LevelSetOctree();
    LevelSetOctree( VolOctree & );

    double                                      computeSizeNarrowBand( LevelSetObject * );
    double                                      updateSizeNarrowBand( const std::vector<adaption::Info> & );
    void                                        update( LevelSetObject *, const std::vector<adaption::Info> & ) ;

};

class LevelSetObject{

    friend  LevelSet ;
    friend  LevelSetKernel ;
    //friend  LevelSetCartesian;
    //friend  LevelSetOctree ;

    private:
    int                                         m_id;           /**< identifier of object */

    public:
    virtual ~LevelSetObject();
    LevelSetObject(int);

    virtual int                                 getId() const ;
    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;
    virtual LevelSetObject*                     clone() const = 0;

    protected:
    virtual void                                computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &)=0 ;
    virtual void                                updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool &)=0 ;

    virtual void                                dumpDerived( std::fstream &) =0 ;
    virtual void                                restoreDerived( std::fstream &) =0 ;

    virtual void                                seedSign( LevelSetKernel *, long &, double & )const =0;

    void                                        dump( std::fstream &) ;
    void                                        restore( std::fstream &) ;

# if BITPIT_ENABLE_MPI
    virtual void                                finalizeMPI() ;
    virtual void                                writeCommunicationBuffer( const std::vector<long> &, OBinaryStream &, OBinaryStream & ) =0 ;
    virtual void                                readCommunicationBuffer( const long &, IBinaryStream & ) =0 ;
# endif 

};

class LevelSetSegmentation : public LevelSetObject {

    private:
    struct SegInfo{
        std::unordered_set<long>                m_segments ;
        long                                    m_support ;
        bool                                    m_checked ;

        SegInfo( ) ;
        SegInfo( const std::unordered_set<long> & ) ;
        SegInfo( const std::unordered_set<long> &, const long & ) ;
    };

    int                                         m_dimension ;               /**< number of space dimensions */

    SurfUnstructured*                           m_segmentation;             /**< surface segmentation */
    std::unordered_map< long, std::vector< std::array<double,3>> > m_vertexNormal;            /**< vertex normals */
    PiercedVector<SegInfo>                      m_seg;                      /**< cell -> segment association information */


    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int, SurfUnstructured*);
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    const std::unordered_set<long> &            getSimplexList(const long &) ;
    const long &                                getSupportSimplex(const long &) ;
    bool                                        isInNarrowBand( const long &) ;

    void                                        dumpDerived( std::fstream &) ;
    void                                        restoreDerived( std::fstream &) ;

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;

    protected:
    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;
    void                                        lsFromSimplex( LevelSetKernel *, const double &, const bool &, bool filter = false ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;
    void                                        seedSign( LevelSetKernel *, long &, double &) const ;

    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &);
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool & ) ;
    void                                        updateSimplexToCell( LevelSetOctree *, const std::vector<adaption::Info> &, const double & ) ;

    void                                        associateSimplexToCell( LevelSetCartesian *, const double &);
    void                                        associateSimplexToCell( LevelSetOctree *, const double &);

# if BITPIT_ENABLE_MPI
    void                                        writeCommunicationBuffer( const std::vector<long> &, OBinaryStream &, OBinaryStream & ) ;
    void                                        readCommunicationBuffer( const long &, IBinaryStream & ) ;
# endif

};

}

#endif /* LEVELSET_HPP */
