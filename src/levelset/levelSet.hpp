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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

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

class LSObject ;
class LevelSetSegmentation ;

namespace levelSetDefaults{
    const double                            VALUE = 1.e18 ;             /**< Default value for levelset function */
    const std::array<double,3>              GRADIENT = {{0.,0.,0.}};    /**< Default value for levelset gradient */
    const short                             SIGN = 1;                   /**< Default value for the sign */
    const int                               OBJECT = -1 ;               /**< Default value for closest object id */
    const std::set<long>                    LIST = { } ;                /**< Default value for closest segment */
    const long                              ELEMENT = -1 ;              /**< Default value for segmments in narrow band */
};

class LevelSet{

    friend  LSObject ;
    friend  LevelSetSegmentation ;

//    public:
//    static double                               NULL_VALUE ;    /**< Default value for levelset function */
//    static std::array<double,3>                 NULL_GRADIENT ; /**< Default value for levelset gradient */
//    static int                                  NULL_OBJECT ;   /**< Default value for closest object id */

    protected:
    struct LSInfo{
        double                                  value ;         /**< Levelset value */
        std::array<double,3>                    gradient ;      /**< Levelset gradient */
        int                                     object ;        /**< Id of closest object */
        short                                   active ;        /**< Flag for nodes for fast marching */
    };

    PiercedVector<LSInfo>                       info ;          /**< Levelset information for each cell */

    double                                      RSearch;        /**< Size of narrow band */
    bool                                        m_userRSearch;  /**< Flag if user has set size of narrow band (default=false)  */

    bool                                        signedDF;       /**< Flag for sigend/unsigned distance function (default = true) */
    bool                                        propagateS;     /**< Flag for sign propagation from narrow band (default = false) */
    bool                                        propagateV;     /**< Flag for value propagation from narrow band (default = false) */

    VolumeKernel*                               m_mesh ;        /**< Pointer to underlying mesh*/

    protected:
    LevelSet() ;
    LevelSet( VolumeKernel *) ;

    virtual void                                computeSizeNarrowBand( LSObject * )=0;

    virtual double                              updateSizeNarrowBand( std::vector<Adaption::Info> & )=0;
    void                                        clearAfterAdaption( std::vector<Adaption::Info> &, double & ) ;

    void                                        propagateSign( LSObject *) ;
    void                                        propagateValue( LSObject *) ;

    void                                        solveEikonal( double, double );
    virtual double                              updateEikonal( double, double, const long & ) ;

    std::array<double,3>                        computeGradientUpwind(const long &) ;
    std::array<double,3>                        computeGradientCentral(const long &) ;


    public:
    virtual ~LevelSet() ;

    double                                      getLS(const long &) const;
    std::array<double,3>                        getGradient(const long &) const ;

    short                                       getSign(const long &) const;
    double                                      getSizeNarrowBand() const;

    void                                        setSizeNarrowBand(double) ;
    void                                        setSign(bool);
    void                                        setPropagateSign(bool) ;
    void                                        setPropagateValue(bool) ;

    bool                                        isInNarrowBand(const long &) const;

    virtual void                                compute( LSObject * )=0 ;
    virtual void                                update( LSObject *, std::vector<Adaption::Info> & )=0 ;

};

class LevelSetCartesian : public LevelSet{

    friend  LSObject ;
    friend  LevelSetSegmentation ;

    private:
    VolCartesian*                               m_cmesh ;       /**< Pointer to underlying cartesian mesh*/

    private:
    void                                        computeSizeNarrowBand( LSObject * );
    double                                      updateSizeNarrowBand( std::vector<Adaption::Info> & );
    double                                      updateEikonal( double, double, const long & ) ; 

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( VolCartesian & );

    void                                        compute( LSObject * ) ;
    void                                        update( LSObject *, std::vector<Adaption::Info> & ) ;
};

class LevelSetOctree : public LevelSet{

    friend  LSObject ;
    friend  LevelSetSegmentation ;

    private:
    VolOctree*                                  m_omesh ;       /**< Pointer to underlying octree mesh*/

    private:
    void                                        computeSizeNarrowBand( LSObject * );
    double                                      updateSizeNarrowBand( std::vector<Adaption::Info> & );

    double                                      computeRSearchFromLevel( uint8_t ) ;
    int                                         computeLevelFromRSearch( double ) ;

    public:
    virtual ~LevelSetOctree();
    LevelSetOctree( VolOctree & );

    void                                        compute( LSObject * ) ;
    void                                        update( LSObject *, std::vector<Adaption::Info> & ) ;

};

class LSObject{

    friend  LevelSet ;
    friend  LevelSetCartesian;
    friend  LevelSetOctree ;

    private:
    int                                         m_id;           /**< identifier of object */

    protected:
    virtual ~LSObject();
    LSObject(int);
    virtual LSObject*                           clone() const = 0;

    virtual int                                 getId() const ;

    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;
    virtual void                                seedSign( LevelSet *, long &, double & )const =0;

    virtual void                                computeLSInNarrowBand( LevelSetCartesian *)=0 ;
    virtual void                                computeLSInNarrowBand( LevelSetOctree *)=0 ;

    virtual void                                updateLSInNarrowBand( LevelSetOctree *, std::vector<Adaption::Info> &, double &)=0 ;
};

class LevelSetSegmentation : public LSObject {

    public:
//    static std::set<long>                       NULL_LIST ;     /**< Default value for closest segment */
//    static long                                 NULL_ELEMENT ;  /**< Default value for segmments in narrow band */

    private:
    struct SegData{
        std::set<long>                          m_segments ;
        long                                    m_support ;
        bool                                    m_checked ;

        SegData( ) ;
        SegData( const std::set<long> & ) ;
        SegData( const std::set<long> &, const long & ) ;

    };

    protected:
    SurfUnstructured                            *stl;           /**< surface Triangulation */
    double                                      abs_tol;        /**< tolerance used for calculating normals */

    PiercedVector<SegData>                      m_segInfo ;     /**< cell segment association information */


    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int, SurfUnstructured*);
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    const std::set<long> &                      getSimplexList(const long &) ;
    const long &                                getSupportSimplex(const long &) ;
    bool                                        isInNarrowBand( const long &) ;

    protected:
    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;
    void                                        lsFromSimplex( LevelSet *, const double &, bool filter = false ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    bool                                        seedNarrowBand( LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;
    void                                        seedSign( LevelSet *, long &, double &) const ;

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;

    void                                        computeLSInNarrowBand( LevelSetCartesian *);
    void                                        associateSimplexToCell( LevelSetCartesian *);

    void                                        computeLSInNarrowBand( LevelSetOctree *);
    void                                        associateSimplexToCell( LevelSetOctree *);
    void                                        updateLSInNarrowBand( LevelSetOctree *, std::vector<Adaption::Info> &, double & ) ;
    void                                        updateSimplexToCell( LevelSetOctree *, std::vector<Adaption::Info> &, double & ) ;

};

}

#endif /* LEVELSET_HPP */
