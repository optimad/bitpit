# ifndef __CLASS_LEVELSET_STL_HPP__
# define __CLASS_LEVELSET_STL_HPP__


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

# include "bitpit.hpp"

/*!
 *  @class  LevelSet
 *  @brief  Level Set on Stl Manager Class
 *
 *  LevelSet is a user interface class. One user should (read can...) work only
 *  with this class and its methods to maintain a signed distance function (Sdf)
 *  computed from a piece-wise linear approximation of a d manifold in a 3D Euclidean
 *  space. Sdf is computed in a narrow band of at least 2 mesh cell centers
 *  around the geometry.
 *
 *
 */

namespace bitpit{
class LSObject ;
class LevelSetSegmentation ;


class LevelSet{

    friend  bitpit::LSObject ;
    friend  bitpit::LevelSetSegmentation ;

    protected:
    struct LSInfo{
        double                                  value ;         /**< Levelset value */
        std::array<double,3>                    gradient ;      /**< Levelset gradient */
        int                                     object ;        /**< Id of closest object */
        short                                   active ;        /**< Flag for nodes for fast marching */
    };

    bitpit::PiercedVector<LSInfo>               info ;          /**< Levelset information for each cell */

    double                                      RSearch;        /**< Size of narrow band */

    bool                                        signedDF;       /**< Flag for sigend/unsigned distance function ( default = true) */
    bool                                        propagateS;     /**< Flag for sign propagation from narrow band (default = false) */
    bool                                        propagateV;     /**< Flag for value propagation from narrow band (default = false) */

    bitpit::VolumeKernel*                       m_mesh ;        /**< Pointer to underlying mesh*/

    protected:
    LevelSet() ;
    LevelSet( bitpit::VolumeKernel *) ;

    virtual void                                compute( bitpit::LSObject * )=0 ;
    virtual void                                computeSizeNarrowBand( bitpit::LSObject * )=0;

    virtual void                                update( bitpit::LSObject *, std::vector<bitpit::Adaption::Info> & )=0 ;
    virtual double                              updateSizeNarrowBand( std::vector<bitpit::Adaption::Info> & )=0;
    void                                        clearAfterRefinement( std::vector<bitpit::Adaption::Info> & ) ;

    void                                        propagateSign( bitpit::LSObject *) ;
    void                                        propagateValue( bitpit::LSObject *) ;

    void                                        solveEikonal( double, double );
    virtual double                              updateEikonal( double, double, const long & ) ;

    std::array<double,3>                        computeGradientUpwind(const long &) ;
    std::array<double,3>                        computeGradientCentral(const long &) ;


    public:
    virtual ~LevelSet() ;

    double                                      getLS(const long &);
    std::array<double,3>                        getGradient(const long &) ;

    double                                      getSizeNarrowBand() ;

    void                                        setSizeNarrowBand(double) ;
    void                                        setSign(bool);
    void                                        setPropagateSign(bool) ;
    void                                        setPropagateValue(bool) ;


    bool                                        isInNarrowBand(const long &) ;

};

class LevelSetCartesian : public bitpit::LevelSet{

    friend  bitpit::LSObject ;
    friend  bitpit::LevelSetSegmentation ;

    private:
    bitpit::VolCartesian*                       m_cmesh ;

    private:
    void                                        computeSizeNarrowBand( bitpit::LSObject * );
    double                                      updateSizeNarrowBand( std::vector<bitpit::Adaption::Info> & );
    double                                      updateEikonal( double, double, const long & ) ; 

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( bitpit::VolCartesian & );

    void                                        compute( bitpit::LSObject * ) ;
    void                                        update( bitpit::LSObject *, std::vector<bitpit::Adaption::Info> & ) ;
};

class LevelSetOctree : public bitpit::LevelSet{

    friend  bitpit::LSObject ;
    friend  bitpit::LevelSetSegmentation ;

    private:
    bitpit::VolOctree*                          m_omesh ;

    private:
    void                                        computeSizeNarrowBand( bitpit::LSObject * );
    double                                      updateSizeNarrowBand( std::vector<bitpit::Adaption::Info> & );

    double                                      computeRSearchFromLevel( uint8_t ) ;
    int                                         computeLevelFromRSearch( double ) ;

    public:
    virtual ~LevelSetOctree();
    LevelSetOctree( bitpit::VolOctree & );

    void                                        compute( bitpit::LSObject * ) ;
    void                                        update( bitpit::LSObject *, std::vector<bitpit::Adaption::Info> & ) ;

};

class LSObject{

    friend  bitpit::LevelSet ;
    friend  bitpit::LevelSetCartesian;
    friend  bitpit::LevelSetOctree ;

    private:
    int                                         m_id;

    protected:
    virtual ~LSObject();
    LSObject(int);
    virtual LSObject*                           clone() const = 0;

    virtual int                                 getId() const ;

    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;
    virtual void                                seedSign( bitpit::LevelSet *, long &, double & )const =0;

    virtual void                                computeLSInNarrowBand( bitpit::LevelSetCartesian *)=0 ;
    virtual void                                computeLSInNarrowBand( bitpit::LevelSetOctree *)=0 ;

    virtual void                                updateLSInNarrowBand( bitpit::LevelSetCartesian *, std::vector<bitpit::Adaption::Info> &, double &)=0 ;
    virtual void                                updateLSInNarrowBand( bitpit::LevelSetOctree *, std::vector<bitpit::Adaption::Info> &, double &)=0 ;
};

/*!
 *  @ingroup    LevelSet
 *  @class      LevelSetSegmentation
 *  @brief      LevelSet of surface segmentations
 *
 *  LevelSetSegmentation provides specific methods for calculating distances with repect to one surface segmentations.
 */
class LevelSetSegmentation : public bitpit::LSObject {

    public:
    static std::set<long>                       NULL_LIST ;
    static long                                 NULL_ELEMENT ;

    private:
    struct SegData{
        std::set<long>                          m_segments ;
        long                                    m_support ;

        SegData( ) ;
        SegData( const std::set<long> & ) ;
        SegData( const std::set<long> &, const long & ) ;

    };

    protected:
    bitpit::SurfUnstructured                    *stl;           /**< surface Triangulation */
    double                                      abs_tol;        /**< tolerance used for calculating normals */

    bitpit::PiercedVector<SegData>              m_segInfo ;


    public:
    virtual ~LevelSetSegmentation();
    LevelSetSegmentation(int, bitpit::SurfUnstructured*);
    LevelSetSegmentation(const LevelSetSegmentation&);

    LevelSetSegmentation*                       clone() const ;

    const std::set<long> &                      getSimplexList(const long &) ;
    const long &                                getSupportSimplex(const long &) ;
    bool                                        isInNarrowBand( const long &) ;

    protected:
    std::vector<std::array<double,3>>           getSimplexVertices( const long & ) const ;
    void                                        lsFromSimplex( bitpit::LevelSet *, const double &, bool filter = false ) ;
    void                                        infoFromSimplex(const std::array<double,3> &, const long &, double &, double &, std::array<double,3> &,std::array<double,3> &) const ;

    bool                                        seedNarrowBand( bitpit::LevelSetCartesian *, std::vector<std::array<double,3>> &, std::vector<int> &) ;
    void                                        seedSign( bitpit::LevelSet *, long &, double &) const ;

    void                                        getBoundingBox( std::array<double,3> &, std::array<double,3> & ) const ;

    void                                        computeLSInNarrowBand( bitpit::LevelSetCartesian *);
    void                                        associateSimplexToCell( bitpit::LevelSetCartesian *);
    void                                        updateLSInNarrowBand( bitpit::LevelSetCartesian *, std::vector<bitpit::Adaption::Info> &, double & ) ;
    void                                        updateSimplexToCell( bitpit::LevelSetCartesian *, std::vector<bitpit::Adaption::Info> &, double & ) ;

    void                                        computeLSInNarrowBand( bitpit::LevelSetOctree *);
    void                                        associateSimplexToCell( bitpit::LevelSetOctree *);
    void                                        updateLSInNarrowBand( bitpit::LevelSetOctree *, std::vector<bitpit::Adaption::Info> &, double & ) ;
    void                                        updateSimplexToCell( bitpit::LevelSetOctree *, std::vector<bitpit::Adaption::Info> &, double & ) ;

};

}

#endif /* LEVELSET_HPP */
