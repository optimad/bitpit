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
# include <unordered_map>
# include <memory>

# include "levelSetCommon.hpp"

namespace bitpit{

namespace adaption{
    class Info;
}
class VolumeKernel;
class VolCartesian;
class VolOctree;
class SurfaceKernel;
class SurfUnstructured;

class LevelSetKernel;
class LevelSetObject;

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

}

#endif
