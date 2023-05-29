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

# ifndef __BITPIT_LEVELSET_HPP__
# define __BITPIT_LEVELSET_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_map>
# include <unordered_set>
# include <memory>

# include "levelSetCommon.hpp"

# include "bitpit_IO.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}
class VolumeKernel;
class SurfaceKernel;
class SurfUnstructured;

class LevelSetKernel;
class LevelSetObject;

class LevelSet{

private:
    std::unique_ptr<LevelSetKernel>                             m_kernel ;            /**< LevelSet computational kernel */

    LevelSetFillIn          m_expectedFillIn;          /**< Expected fill-in for data structures */

    bool                    m_signedDistance;          /**< Flag for signed/unsigned distance (default = true) */
    bool                    m_propagateSign;           /**< Flag for sign propagation from narrow band (default = false) */

    IndexGenerator<long>                                        m_objectIdentifierGenerator; /**< Object identifier generator */
    std::unordered_map<int,std::unique_ptr<LevelSetObject>>     m_objects ;                  /**< Objects defining the boundaries */
    std::vector<int>                                            m_objectsProcessingOrder ;   /**< Processing order of objects */

    std::unordered_set<int>                                     m_cacheFilledObjects ;       /**< Objects whose cache has been filled */

    int                     registerObject( std::unique_ptr<LevelSetObject> && ) ;
    bool                    unregisterObject(int id, bool force);

    void                    setObjectProcessingOrder(int) ;
    void                    unsetObjectProcessingOrder(int) ;

    void                    incrementObjectsReferenceCount(int parentId) ;
    void                    decrementObjectsReferenceCount(int parentId) ;

public:
    LevelSet(LevelSetFillIn expectedFillIn = LevelSetFillIn::SPARSE) ;

    LevelSet(LevelSet&& other) = default;

    void                    clear();

    void                    setMesh( VolumeKernel* mesh ) ;

    int                     addObject( std::unique_ptr<SurfaceKernel> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( SurfaceKernel *, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( std::unique_ptr<SurfUnstructured> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( SurfUnstructured *, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( LevelSetBooleanOperation, int, int, int id=levelSetDefaults::OBJECT ) ;
    int                     addObject( LevelSetBooleanOperation, const std::vector<int> &, int id=levelSetDefaults::OBJECT ) ;
    int                     addObject( const std::unordered_set<long> &, int id=levelSetDefaults::OBJECT ) ;
    int                     addObject( const std::vector<long> &, long, bool, int id=levelSetDefaults::OBJECT ) ;
    int                     addObject( std::unique_ptr<LevelSetObject> && ) ;

    void                    removeObjects();
    bool                    removeObject(int);

    bool                    isObjectRemovable(int);

    LevelSetObject &  getObject( int ) const ;
    LevelSetObject *  getObjectPtr( int ) const ;
    std::vector<LevelSetObject *> getObjectPtrs( ) const ;

    template<typename T>
    T &                     getObject( int ) const ;
    template<typename T>
    T *                     getObjectPtr( int ) const ;
    template<typename T>
    std::vector<T *>        getObjectPtrs( ) const ;

    int                     getObjectCount( ) const ;
    std::vector<int>        getObjectIds( ) const ;

    void                    setSign(bool);

    void                    fillCache() ;
    void                    fillCache(int id);
    void                    fillCache(const std::vector<int> &ids);
    void                    updateCache(const std::vector<adaption::Info> &adaptionData) ;
    void                    clearCache(bool release) ;

    BITPIT_DEPRECATED(void   setSizeNarrowBand(double));
    BITPIT_DEPRECATED(double getSizeNarrowBand() const);

    BITPIT_DEPRECATED(void   setPropagateSign(bool));

    BITPIT_DEPRECATED(void   dump( std::ostream &));
    BITPIT_DEPRECATED(void   restore( std::istream &));

    BITPIT_DEPRECATED(void   compute());
    BITPIT_DEPRECATED(void   compute(int id));
    BITPIT_DEPRECATED(void   compute(const std::vector<int> &ids));

    BITPIT_DEPRECATED(void   update(const std::vector<adaption::Info> &adaptionData));
    BITPIT_DEPRECATED(void   update(const std::vector<adaption::Info> &adaptionData, int id));
    BITPIT_DEPRECATED(void   update(const std::vector<adaption::Info> &adaptionData, const std::vector<int> &ids));

# if BITPIT_ENABLE_MPI
    BITPIT_DEPRECATED(void  partition( const std::vector<adaption::Info> & ));
# endif
};

}

// Template implementation
#include "levelSet.tpp"

#endif
