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

    bool                    m_forceSignPropagation;    /**< Flag for forcing sign propagation from narrow band (default = false) */
    bool                    m_signPropagationEnabled;  /**< Flag for sign propagation from narrow band (default = false) */

    double                  m_narrowBandSize;          /**< Size of the narrow band. Regardless of the specified size, the narrow band
                                                            will always contain the intersected cells and their neighbours */

    IndexGenerator<long>                                        m_objectIdentifierGenerator; /**< Object identifier generator */
    std::unordered_map<int,std::unique_ptr<LevelSetObject>>     m_objects ;                  /**< Objects defining the boundaries */
    std::vector<int>                                            m_orderedObjectsIds ;        /**< Object ids sorted according to the order in which they should be processed */

    int                     registerObject( std::unique_ptr<LevelSetObject> && ) ;
    bool                    unregisterObject(int id, bool force);

    void                    registerObjectId(int) ;
    void                    unregisterObjectId(int) ;

    void                    incrementObjectsReferenceCount(int parentId) ;
    void                    decrementObjectsReferenceCount(int parentId) ;

public:
    LevelSet(LevelSetFillIn expectedFillIn = LevelSetFillIn::SPARSE) ;

    LevelSet(LevelSet&& other) = default;

    void                    clear();
    void                    update(const std::vector<adaption::Info> &adaptionData) ;

    double                  getNarrowBandSize() const;
    void                    setNarrowBandSize(double size = 0);

    void                    setMesh( VolumeKernel* mesh ) ;

    template<typename LevelSetSourceObject = LevelSetObject>
    int                     addObjectComplement( int, int id=levelSetDefaults::OBJECT ) ;

    template<typename LevelSetSourceObject = LevelSetObject>
    int                     addObject( LevelSetBooleanOperation, int, int, int id=levelSetDefaults::OBJECT ) ;
    template<typename LevelSetSourceObject = LevelSetObject>
    int                     addObject( LevelSetBooleanOperation, const std::vector<int> &, int id=levelSetDefaults::OBJECT ) ;

    int                     addObject( std::unique_ptr<SurfaceKernel> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( SurfaceKernel *, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( std::unique_ptr<SurfUnstructured> &&, double, int id = levelSetDefaults::OBJECT ) ;
    int                     addObject( SurfUnstructured *, double, int id = levelSetDefaults::OBJECT ) ;
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

    void                    dump( std::ostream &) const;
    void                    restore( std::istream &);

    BITPIT_DEPRECATED_FOR(void setSizeNarrowBand(double size), void setNarrowBandSize(double size));
    BITPIT_DEPRECATED_FOR(double getSizeNarrowBand() const, double getNarrowBandSize() const);

    BITPIT_DEPRECATED(void   setSign(bool));
    BITPIT_DEPRECATED(void   setPropagateSign(bool));

    BITPIT_DEPRECATED_FOR(void compute(), void evaluate());
    BITPIT_DEPRECATED_FOR(void compute(int id), void evaluate(int id));
    BITPIT_DEPRECATED_FOR(void compute(const std::vector<int> &ids), void evaluate(const std::vector<int> &ids));

    BITPIT_DEPRECATED_FOR(void update(const std::vector<adaption::Info> &adaptionData, int id), void update(const std::vector<adaption::Info> &adaptionData));
    BITPIT_DEPRECATED_FOR(void update(const std::vector<adaption::Info> &adaptionData, const std::vector<int> &ids), void update(const std::vector<adaption::Info> &adaptionData));

# if BITPIT_ENABLE_MPI
    BITPIT_DEPRECATED(void  partition( const std::vector<adaption::Info> & ));
# endif
};

}

// Template implementation
#include "levelSet.tpp"

#endif
