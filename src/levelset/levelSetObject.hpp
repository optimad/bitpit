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

# ifndef __BITPIT_LEVELSET_OBJECT_HPP__
# define __BITPIT_LEVELSET_OBJECT_HPP__

// Standard Template Library
# include <iostream>
# include <array>
# include <vector>
# include <unordered_map>
# include <unordered_set>

# include "bitpit_IO.hpp"
# if BITPIT_ENABLE_MPI
# include "bitpit_communications.hpp"
# endif
# include "bitpit_containers.hpp"
# include "levelSetCommon.hpp"
# include "levelSetCache.hpp"
# include "levelSetKernel.hpp"
# include "levelSetCartesianKernel.hpp"
# include "levelSetOctreeKernel.hpp"
# include "levelSetUnstructuredKernel.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}
class SendBuffer;
class RecvBuffer;

class LevelSetObject : public VTKBaseStreamer {

friend class LevelSet;

public:
    typedef LevelSetCachedKernel::CellCacheKey CellCacheKey;
    typedef LevelSetCachedKernel::CellCache CellCache;
    typedef LevelSetCachedKernel::CellCacheCollection CellCacheCollection;

    template<typename value_t>
    using CellValueCache = LevelSetCachedKernel::CellValueCache<value_t>;

    virtual ~LevelSetObject();

    virtual LevelSetObject*                     clone() const =0;

    virtual const LevelSetKernel *              getKernel() const;

    void                                        setFieldCache(const LevelSetFieldset &fieldset, LevelSetCacheMode cacheMode);
    virtual void                                setFieldCache(LevelSetField field, LevelSetCacheMode cacheMode);
    void                                        fillCache();
    void                                        updateCache(const std::vector<adaption::Info> &adaptionData);
    void                                        clearCache(bool release = false);

    virtual LevelSetFieldset                    getSupportedFields() const;

    int                                         getId() const ;
    virtual bool                                isPrimary() const ;

    std::size_t                                 getReferenceCount() const ;

    LevelSetIntersectionStatus                  intersectSurface(long, LevelSetIntersectionMode=LevelSetIntersectionMode::FAST_FUZZY) const;

    virtual void                                setNarrowBandSize(double size);
    virtual double                              getNarrowBandSize() const;
    virtual bool                                isCellInNarrowBand(long id) const;
    virtual bool                                isInNarrowBand(const std::array<double,3> &point) const;

    virtual short                               evalCellSign(long id) const;
    virtual double                              evalCellValue(long id, bool signedLevelSet) const;
    virtual std::array<double,3>                evalCellGradient(long id, bool signedLevelSet) const;
    virtual std::array<double,3>                evalCellProjectionPoint(long id) const;

    virtual short                               evalSign(const std::array<double,3> &point) const;
    virtual double                              evalValue(const std::array<double,3> &point, bool signedLevelSet) const;
    virtual std::array<double,3>                evalGradient(const std::array<double,3> &point, bool signedLevelSet) const;
    virtual std::array<double,3>                evalProjectionPoint(const std::array<double,3> &point) const;

    void                                        enableVTKOutput(LevelSetWriteField field, bool enable=true);
    void                                        enableVTKOutput(const LevelSetFieldset &fieldset, bool enable=true);
    void                                        enableVTKOutput(LevelSetField field, bool enable=true);
    void                                        enableVTKOutput(LevelSetWriteField fieldset, const std::string &objectName, bool enable=true);
    void                                        enableVTKOutput(const LevelSetFieldset &fieldset, const std::string &objectName, bool enable=true);
    void                                        enableVTKOutput(LevelSetField field, const std::string &objectName, bool enable=true);
    void                                        flushData(std::fstream &, const std::string &, VTKFormat) override;

    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   computeProjectionPoint(long cellId) const);

    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   computeProjectionPoint(const std::array<double,3> &point) const);

    BITPIT_DEPRECATED(short                               getSign(long cellId) const);
    BITPIT_DEPRECATED(double                              getValue(long cellId) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   getGradient(long cellId) const);

    BITPIT_DEPRECATED(LevelSetInfo                        getLevelSetInfo(long cellId) const);
    BITPIT_DEPRECATED(double                              getLS(long cellId) const);

    BITPIT_DEPRECATED(double                              getSizeNarrowBand() const);

protected:
    LevelSetKernel*                             m_kernel;           /**< Levelset kernel */

    double                                      m_narrowBandSize; //!< Size of narrow band
    std::size_t                                 m_cellNarrowBandCacheId; //!< Id of the cache that will keep track if a cell is inside the narrow band

    bool                                        m_defaultSignedLevelSet;

    LevelSetFieldMap<std::string>               m_enabledOutputFields;

    LevelSetObject(int);
    LevelSetObject(const LevelSetObject &other);
    LevelSetObject(LevelSetObject &&other);

    virtual void                                setKernel(LevelSetKernel *);
    virtual LevelSetKernel *                    getKernel();

    void                                        setDefaultLevelSetSigned(bool signedLevelSet);

    void                                        dump(std::ostream &);
    virtual void                                _dump(std::ostream &);
    void                                        restore(std::istream &);
    virtual void                                _restore( std::istream &);

    virtual bool                                _isCellInNarrowBand(long id, bool checkNeighbours, double *maximumDistance = nullptr) const;

    virtual short                               _evalCellSign(long id) const = 0;
    virtual double                              _evalCellValue(long id, bool signedLevelSet) const = 0;
    virtual std::array<double,3>                _evalCellGradient(long id, bool signedLevelSet) const = 0;

    virtual short                               _evalSign(const std::array<double,3> &point) const = 0;
    virtual double                              _evalValue(const std::array<double,3> &point, bool signedLevelSet) const = 0;
    virtual std::array<double,3>                _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const = 0;

# if BITPIT_ENABLE_MPI
    void                                        exchangeGhosts() ;
    void                                        startExchange( const std::unordered_map<int,std::vector<long>> &, DataCommunicator * );
    void                                        completeExchange( const std::unordered_map<int,std::vector<long>> &, DataCommunicator * );

    virtual void                                writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ;
    virtual void                                _writeCommunicationBuffer(const std::vector<long>&, SendBuffer&)  ;
    virtual void                                readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ;
    virtual void                                _readCommunicationBuffer(const std::vector<long>&, RecvBuffer&)  ;
# endif

    virtual void                                fillFullCellCaches();
    virtual void                                fillFullCellCaches(const std::vector<adaption::Info> &adaptionData);
    virtual void                                fillFullCellSignCache();
    virtual void                                fillFullCellSignCache(const std::vector<adaption::Info> &adaptionData);
    virtual void                                fillNarrowBandCellCaches();
    virtual void                                fillNarrowBandCellCaches(const std::vector<adaption::Info> &adaptionData);

    LevelSetFieldset                            getCachedFields(LevelSetCacheMode cacheMode) const;

    template<typename value_t>
    CellValueCache<value_t> *                   getFieldCellCache(LevelSetField field) const;
    CellCache *                                 getFieldCellCache(LevelSetField field) const;
    LevelSetCacheMode                           getFieldCellCacheMode(LevelSetField field) const;
    template<typename value_t>
    std::size_t                                 registerFieldCellCache(LevelSetField field, LevelSetCacheMode cacheMode);
    void                                        unregisterFieldCellCache(LevelSetField field);
    virtual void                                fillFieldCellCache(long id, LevelSetField field, double searchRadius = std::numeric_limits<double>::max());

    template<typename value_t>
    CellValueCache<value_t> *                   getCellCache(std::size_t cacheId) const;
    CellCache *                                 getCellCache(std::size_t cacheId) const;
    template<typename value_t>
    std::size_t                                 registerCellCache(LevelSetCacheMode cacheMode);
    void                                        unregisterCellCache(std::size_t cacheId);

    virtual void                                flushField(LevelSetField field, std::fstream &stream, VTKFormat format) const;

    BITPIT_DEPRECATED(void                      setSizeNarrowBand(double));

private:
    int                                          m_id;           /**< identifier of object */

    std::size_t                                  m_nReferences;

    mutable std::unique_ptr<CellCacheCollection> m_cellCacheCollection; //!< Cell cache collection
    std::vector<std::size_t>                     m_cellFieldCacheIds; //!< Ids of the cell caches
    std::vector<LevelSetCacheMode>               m_cellFieldCacheModes; //!< Modes of the cell caches

    void                                         setId(int id);

    std::size_t                                  incrementReferenceCount();
    std::size_t                                  decrementReferenceCount();

    void                                         fillNarrowBandCellCaches(LevelSetKernel *levelsetKernel);
    void                                         fillNarrowBandCellCaches(LevelSetKernel *levelsetKernel, const std::vector<adaption::Info> &adaptionData);
    void                                         fillNarrowBandCellCaches(LevelSetCartesianKernel *levelsetKernel);

};

}

// Include template implementations
#include "levelSetObject.tpp"

#endif
