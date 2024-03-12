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

template<typename SourceLevelSetObject, typename BaseLevelSetObject>
friend class LevelSetProxyObject;

public:
    typedef LevelSetCachedKernel::CellCacheCollection CellCacheCollection;

    virtual ~LevelSetObject();

    virtual LevelSetObject*                     clone() const =0;

    virtual const LevelSetKernel *              getKernel() const;

    virtual LevelSetFieldset                    getSupportedFields() const;

    int                                         getId() const ;
    virtual bool                                isPrimary() const ;

    virtual bool                                empty() const = 0;

    std::size_t                                 getReferenceCount() const ;

    void                                        update(const std::vector<adaption::Info> &adaptionData);

    LevelSetBulkEvaluationMode                  getCellBulkEvaluationMode() const;
    void                                        setCellBulkEvaluationMode(LevelSetBulkEvaluationMode evaluationMode);

    LevelSetCacheMode                           getFieldCellCacheMode(LevelSetField field) const;
    void                                        enableFieldCellCache(LevelSetField field, LevelSetCacheMode cacheMode);
    void                                        disableFieldCellCache(LevelSetField field);

    double                                      getNarrowBandSize() const;
    virtual bool                                isCellInNarrowBand(long id) const;
    virtual bool                                isInNarrowBand(const std::array<double,3> &point) const;

    LevelSetIntersectionStatus                  intersectSurface(long, LevelSetIntersectionMode=LevelSetIntersectionMode::FAST_FUZZY) const;

    virtual short                               evalCellSign(long id) const;
    virtual double                              evalCellValue(long id, bool signedLevelSet) const;
    virtual std::array<double,3>                evalCellGradient(long id, bool signedLevelSet) const;
    virtual std::array<double,3>                evalCellProjectionPoint(long id) const;

    virtual short                               evalSign(const std::array<double,3> &point) const;
    virtual double                              evalValue(const std::array<double,3> &point, bool signedLevelSet) const;
    virtual std::array<double,3>                evalGradient(const std::array<double,3> &point, bool signedLevelSet) const;
    virtual std::array<double,3>                evalProjectionPoint(const std::array<double,3> &point) const;

    void                                        enableVTKOutput(const LevelSetFieldset &fieldset, bool enable=true);
    void                                        enableVTKOutput(const LevelSetFieldset &fieldset, const std::string &objectName, bool enable=true);
    void                                        enableVTKOutput(LevelSetField field, bool enable=true);
    void                                        enableVTKOutput(LevelSetField field, const std::string &objectName, bool enable=true);
    void                                        enableVTKOutput(LevelSetWriteField field, bool enable=true);
    void                                        enableVTKOutput(LevelSetWriteField fieldset, const std::string &objectName, bool enable=true);
    void                                        flushData(std::fstream &, const std::string &, VTKFormat) override;

    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   computeProjectionPoint(long cellId) const);
    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   computeVertexProjectionPoint(long vertexId) const);

    BITPIT_DEPRECATED(std::array<double BITPIT_COMMA 3>   computeProjectionPoint(const std::array<double,3> &point) const);

    BITPIT_DEPRECATED_FOR(short                               getSign(long cellId) const, short evalCellSign(long id) const);
    BITPIT_DEPRECATED_FOR(double                              getValue(long cellId) const, double evalCellValue(long id, bool signedLevelSet) const);
    BITPIT_DEPRECATED_FOR(std::array<double BITPIT_COMMA 3>   getGradient(long cellId) const, std::array<double BITPIT_COMMA 3> evalCellGradient(long id, bool signedLevelSet) const);

    BITPIT_DEPRECATED(LevelSetInfo                        getLevelSetInfo(long cellId) const);
    BITPIT_DEPRECATED(double                              getLS(long cellId) const);

    BITPIT_DEPRECATED(double                              getSizeNarrowBand() const);

protected:
    template<typename data_t>
    using CellCacheEntry = typename CellCacheCollection::ValueCache<data_t>::Entry;

    static const bool                           CELL_CACHE_IS_SIGNED;
    static const LevelSetIntersectionMode       CELL_LOCATION_INTERSECTION_MODE;

    LevelSetKernel*                             m_kernel;           /**< Levelset kernel */

    bool                                        m_defaultSignedLevelSet;

    LevelSetFieldMap<std::string>               m_enabledOutputFields;

    double                                      m_narrowBandSize; //!< Size of narrow band

    std::size_t                                 m_cellLocationCacheId; //!< Id of the cache that will keep track if cell zones
    std::size_t                                 m_cellPropagatedSignCacheId; //!< Id of the cache that will keep track if cell propagated sign

    LevelSetObject(int);
    LevelSetObject(const LevelSetObject &other);
    LevelSetObject(LevelSetObject &&other);

    void                                        setDefaultLevelSetSigndness(bool signedLevelSet);

    virtual void                                setKernel(LevelSetKernel *);
    virtual LevelSetKernel *                    getKernel();

    void                                        evaluate();
    void                                        update();

    LevelSetZone                                getCellZone(long id) const;
    LevelSetCellLocation                        getCellLocation(long id) const;

    void                                        setNarrowBandSize(double size);
    void                                        evaluateCellNarrowBandData();
    void                                        updateCellNarrowBandData(const std::vector<adaption::Info> &adaptionData);
    void                                        destroyCellNarrowBandData();

    virtual void                                fillCellLocationCache();
    virtual void                                fillCellLocationCache(const std::vector<adaption::Info> &adaptionData);
    virtual LevelSetCellLocation                fillCellGeometricNarrowBandLocationCache(long id);
    virtual std::size_t                         createCellLocationCache(std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID);
    void                                        destroyCellLocationCache();

    void                                        evaluateCellBulkData();
    void                                        updateCellBulkData(const std::vector<adaption::Info> &adaptionData);
    void                                        destroyCellBulkData();

    virtual void                                fillCellPropagatedSignCache();
    virtual std::size_t                         createCellPropagatedSignCache(std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID);
    void                                        destroyCellPropagatedSignCache();

    virtual void                                dump(std::ostream &);
    virtual void                                restore(std::istream &);

    virtual LevelSetIntersectionStatus          _intersectSurface(long, double distance, LevelSetIntersectionMode=LevelSetIntersectionMode::FAST_FUZZY) const;

    virtual short                               _evalCellSign(long id) const = 0;
    virtual double                              _evalCellValue(long id, bool signedLevelSet) const = 0;
    virtual std::array<double,3>                _evalCellGradient(long id, bool signedLevelSet) const = 0;

    virtual short                               _evalSign(const std::array<double,3> &point) const;
    virtual double                              _evalValue(const std::array<double,3> &point, bool signedLevelSet) const = 0;
    virtual std::array<double,3>                _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const = 0;

    short                                       evalValueSign(double value) const;

# if BITPIT_ENABLE_MPI
    void                                        startCellCacheExchange( const std::unordered_map<int,std::vector<long>> &recvCellIds, std::size_t cacheIds, DataCommunicator * ) const;
    void                                        startCellCachesExchange( const std::unordered_map<int,std::vector<long>> &recvCellIds, const std::vector<std::size_t> &cacheIds, DataCommunicator * ) const;
    void                                        completeCellCacheExchange( const std::unordered_map<int,std::vector<long>> &sendCellIds, std::size_t cacheIds, DataCommunicator * );
    void                                        completeCellCachesExchange( const std::unordered_map<int,std::vector<long>> &sendCellIds, const std::vector<std::size_t> &cacheIds, DataCommunicator * );
# endif

    void                                        adaptCellCaches(const std::vector<adaption::Info> &adaptionData);

    void                                        clearCellCache(std::size_t cacheId, bool release);
    void                                        pruneCellCache(std::size_t cacheId, const std::vector<long> &cellIds);

    std::vector<long>                           evalCellCacheFillIds(LevelSetZone zone, LevelSetCacheMode cacheMode) const;
    std::vector<long>                           evalCellCacheFillIds(LevelSetZone zone, LevelSetCacheMode cacheMode, const std::vector<adaption::Info> &adaptionData) const;
    std::vector<long>                           evalCellOnDemandCacheFillIds(LevelSetZone zone) const;
    std::vector<long>                           evalCellOnDemandCacheFillIds(LevelSetZone zone, const std::vector<adaption::Info> &adaptionData) const;
    std::vector<long>                           evalCellNarrowBandCacheFillIds(LevelSetZone zone) const;
    std::vector<long>                           evalCellNarrowBandCacheFillIds(LevelSetZone zone,   const std::vector<adaption::Info> &adaptionData) const;
    std::vector<long>                           evalCellFullCacheFillIds(LevelSetZone zone) const;
    std::vector<long>                           evalCellFullCacheFillIds(LevelSetZone zone, const std::vector<adaption::Info> &adaptionData) const;

    std::vector<long>                           evalCellCacheStaleIds(const std::vector<adaption::Info> &adaptionData) const;

    template<typename value_t, typename evaluator_t, typename fallback_t>
    value_t                                     evalCellFieldCached(LevelSetField field, long id, const evaluator_t &evaluator, const fallback_t &fallback) const;
    template<typename value_t, typename evaluator_t, typename fallback_t>
    value_t                                     evalCellField(LevelSetField field, long id, const evaluator_t &evaluator, const fallback_t &fallback) const;

    void                                        fillFieldCellCaches(LevelSetZone zone, const std::vector<LevelSetField> &fields);
    void                                        fillFieldCellCaches(LevelSetZone zone, const std::vector<LevelSetField> &fields, const std::vector<adaption::Info> &adaptionData);
    void                                        fillFieldCellCache(LevelSetField field, const std::vector<long> &cellIds);
    virtual void                                fillFieldCellCache(LevelSetField field, long id);
    template<typename value_t>
    void                                        fillFieldCellCache(LevelSetField field, long id, const value_t &value) const;

    template<typename value_t>
    CellCacheCollection::ValueCache<value_t> *  getFieldCellCache(LevelSetField field) const;
    CellCacheCollection::Cache *                getFieldCellCache(LevelSetField field) const;
    std::size_t                                 getFieldCellCacheId(LevelSetField field) const;
    virtual std::size_t                         createFieldCellCache(LevelSetField field, std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID);
    template<typename value_t>
    std::size_t                                 createFieldCellCache(LevelSetField field, std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID);
    virtual void                                destroyFieldCellCache(LevelSetField field);



    template<typename value_t>
    CellCacheCollection::ValueCache<value_t> *  getCellCache(std::size_t cacheId) const;
    CellCacheCollection::Cache *                getCellCache(std::size_t cacheId) const;
    template<typename value_t>
    std::size_t                                 createCellCache(LevelSetFillIn expectedFillIn, std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID);
    void                                        destroyCellCache(std::size_t cacheId);

    bool                                        hasVTKOutputData(LevelSetField field, const std::string &objectName) const;
    void                                        removeVTKOutputData(LevelSetField field, const std::string &objectName);
    virtual void                                addVTKOutputData(LevelSetField field, const std::string &objectName);
    std::string                                 getVTKOutputDataName(LevelSetField field, const std::string &objectName) const;
    virtual std::string                         getVTKOutputFieldName(LevelSetField field) const;
    virtual void                                flushVTKOutputData(std::fstream &stream, VTKFormat format,
                                                                   LevelSetField field) const;
    template<typename value_t, typename evaluator_t, typename fallback_t>
    void                                        flushVTKOutputData(std::fstream &stream, VTKFormat format, LevelSetField field,
                                                                   const evaluator_t evluator, const fallback_t fallback) const;

    BITPIT_DEPRECATED(void                      setSizeNarrowBand(double));

private:
    int                                          m_id;           /**< identifier of object */

    std::size_t                                  m_nReferences;

    LevelSetBulkEvaluationMode                   m_cellBulkEvaluationMode; //!< Evaluation mode for cell data in the bulk

    mutable std::unique_ptr<CellCacheCollection> m_cellCacheCollection; //!< Cell cache collection

    std::vector<LevelSetCacheMode>               m_cellFieldCacheModes; //!< Mode of the cell cache for the fields
    std::vector<std::size_t>                     m_cellFieldCacheIds; //!< Ids of the field cell caches

    void                                         setId(int id);

    std::size_t                                  incrementReferenceCount();
    std::size_t                                  decrementReferenceCount();

};

}

// Include template implementations
#include "levelSetObject.tpp"

#endif
