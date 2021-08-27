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

# ifndef __BITPIT_LEVELSET_CACHED_HPP__
# define __BITPIT_LEVELSET_CACHED_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_set>

# include "bitpit_containers.hpp"

# include "levelSetBoundedObject.hpp"
# include "levelSetKernel.hpp"
# include "levelSetSignedObject.hpp"
# include "levelSetStorage.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class SendBuffer;
class RecvBuffer;

class LevelSetObject;

template<typename storage_manager_t>
class LevelSetNarrowBandCacheBase : public virtual storage_manager_t
{
    public:
    using typename storage_manager_t::Kernel;
    using typename storage_manager_t::KernelIterator;

    template<typename T>
    using Storage = typename storage_manager_t::template Storage<T>;

    bool isDirty() const = delete;
    void setDirty(bool dirty) = delete;

    virtual double &                            getValue(const KernelIterator &itr) = 0;
    virtual double                              getValue(const KernelIterator &itr) const = 0;

    virtual std::array<double, 3> &             getGradient(const KernelIterator &itr) = 0;
    virtual const std::array<double, 3> &       getGradient(const KernelIterator &itr) const = 0;

    void                                        set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient);

    void                                        swap(LevelSetNarrowBandCacheBase &other) noexcept;

    protected:
    Storage<double>                            *m_values;    /** Levelset values of the cells inside the narrow band */
    Storage<std::array<double, 3>>             *m_gradients; /** Levelset gradient of the cells inside the narrow band */

    LevelSetNarrowBandCacheBase();

};

template<typename storage_manager_t>
class LevelSetNarrowBandCache : public virtual storage_manager_t, public virtual LevelSetNarrowBandCacheBase<storage_manager_t>
{

};

template<>
class LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager> : public virtual LevelSetExternalPiercedStorageManager, public virtual LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>
{

public:
    LevelSetNarrowBandCache(Kernel *kernel);

    KernelIterator insert(long id, bool sync = true) override;
    void erase(long id, bool sync = true) override;

    bool contains(long id) const override;

    KernelIterator find(long id) const override;
    KernelIterator rawFind(std::size_t) const override;

    double &                            getValue(const KernelIterator &itr) override;
    double                              getValue(const KernelIterator &itr) const override;

    std::array<double, 3> &             getGradient(const KernelIterator &itr) override;
    const std::array<double, 3> &       getGradient(const KernelIterator &itr) const override;

    void swap(LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager> &other) noexcept;

protected:
    Storage<char> *m_narrowBandFlag; //! Flag that defines if the entry is inside the narrow band

    void clearKernel() override;

    void dumpKernel(std::ostream &stream) override;
    void restoreKernel(std::istream &stream) override;

};

template<>
class LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager> : public virtual LevelSetInternalPiercedStorageManager, public virtual LevelSetNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>
{

public:
    LevelSetNarrowBandCache();

    double &                            getValue(const KernelIterator &itr) override;
    double                              getValue(const KernelIterator &itr) const override;

    std::array<double, 3> &             getGradient(const KernelIterator &itr) override;
    const std::array<double, 3> &       getGradient(const KernelIterator &itr) const override;

    void swap(LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager> &other) noexcept;

};

template<>
class LevelSetNarrowBandCache<LevelSetDirectStorageManager> : public virtual LevelSetDirectStorageManager, public virtual LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>
{

public:
    LevelSetNarrowBandCache(std::size_t nItems);

    KernelIterator insert(long id, bool sync = true) override;
    void erase(long id, bool sync = true) override;

    bool contains(long id) const override;

    KernelIterator find(long id) const override;
    KernelIterator rawFind(std::size_t) const override;

    double &                            getValue(const KernelIterator &itr) override;
    double                              getValue(const KernelIterator &itr) const override;

    std::array<double, 3> &             getGradient(const KernelIterator &itr) override;
    const std::array<double, 3> &       getGradient(const KernelIterator &itr) const override;

    void swap(LevelSetNarrowBandCache<LevelSetDirectStorageManager> &other) noexcept;

protected:
    Storage<char> *m_narrowBandFlag; //! Flag that defines if the entry is inside the narrow band

    void clearKernel() override;

    void dumpKernel(std::ostream &stream) override;
    void restoreKernel(std::istream &stream) override;

};

template<typename narrow_band_cache_t>
class LevelSetCachedObjectInterface : public virtual LevelSetObjectInterface {

public:
    narrow_band_cache_t * initializeNarrowBandCache();

    narrow_band_cache_t * getNarrowBandCache();
    const narrow_band_cache_t * getNarrowBandCache() const;

    void clearNarrowBandCache();

    void dumpNarrowBandCache(std::ostream &stream);
    void restoreNarrowBandCache(std::istream &stream);

    bool isInNarrowBand(long id) const override;

    void swap(LevelSetCachedObjectInterface &other) noexcept;

protected:
    std::shared_ptr<narrow_band_cache_t> m_narrowBandCache; //! Narrow band cache

    LevelSetCachedObjectInterface() = default;

};

template<typename narrow_band_cache_t>
class LevelSetNarrowBandCacheFactory
{

public:
    static std::shared_ptr<narrow_band_cache_t> create(LevelSetCachedObjectInterface<narrow_band_cache_t> *object);

};

template<>
class LevelSetNarrowBandCacheFactory<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>
{

public:
    static std::shared_ptr<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>> create(LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>> *object);

};

template<>
class LevelSetNarrowBandCacheFactory<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>
{

public:
    static std::shared_ptr<LevelSetNarrowBandCache<LevelSetDirectStorageManager>> create(LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetDirectStorageManager>> *object);

};

template<typename narrow_band_cache_t>
class LevelSetCachedObject : public LevelSetObject, public LevelSetCachedObjectInterface<narrow_band_cache_t>, public LevelSetSignedObjectInterface {

    protected:
    void                                        _clear( ) override ;

    void                                        _clearAfterMeshAdaption(const std::vector<adaption::Info> & ) override ;

    void                                        _dump( std::ostream &) override ;
    void                                        _restore( std::istream &) override ;

# if BITPIT_ENABLE_MPI
    void                                        _writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) override ;
    void                                        _readCommunicationBuffer( const std::vector<long> &, RecvBuffer & )  override;
# endif 

    std::shared_ptr<LevelSetSignStorage>        createSignStorage() override;

    public:
    LevelSetCachedObject(int);

    void                                        setKernel(LevelSetKernel *) override ;

    LevelSetInfo                                getLevelSetInfo(long ) const override ;
    double                                      getLS(long ) const override ;
    double                                      getValue(long ) const override ;
    short                                       getSign(long ) const override ;
    std::array<double,3>                        getGradient(long ) const override ;

};

}

// Include template implementations
# include "levelSetCachedObject.tpp"

// Explicit instantization
#ifndef __BITPIT_LEVELSET_CACHED_SRC__
namespace bitpit {

extern template class LevelSetNarrowBandCacheBase<LevelSetExternalPiercedStorageManager>;
extern template class LevelSetNarrowBandCacheBase<LevelSetInternalPiercedStorageManager>;
extern template class LevelSetNarrowBandCacheBase<LevelSetDirectStorageManager>;

extern template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
extern template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
extern template class LevelSetCachedObjectInterface<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>;

extern template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
extern template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetInternalPiercedStorageManager>>;
extern template class LevelSetCachedObject<LevelSetNarrowBandCache<LevelSetDirectStorageManager>>;

}
#endif

#endif 
