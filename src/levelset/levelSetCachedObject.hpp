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
# include "levelSetSignedObject.hpp"
# include "levelSetStorage.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class SendBuffer;
class RecvBuffer;

class LevelSetObject;

class LevelSetNarrowBandCache : public LevelSetInternalPiercedStorageManager
{
    protected:
    Storage<double>                            *m_values;    /** Levelset values of the cells inside the narrow band */
    Storage<std::array<double, 3>>             *m_gradients; /** Levelset gradient of the cells inside the narrow band */

    public:
    LevelSetNarrowBandCache();

    bool isDirty() const = delete;
    void setDirty(bool dirty) = delete;

    double                                      getValue(const KernelIterator &itr) const;
    const std::array<double, 3> &               getGradient(const KernelIterator &itr) const;

    void                                        set(const KernelIterator &itr, double value, const std::array<double, 3> &gradient);

    void                                        swap(LevelSetNarrowBandCache &other) noexcept;

};

class LevelSetCachedObjectInterface : public virtual LevelSetObjectInterface {

public:
    LevelSetNarrowBandCache * initializeNarrowBandCache();

    virtual LevelSetNarrowBandCache * getNarrowBandCache();
    virtual const LevelSetNarrowBandCache * getNarrowBandCache() const;

    void clearNarrowBandCache();

    void dumpNarrowBandCache(std::ostream &stream);
    void restoreNarrowBandCache(std::istream &stream);

    bool isInNarrowBand(long id) const override;

    void swap(LevelSetCachedObjectInterface &other) noexcept;

protected:
    std::shared_ptr<LevelSetNarrowBandCache> m_narrowBandCache; //! Narrow band cache

    virtual std::shared_ptr<LevelSetNarrowBandCache> createNarrowBandCache() = 0;

};

class LevelSetCachedObject : public LevelSetObject, public LevelSetCachedObjectInterface, public LevelSetSignedObjectInterface {

    protected:

    void                                        _clear( ) override ;

    void                                        _clearAfterMeshAdaption(const std::vector<adaption::Info> & ) override ;

    void                                        _dump( std::ostream &) override ;
    void                                        _restore( std::istream &) override ;

# if BITPIT_ENABLE_MPI
    void                                        _writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) override ;
    void                                        _readCommunicationBuffer( const std::vector<long> &, RecvBuffer & )  override;
# endif 

    std::shared_ptr<LevelSetNarrowBandCache>    createNarrowBandCache() override;

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

#endif 
