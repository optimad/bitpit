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
# include "levelSetSignPropagator.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class SendBuffer;
class RecvBuffer;

class LevelSetObject;

class LevelSetCachedObject : public LevelSetObject, public LevelSetBoundedObject, public LevelSetSignStorage {

    protected:

    typedef PiercedKernel<long>                 NarrowBandKernel;
    typedef NarrowBandKernel::const_iterator    NarrowBandIterator;

    template<typename T>
    using NarrowBandStorage = PiercedStorage<T, long>;

    NarrowBandKernel                            m_narrowBandKernel;    /** Pierced kernel for storing cell information in the narrow band */
    NarrowBandStorage<double>                   m_narrowBandValues;    /** Levelset values of the cells inside the narrow band */
    NarrowBandStorage<std::array<double, 3>>    m_narrowBandGradients; /** Levelset gradient of the cells inside the narrow band */

    void                                        _clear( ) override ;

    void                                        _clearAfterMeshAdaption(const std::vector<adaption::Info> & ) override ;

    void                                        _dump( std::ostream &) override ;
    void                                        _restore( std::istream &) override ;

# if BITPIT_ENABLE_MPI
    void                                        _writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) override ;
    void                                        _readCommunicationBuffer( const std::vector<long> &, RecvBuffer & )  override;
# endif 

    NarrowBandIterator                          createNarrowBandEntry(long id, bool sync = true);
    void                                        eraseNarrowBandEntry(long id, bool sync = true);
    bool                                        containsNarrowBandEntry(long id) const;
    NarrowBandIterator                          getNarrowBandEntryIterator(long id) const;
    void                                        setNarrowBandEntry(NarrowBandIterator itr, double distance, const std::array<double, 3> &gradient);

    void                                        syncNarrowBandStorages();

    # if BITPIT_ENABLE_MPI
    virtual std::size_t                         getNarrowBandEntryBinarySize() const;
    virtual void                                writeNarrowBandEntryCommunicationBuffer( NarrowBandIterator narrowBandItr, SendBuffer &dataBuffer ) ;
    virtual void                                readNarrowBandEntryCommunicationBuffer( NarrowBandIterator narrowBandItr, RecvBuffer &dataBuffer ) ;
#endif

    public:
    LevelSetCachedObject(int);

    LevelSetInfo                                getLevelSetInfo(long ) const override ;
    double                                      getLS(long ) const override ;
    double                                      getValue(long ) const override ;
    short                                       getSign(long ) const override ;
    std::array<double,3>                        getGradient(long ) const override ;

    bool                                        isInNarrowBand(long id) const override;

};



}

#endif 
