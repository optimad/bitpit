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

# ifndef __BITPIT_LEVELSET_SIGN_PROPAGATOR_HPP__
# define __BITPIT_LEVELSET_SIGN_PROPAGATOR_HPP__

# include "bitpit_containers.hpp"
# include "bitpit_patchkernel.hpp"

namespace bitpit{

class LevelSetSignPropagator;

class LevelSetSignStorage {

friend class LevelSetSignPropagator;

public:
    typedef signed char Sign;

    static const Sign SIGN_UNDEFINED;
    static const Sign SIGN_NEGATIVE;
    static const Sign SIGN_ZERO;
    static const Sign SIGN_POSITIVE;

    bool isStoredSignDirty() const;

    Sign getStoredSign(long id) const;
    Sign getStoredSign(const VolumeKernel::CellConstIterator &itr) const;

    Sign rawGetStoredSign(std::size_t rawIndex) const;

protected:
    LevelSetSignStorage();

    bool isSignStorageInitialized() const;
    void initializeSignStorage(PiercedKernel<long> *cellKernel);
    void clearSignStorage(bool release = true);

    void setStoredSignDirty(bool dirty);

    void setStoredSign(Sign sign);
    void setStoredSign(const VolumeKernel::CellConstIterator &itr, Sign sign);

    void dumpStoredSign(std::ostream &stream);
    void restoreStoredSign(std::istream &stream);

private:
    bool m_dirty; /** Check if the storage is dirty */
    PiercedStorage<signed char> m_storage; /** Storage for the levelset sign */

};

class LevelSetSignPropagator {

public:
    LevelSetSignPropagator(VolumeKernel *mesh);

    void execute(const LevelSetObject *object, LevelSetSignStorage *storage);
    void execute(const std::vector<adaption::Info> &adaptionData, const LevelSetObject *object, LevelSetSignStorage *storage);

private:
    typedef signed char PropagationState;

    static const PropagationState STATE_EXTERNAL;
    static const PropagationState STATE_WAITING;
    static const PropagationState STATE_REACHED;

    VolumeKernel *m_mesh;

    long m_nWaiting;
    long m_nExternal;
    LevelSetSignStorage::Sign m_externalSign;
    PiercedStorage<PropagationState, long> m_propagationStates;

    void propagate(const LevelSetObject *object, LevelSetSignStorage *storage);

    void initializePropagation(const LevelSetObject *object);
    void executeSeedPropagation(const std::vector<std::size_t> &rawSeeds, LevelSetSignStorage *storage);
    void finalizePropagation();

    void setSign(const VolumeKernel::CellConstIterator &cellItr, LevelSetSignStorage::Sign sign, LevelSetSignStorage *storage);

};

}

#endif 
