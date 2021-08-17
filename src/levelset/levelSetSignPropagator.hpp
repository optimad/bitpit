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

# include "bitpit_patchkernel.hpp"

# include "levelSetSignedObject.hpp"

namespace bitpit{

class LevelSetSignPropagator {

public:
    LevelSetSignPropagator(VolumeKernel *mesh);

    void execute(LevelSetSignedObjectInterface *object);
    void execute(const std::vector<adaption::Info> &adaptionData, LevelSetSignedObjectInterface *object);

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

    void propagate(const LevelSetObjectInterface *object, LevelSetSignStorage *storage);

    void initializePropagation(const LevelSetObjectInterface *object);
    void executeSeedPropagation(const std::vector<std::size_t> &rawSeeds, LevelSetSignStorage *storage);
    void finalizePropagation(LevelSetSignStorage *storage);

    void setSign(std::size_t cellRawId, LevelSetSignStorage::Sign sign, LevelSetSignStorage *storage);

};

}

#endif 
