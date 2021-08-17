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

# ifndef __BITPIT_LEVELSET_SIGNED_OBJECT_HPP__
# define __BITPIT_LEVELSET_SIGNED_OBJECT_HPP__

# include "bitpit_containers.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetObject.hpp"
# include "levelSetStorage.hpp"

namespace bitpit{

class LevelSetSignStorage : public LevelSetExternalPiercedStorageManager
{

public:
    typedef signed char Sign;

    static const Sign SIGN_UNDEFINED;
    static const Sign SIGN_NEGATIVE;
    static const Sign SIGN_ZERO;
    static const Sign SIGN_POSITIVE;

    LevelSetSignStorage(PiercedKernel<long> *kernel);

    Sign at(const KernelIterator &itr) const;
    Sign & at(const KernelIterator &itr);

    void fill(Sign sign);

    void swap(LevelSetSignStorage &other) noexcept;

protected:
    Storage<Sign> *m_signs; //! Levelset signs of the cells

};

class LevelSetSignedObjectInterface : public virtual LevelSetObjectInterface {

public:
    LevelSetSignStorage * initializeSignStorage();

    LevelSetSignStorage * getSignStorage();
    const LevelSetSignStorage * getSignStorage() const;

    bool isSignStorageDirty() const;
    void setSignStorageDirty(bool available);

    void clearSignStorage();

    void dumpSignStorage(std::ostream &stream);
    void restoreSignStorage(std::istream &stream);

    void swap(LevelSetSignedObjectInterface &other) noexcept;

protected:
    std::shared_ptr<LevelSetSignStorage> m_signStorage; //! Storage for levelset signs.

    virtual std::shared_ptr<LevelSetSignStorage> createSignStorage() = 0;

};

}

#endif
