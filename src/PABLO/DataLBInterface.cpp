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

#include <DataLBInterface.hpp>

namespace bitpit {

std::size_t DummyDataLBImpl::size(const uint32_t e) const
{
    BITPIT_UNUSED(e);

    return 0;
}

std::size_t DummyDataLBImpl::fixedSize() const
{
    return 0;
}

void DummyDataLBImpl::move(const uint32_t from, const uint32_t to)
{
    BITPIT_UNUSED(from);
    BITPIT_UNUSED(to);
}

void DummyDataLBImpl::assign(uint32_t stride, uint32_t length)
{
    BITPIT_UNUSED(stride);
    BITPIT_UNUSED(length);
}

void DummyDataLBImpl::resize(uint32_t newSize)
{
    BITPIT_UNUSED(newSize);
}

void DummyDataLBImpl::resizeGhost(uint32_t newSize)
{
    BITPIT_UNUSED(newSize);
}

void DummyDataLBImpl::shrink()
{
}

};
