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

#include "proxyVector.hpp"

namespace bitpit{

// Some commonly used ProxyVectors are instantiated explicitly
template class ProxyVector<int, true>;
template class ProxyVector<long, true>;
template class ProxyVector<double, true>;

template class ProxyVector<int, false>;
template class ProxyVector<long, false>;
template class ProxyVector<double, false>;

template class ProxyVector<const int, true>;
template class ProxyVector<const long, true>;
template class ProxyVector<const double, true>;

template class ProxyVector<const int, false>;
template class ProxyVector<const long, false>;
template class ProxyVector<const double, false>;

}
