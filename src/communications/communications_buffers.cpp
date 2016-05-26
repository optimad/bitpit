/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include "communications_buffers.hpp"

namespace bitpit {

/*!
    \ingroup communications
    \class SendBuffer

    @brief Buffer to be used for send communications.
*/

/*!
    Creates a new send buffer
*/
SendBuffer::SendBuffer(size_t capacity, bool doubleBuffer)
    : CommunicationBuffer(capacity, doubleBuffer)
{
}

/*!
    \ingroup communications
    \class RecvBuffer

    @brief Buffer to be used for receive communications.
*/

/*!
    Creates a new receive buffer
*/
RecvBuffer::RecvBuffer(size_t capacity, bool doubleBuffer)
    : CommunicationBuffer(capacity, doubleBuffer)
{
}

}
