/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include "communications_tags.hpp"

namespace bitpit {

/*!
	\class CommunicationTags
	\ingroup communications

	\brief The CommunicationTags oversee the handling of the communication
	tags.
*/

/*
    Initialize communication tags instance.
*/
std::unique_ptr<CommunicationTags> CommunicationTags::m_instance = nullptr;

/*!
    Returns an instance of the communication tag manager.

    \result An instance of the communication tag manager.
*/
CommunicationTags & CommunicationTags::instance()
{
    if (!m_instance) {
        m_instance = std::unique_ptr<CommunicationTags>(new CommunicationTags());
    }

    return *m_instance;
}

/*!
	Creates a new communication tag manager.
*/
CommunicationTags::CommunicationTags()
{
}

// Communication manager global functions
namespace communications {

    /*!
        Returns the communication tag manager.

        \result The communication tag manager.
    */
    CommunicationTags & tags()
    {
        return CommunicationTags::instance();
    }

}

}
