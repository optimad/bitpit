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
 *  as published by by the Free Software Foundation.
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

#include "VTK.hpp"

namespace bitpit{

/*!
 * @ingroup     VisualizationToolKit
 * @{
 *
 * @class       VTKFieldMetaData
 * @brief       A class used in the interface VTK::getMetaData( std::string )
 *
 */

/*! 
 * Constructor with all members set
 * @param[in] size the entire size of field date, e.g. in case of a vector field on nodes size = NNodes x 3 
 * @param[in] type the type of the basic data used
 */
VTKFieldMetaData::VTKFieldMetaData( uint64_t size, const std::type_info &type): m_size(size), m_type(type){
};

/*! 
 * Get the size of field
 * @return size of field
 */
uint64_t VTKFieldMetaData::getSize() const{
    return m_size;
}

/*! 
 * Get the type of data used in field
 * @return data type of field
 */
const std::type_info& VTKFieldMetaData::getType() const{
    return m_type;
}

/*! 
 * @} 
 */

}
