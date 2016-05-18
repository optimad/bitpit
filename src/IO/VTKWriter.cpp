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

#include "VTK.hpp"

namespace bitpit{


/*!
 * @ingroup    VisualizationToolKit
 * @{
 */

/*!
 * @class VTKBaseContainer
 * @brief An interface class to all containers that are batively supported by VTK
 */

/*!
 * @class VTKBaseWriter
 * @brief The base class to be used to derive VTK Writers form
 */

/*!
 * Reads data from stream 
 * @param[in] str file stream for writing
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKBaseWriter::flushData( std::fstream &str, std::string name, VTKFormat format){

    BITPIT_UNUSED(str) ;
    BITPIT_UNUSED(name) ;
    BITPIT_UNUSED(format) ;

    return;
};

/*!
 * Reads data from stream 
 * @param[in] str file stream for writing
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKBaseWriter::absorbData( std::fstream &str, std::string name, VTKFormat format, uint64_t entries, uint8_t components){

    BITPIT_UNUSED(str) ;
    BITPIT_UNUSED(name) ;
    BITPIT_UNUSED(format) ;
    BITPIT_UNUSED(entries) ;
    BITPIT_UNUSED(components) ;

    return;
};

/*!
 * @class VTKNativeWriter
 * @brief A VTK writer which supports natively std::vector
 */

/*!
 * Destructor
 */
VTKNativeWriter::~VTKNativeWriter(){

    for( auto &field : m_field){
        delete field.second ;
    };

    m_field.clear() ;
    owner = NULL ;

};

/*!
 * Constructor
 */
VTKNativeWriter::VTKNativeWriter(VTK& owner_) :owner(&owner_){

};

/*!
 * Removes a field from writer
 * @param[in] name name of field
 */
void VTKNativeWriter::removeData( std::string name){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr != m_field.end()){
        m_field.erase( fieldItr) ;
    }

    owner->removeData(name) ;

    return;
};

/*!
 * Writes data to stream 
 * @param[in] str file stream for writing
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKNativeWriter::flushData( std::fstream &str, std::string name, VTKFormat format){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr != m_field.end()){
        fieldItr->second->flushData(str,format) ;
    }

    return;
};

/*!
 * Reads data from stream 
 * @param[in] str file stream for reading
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKNativeWriter::absorbData( std::fstream &str, std::string name, VTKFormat format, uint64_t entries, uint8_t components){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr != m_field.end()){
        fieldItr->second->absorbData( str, format, entries, components ) ;
    }

    return;
};
/*!
 * @}
 */

}
