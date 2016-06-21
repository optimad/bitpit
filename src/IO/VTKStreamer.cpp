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
 * @ingroup VisualizationToolKit
 * @interface VTKBaseContainer
 * @brief An interface class to all containers that are batively supported by VTK
 */

/*!
 * Destructor
 */
VTKBaseContainer::~VTKBaseContainer( ){

};

/*!
 * @ingroup VisualizationToolKit
 * @interface VTKBaseStreamer
 * @brief The base class to be used to derive VTK streamers form
 */

/*!
 * Reads data from stream 
 * @param[in] str file stream for writing
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKBaseStreamer::flushData( std::fstream &str, std::string name, VTKFormat format){

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
 * @param[in] entries total number of entries to be read
 * @param[in] components size of grouping (e.g. =3 for vectors)
 */
void VTKBaseStreamer::absorbData( std::fstream &str, std::string name, VTKFormat format, uint64_t entries, uint8_t components){

    BITPIT_UNUSED(str) ;
    BITPIT_UNUSED(name) ;
    BITPIT_UNUSED(format) ;
    BITPIT_UNUSED(entries) ;
    BITPIT_UNUSED(components) ;

    return;
};

/*!
 * @ingroup VisualizationToolKit
 * @class VTKNativeStreamer
 * @brief In VTKNativeStreamer all instances of classes derived from VTKBaseConatiner are stored. 
 * Right now only std::vector is supported.
 */

/*!
 * Destructor
 */
VTKNativeStreamer::~VTKNativeStreamer(){

    m_field.clear() ;
};

/*!
 * Constructor
 */
VTKNativeStreamer::VTKNativeStreamer( ) {

};

/*!
 * Removes a field from streamer
 * @param[in] name name of field
 */
void VTKNativeStreamer::removeData( std::string name){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr != m_field.end()){
        m_field.erase( fieldItr) ;
    }

    return;
};

/*!
 * Writes data to stream 
 * @param[in] str file stream for writing
 * @param[in] name name of field
 * @param[in] format ASCII or BINARY format
 */
void VTKNativeStreamer::flushData( std::fstream &str, std::string name, VTKFormat format){

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
 * @param[in] entries total number of entries to be read
 * @param[in] components size of groups
 */
void VTKNativeStreamer::absorbData( std::fstream &str, std::string name, VTKFormat format, uint64_t entries, uint8_t components){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr != m_field.end()){
        fieldItr->second->absorbData( str, format, entries, components ) ;
    }

    return;
};

}
