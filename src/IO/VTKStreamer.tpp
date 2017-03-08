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

namespace bitpit{

/*!
 * @ingroup VisualizationToolKit
 * @class VTKVectorContainer
 * @brief Implementation of VTKBaseContainer in order to support natively std::vector in VTK
 */

/*!
 * Destructor
 */
template<class T>
VTKVectorContainer<T>::~VTKVectorContainer( ){
    m_ptr = NULL;
}

/*!
 * Copy constructor
 */
template<class T>
VTKVectorContainer<T>::VTKVectorContainer( const VTKVectorContainer &other)
    : VTKBaseContainer(other) {
    m_ptr = other.m_ptr;
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
template<class T>
VTKVectorContainer<T> * VTKVectorContainer<T>::clone() const {
    return new VTKVectorContainer<T>( *this );
}

/*!
 * Constructor assigns container reference to internal pointer
 * @tparam T type of std::vector<>
 * @param[in] data data container
 */
template<class T>
VTKVectorContainer<T>::VTKVectorContainer( std::vector<T> &data){
    m_ptr = &data;
}

/*!
 * Writes data to file stream
 * @tparam T type of std::vector<>
 * @param[in] str file stream
 * @param[in] format VTKFormat::ASCII or VTKFormat::APPENDED
 */
template<class T>
void VTKVectorContainer<T>::flushData( std::fstream &str, VTKFormat format){

    if( format==VTKFormat::ASCII){
        genericIO::flushASCII( str, 8, *m_ptr) ;

    } else if( format==VTKFormat::APPENDED) {
        genericIO::flushBINARY( str, *m_ptr) ;

    }

}

/*!
 * Reads data from file stream
 * @tparam T type of std::vector<>
 * @param[in] str file stream
 * @param[in] format VTKFormat::ASCII or VTKFormat::APPENDED
 * @param[in] entries total number of entries to be read
 * @param[in] components size of inner chunks
 */
template<class T>
void VTKVectorContainer<T>::absorbData( std::fstream &str, VTKFormat format, uint64_t entries, uint8_t components){

    resize( std::is_fundamental<T>{}, entries, components) ;

    if( format==VTKFormat::ASCII){
        genericIO::absorbASCII( str, *m_ptr) ;

    } else if( format==VTKFormat::APPENDED) {
        genericIO::absorbBINARY( str, *m_ptr) ;

    }

}

/*!
 * Overloaded version of resize for std::vector<fundamental_type>
 * @tparam T type of std::vector<>
 * @param[in] entries total number of entries to be read
 * @param[in] components size of inner chunks
 */
template<class T>
void VTKVectorContainer<T>::resize( std::true_type, uint64_t entries, uint8_t components){

    BITPIT_UNUSED(components) ;
    m_ptr->resize(entries) ;

}

/*!
 * Overloaded version of resize for std::vector<non_fundamental_type>
 * @tparam T type of std::vector<>
 * @param[in] entries total number of entries to be read
 * @param[in] components size of inner chunks
 */
template<class T >
void VTKVectorContainer<T>::resize( std::false_type, uint64_t entries, uint8_t components){

    uint64_t elements = entries /components ;
    m_ptr->resize(elements) ;

    for( auto & element : (*m_ptr) ){
        vtk::allocate( element, components) ;
    }

}

/*!
 * Adds data strored in std::vector<> to NativeStreamer
 * @tparam T type of std::vector<>
 * @param[in] name name of data set
 * @param[in] data std::vector containing the data
 */
template<class T>
void VTKNativeStreamer::addData( std::string name, std::vector<T> &data ){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr!=m_field.end() ){
        m_field.erase(fieldItr) ;
    }
    
    std::unique_ptr<VTKBaseContainer> temp = std::unique_ptr<VTKBaseContainer>( new VTKVectorContainer<T>(data) ) ;
    m_field.emplace( name, std::move(temp) ) ;

}

}
