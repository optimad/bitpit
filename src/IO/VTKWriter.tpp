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

namespace bitpit{

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

    return ;

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

    return ;

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

    return ;
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

    return ;
}

/*!
 * Adds data strored in std::vector<> to NativeWriter
 * @tparam T type of std::vector<>
 * @param[in] name name of data set
 * @param[in] data std::vector containing the data
 */
template<class T>
void VTKNativeWriter::addData( std::string name, std::vector<T> &data){

    auto fieldItr = m_field.find(name) ;

    if( fieldItr == m_field.end() ){
        m_field.insert( {{ name, new VTKVectorContainer<T>(data) }}) ;
        owner->setWriter( name, this) ;

    } else {
        log::cout() << "Field with name " << name << " already exists in VTKNativeWriter " << std::endl ;
    }



    return ;

}


}
