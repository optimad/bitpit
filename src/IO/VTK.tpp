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
 * Adds data strored in std::vector<> to NativeStreamer
 * @tparam T type of std::vector<>
 * @param[in] name name of data set
 * @param[in] data std::vector containing the data
 */
template<class T>
VTKField& VTK::addData( std::string name, std::vector<T> &data ){

    VTKField& field= addData( name, &nativeStreamer) ;
    field.setDataType( VTKTypes::whichType(data)) ;

    nativeStreamer.addData(name,data) ;

    return field;

}

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @tparam T type of std::vector<>
 * @param[in]  name_    name of field
 * @param[in]  comp_    type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 * @param[in]  loc_     location of data [VTKLocation::CELL/VTKLocation::POINT]
 */
template<class T>
VTKField& VTK::addData( std::string name_, VTKFieldType comp_,  VTKLocation loc_, std::vector<T> &data ){

    VTKField&   field = addData( name_, data ) ;

    field.setFieldType(comp_) ;
    field.setLocation(loc_) ;

    return field ;

};

/*!
 * Associates the NativeStreamer to a geometrical field
 * @tparam T type of std::vector<>
 * @param[in] fieldEnum which geometrical field 
 * @param[in] data std::vector containing the data
 */
template<class T>
void VTKUnstructuredGrid::setGeomData( VTKUnstructuredField fieldEnum, std::vector<T> &data ){

    int      index = static_cast<int>(fieldEnum) ;
    VTKField& field = geometry[index] ;
    std::string   name = field.getName() ;

    nativeStreamer.addData(name,data) ;

    field.setStreamer( nativeStreamer) ;
    field.setDataType( VTKTypes::whichType(data)) ;

    return ;
}

/*!
 * Associates the NativeStreamer to a geometrical field
 * @tparam T type of std::vector<>
 * @param[in] fieldEnum which geometrical field 
 * @param[in] data std::vector containing the data
 */
template<class T>
void VTKRectilinearGrid::setGeomData( VTKRectilinearField fieldEnum, std::vector<T> &data ){

    int      index = static_cast<int>(fieldEnum) ;
    VTKField& field = geometry[index] ;
    std::string   name = field.getName() ;

    nativeStreamer.addData(name,data) ;

    field.setStreamer( nativeStreamer) ;
    field.setDataType( VTKTypes::whichType(data)) ;

    return ;
}
}
