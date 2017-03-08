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

#include "VTK.hpp"

namespace bitpit{


// Initialize static variables
std::unordered_map<std::type_index, VTKDataType> VTKTypes::m_types = std::unordered_map<std::type_index, VTKDataType>();

// Register standadd types
static VTKDataType VTKTypeIndex_sint   = VTKTypes::registerType<short int>();
static VTKDataType VTKTypeIndex_int    = VTKTypes::registerType<int>();
static VTKDataType VTKTypeIndex_long   = VTKTypes::registerType<long>();
static VTKDataType VTKTypeIndex_llong  = VTKTypes::registerType<long long>();
static VTKDataType VTKTypeIndex_usint  = VTKTypes::registerType<unsigned short int>();
static VTKDataType VTKTypeIndex_uint   = VTKTypes::registerType<unsigned int>();
static VTKDataType VTKTypeIndex_ulong  = VTKTypes::registerType<unsigned long>();
static VTKDataType VTKTypeIndex_ullong = VTKTypes::registerType<unsigned long long>();
static VTKDataType VTKTypeIndex_uint8  = VTKTypes::registerType<uint8_t>();
static VTKDataType VTKTypeIndex_uint16 = VTKTypes::registerType<uint16_t>();
static VTKDataType VTKTypeIndex_uint32 = VTKTypes::registerType<uint32_t>();
static VTKDataType VTKTypeIndex_uint64 = VTKTypes::registerType<uint64_t>();
static VTKDataType VTKTypeIndex_int8   = VTKTypes::registerType<int8_t>();
static VTKDataType VTKTypeIndex_int16  = VTKTypes::registerType<int16_t>();
static VTKDataType VTKTypeIndex_int32  = VTKTypes::registerType<int32_t>();
static VTKDataType VTKTypeIndex_int64  = VTKTypes::registerType<int64_t>();
static VTKDataType VTKTypeIndex_float  = VTKTypes::registerType<float>();
static VTKDataType VTKTypeIndex_double = VTKTypes::registerType<double>();

// Register VTK types
static VTKDataType VTKTypeIndex_VTKElement = VTKTypes::registerType<VTKElementType>();

/*!
 * @ingroup    VisualizationToolKit
 * @{
 *
 * @class VTKTypes
 * @brief VTK data types handling
 */

/*!
 * Calculates the size in bytes of basic types supported by VTK
 * @param[in]  type    type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * @return      size of basic type
 */
uint8_t VTKTypes::sizeOfType( const VTKDataType & type ){

    switch( type){
        case VTKDataType::Int8 : case VTKDataType::UInt8:
           return( sizeof(int8_t) );

        case VTKDataType::Int16 : case VTKDataType::UInt16 :
           return( sizeof(int16_t) ) ;

        case VTKDataType::Int32 : case VTKDataType::UInt32 :
           return( sizeof(int32_t) ) ;

        case VTKDataType::Int64 : case VTKDataType::UInt64 :
           return( sizeof(int64_t) ) ;

        case VTKDataType::Float32 :
           return( sizeof(float) ) ;

        case VTKDataType::Float64 :
           return( sizeof(double) ) ;

        case VTKDataType::UNDEFINED :
           return( 0 ) ;

        default:
           return(0);
    }

}

/*!
 *  Determines the basic VTK type from argument.
 *  @param[in]  typeInfo argument type
 *  @return     basic VTK type
 */
VTKDataType VTKTypes::whichType( const std::type_info & typeInfo ){

    std::type_index index = std::type_index(typeInfo);
    if (m_types.count(index) == 0) {
        return VTKDataType::UNDEFINED;
    }

    return m_types.at(typeInfo);

}

/*!
 * @}
 */

}
