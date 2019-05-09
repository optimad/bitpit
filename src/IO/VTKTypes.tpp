/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include <cassert>
#include <typeinfo>
#include <type_traits>

namespace bitpit{

/*!
 * Registers a type that will be used by the VTK functions
 * @tparam  T is the type that will be registered
 */
template<typename T>
VTKDataType VTKTypes::registerType()
{
    // If the type is already registered return the type
    // currently associated with it
    const std::type_info &typeInfo = typeid(T);
    VTKDataType VTKType = whichType(typeInfo);
    if (VTKType != VTKDataType::UNDEFINED) {
        return VTKType;
    }

    // Find out the type
    int size = 8 * sizeof(T);
    if (typeInfo == typeid(VTKElementType)) {
        if (size == 8) {
            VTKType = VTKDataType::Int8;
        } else if (size == 16) {
            VTKType = VTKDataType::Int16;
        } else if (size == 32) {
            VTKType = VTKDataType::Int32;
        } else if (size == 64) {
            VTKType = VTKDataType::Int64;
        }
    } else if (std::is_floating_point<T>::value) {
        if (size == 32) {
          VTKType = VTKDataType::Float32;
        } else if (size == 64) {
          VTKType = VTKDataType::Float64;
        }
    } else if (std::is_integral<T>::value) {
        if (!std::is_signed<T>::value) {
            if (size == 8) {
                VTKType = VTKDataType::UInt8;
            } else if (size == 16) {
                VTKType = VTKDataType::UInt16;
            } else if (size == 32) {
                VTKType = VTKDataType::UInt32;
            } else if (size == 64) {
                VTKType = VTKDataType::UInt64;
            }
        } else {
            if (size == 8) {
                VTKType = VTKDataType::Int8;
            } else if (size == 16) {
                VTKType = VTKDataType::Int16;
            } else if (size == 32) {
                VTKType = VTKDataType::Int32;
            } else if (size == 64) {
                VTKType = VTKDataType::Int64;
            }
        }
    }
    assert(VTKType != VTKDataType::UNDEFINED);

    return registerType<T>(VTKType);
}

/*!
 * Registers a type that will be used by the VTK functions
 * @tparam  T is the type that will be registered
 * @param[in]   VTKType is the VTK type that will be associated with the specified
 */
template<typename T>
VTKDataType VTKTypes::registerType(VTKDataType VTKType)
{
    // Index associated to the type
    const std::type_info &typeInfo = typeid(T);
    std::type_index index = std::type_index(typeInfo);

    // Add the type to the registered data types
    m_types.insert({{index, VTKType}});

    return VTKType;
}

/*!
 *  Determines the basic VTK type from argument.
 *  @tparam T type of argument
 *  @tparam nesting level of nesting
 *  @return basic VTK type 
 */
template<typename T, int nesting, typename std::enable_if<std::is_pod<T>::value&&!utils::is_iterable<T>::value>::type*>
VTKDataType VTKTypes::whichType( ){
    return whichType(typeid(T)) ;
}

/*!
 *  Determines the basic VTK type from argument.
 *  @tparam T type of argument
 *  @tparam nesting level of nesting
 *  @return basic VTK type 
 */
template<typename T, int nesting, typename std::enable_if<utils::is_iterable<T>::value>::type*>
VTKDataType VTKTypes::whichType( ){
    static_assert(nesting < 1,"Nested containers are not valid VTK data");
    return whichType<typename T::value_type,nesting+1>();
}

}
