#include <assert.h>
#include <typeinfo>
#include <type_traits>

/*!
 * Registers a type that will be used by the VTK functions
 * @tparam  T is the type that will be registered
 * @param   VTKType is the VTK type that will be associated with the specified
 * type
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
 *  @tparam     T           type of argument
 *  @param[in]  dummy       argument type to be deduced
 *  @return     basic VTK type 
 *
 */
template<class T>
VTKDataType VTKTypes::whichType( T dummy ){

    BITPIT_UNUSED(dummy) ;

    const std::type_info &typeInfo = typeid(dummy) ;
    return whichType(typeInfo) ;
};

/*!
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy      argument type to be deduced
 *  @return     basic VTK type 
 */
template<class T>
VTKDataType VTKTypes::whichType( std::vector<T> dummy ){

    BITPIT_UNUSED(dummy) ;

    T dummy2 ;
    return  whichType(dummy2) ;
};

/*!
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy      argument type to be deduced
 *  @return     basic VTK type 
 */
template<class T, size_t d>
VTKDataType VTKTypes::whichType( std::array<T,d> dummy ){

    BITPIT_UNUSED(dummy) ;

    T dummy2 ;
    return  whichType(dummy2) ;
};

