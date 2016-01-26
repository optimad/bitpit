/*!
 * @ingroup    VTKUtils
 * @{
 */


/*!
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy       argument type to be deduced
 *  @return     basic VTK type 
 *
 */
template<class T>
VTKDataType VTKUtils::whichType( T dummy ){

    BITPIT_UNUSED(dummy) ;

    const std::type_info &type = typeid(dummy) ;
    return VTKUtils::whichType( type) ;

//    if(std::is_same<T,int8_t>::value){
//        return VTKDataType::Int8 ;
//    }
//
//    else if(std::is_same<T,uint8_t>::value){
//        return VTKDataType::UInt8 ;
//    }
//
//    else if(std::is_same<T,int16_t>::value){
//        return VTKDataType::Int16 ;
//    }
//
//    else if(std::is_same<T,uint16_t>::value){
//        return VTKDataType::UInt16 ;
//    }
//
//    else if(std::is_same<T,int32_t>::value){
//        return VTKDataType::Int32 ;
//    }
//
//    else if(std::is_same<T,uint32_t>::value){
//        return VTKDataType::UInt32 ;
//    }
//
//    else if(std::is_same<T,int64_t>::value){
//        return VTKDataType::Int64 ;
//    }
//
//    else if(std::is_same<T,uint64_t>::value){
//        return VTKDataType::UInt64 ;
//    }
//
//    else if(std::is_same<T,float>::value){
//        return VTKDataType::Float32 ;
//    }
//
//    else if(std::is_same<T,double>::value){
//        return VTKDataType::Float64 ;
//    }
//
//
//    else{
//        return VTKDataType::UNDEFINED ;
//    };


};

/*!
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy      argument type to be deduced
 *  @return     basic VTK type 
 */
template<class T>
VTKDataType VTKUtils::whichType( std::vector<T> dummy ){

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
VTKDataType VTKUtils::whichType( std::array<T,d> dummy ){

    BITPIT_UNUSED(dummy) ;

    T dummy2 ;
    return  whichType(dummy2) ;
};

/*!
 * @}
 */
