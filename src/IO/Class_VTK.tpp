/*! -------------------------------------------------------
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy       argument type to be deduced
 *  @return     basic VTK type [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64, unsupported ]
 *
 */
template<class T>
std::string VTK::WhichType( T dummy ){

    if(std::is_same<T,int8_t>::value){
        return "Int8" ;
    }

    else if(std::is_same<T,uint8_t>::value){
        return "UInt8" ;
    }

    else if(std::is_same<T,int16_t>::value){
        return "Int16" ;
    }

    else if(std::is_same<T,uint16_t>::value){
        return "UInt16" ;
    }

    else if(std::is_same<T,int32_t>::value){
        return "Int32" ;
    }

    else if(std::is_same<T,uint32_t>::value){
        return "UInt32" ;
    }

    else if(std::is_same<T,int64_t>::value){
        return "Int64" ;
    }

    else if(std::is_same<T,uint64_t>::value){
        return "UInt64" ;
    }

    else if(std::is_same<T,float>::value){
        return "Float32" ;
    }

    else if(std::is_same<T,double>::value){
        return "Float64" ;
    }

    else if(std::is_same<T,double>::value){
        return "Float64" ;
    }


    else{
        return "unsupported" ;
    };


};

/*! -------------------------------------------------------
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy      argument type to be deduced
 *  @return     basic VTK type [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64, unsupported ]
 */
template<class T>
std::string VTK::WhichType( std::vector<T> dummy ){

    T dummy2 ;
    return  WhichType(dummy2) ;
};

/*! -------------------------------------------------------
 *  Determines the basic VTK type from argument.
 *  @tparam     T           type of argument
 *  @param[in]  dummy      argument type to be deduced
 *  @return     basic VTK type [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64, unsupported ]
 */
template<class T, size_t d>
std::string VTK::WhichType( std::array<T,d> dummy ){

    T dummy2 ;
    return  WhichType(dummy2) ;
};
