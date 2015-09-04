template<class T>
string VTK::WhichType( T dummy ){

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

//-----------------------------------------------------------
//-----------------------------------------------------------
template<class T>
string VTK::WhichType( vector<T> dummy ){

    T dummy2 ;
    return  WhichType(dummy2) ;
};

//-----------------------------------------------------------
//-----------------------------------------------------------
template<class T, size_t d>
string VTK::WhichType( std::array<T,d> dummy ){

    T dummy2 ;
    return  WhichType(dummy2) ;
};
