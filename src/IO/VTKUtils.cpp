#include "VTK.hpp"

/*!
 * @ingroup    VTKUtils
 * @{
 *
 * @brief Utility fuctions for VTK
 */

/*!
 * Calculates the size in bytes of basic types supported by VTK
 * @param[in]  type    type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * \return      size of basic type
 */
uint8_t VTKUtils::sizeOfType( const VTKDataType & type ){

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
    };

};

/*!  
 *  Returns number of vertices for a given element type.
 *  Codification of element type according to http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *  @param[in]  t           element type
 *  \return     number of vertices of element type
 */
uint8_t VTKUtils::getNNodeInElement( const VTKElementType & t){
    int myType = static_cast<std::underlying_type<VTKElementType>::type>(t);

    switch (myType){

        case 1: 
            return(1);
        case 3: 
            return(2);
        case 5: case 21: 
            return(3);
        case 8: case 9:  case 10: 
            return(4);
        case 11: case 12: case 23: case 24: 
            return(8);
        case 13: case 22: 
            return(6);
        case 14: 
            return(5) ;
        case 25: 
            return(20);
        default:
           return(-1);
    };


};

/*!
 * Converts a string into a Field information.
 * @param[in]   line_       string to be converted
 * @param[out]  data_       Field information
 * \return      true if successful
 */
bool VTKUtils::convertStringToDataArray( const std::string &line_, VTKField &data_  ){

    std::string type_, name_, code_, comp_, offs_ ;
    int         components_(1), offset_ ;

    VTKFormat    codex ;
    VTKDataType  type ;
    VTKFieldType comp(VTKFieldType::SCALAR) ;

    bool  success(true) ;


    if( Keyword_In_String( line_, "<DataArray ") ){  
        success = success && Get_After_Keyword( line_, "type=", '\"', type_) ;
        success = success && Get_After_Keyword( line_, "Name=", '\"', name_) ;
        success = success && Get_After_Keyword( line_, "format=", '\"', code_) ;

        if( Get_After_Keyword( line_, "NumberOfComponents=", '\"', comp_)  ){
            convert_string( comp_, components_ ) ;
        };

        if(components_==3)
            comp=VTKFieldType::VECTOR ;

        convertStringToEnum( type_, type) ;
        convertStringToEnum( code_, codex) ;

        data_.setType(type) ;
        data_.setName(name_) ;
        data_.setComponents(comp) ;
        data_.setCodification(codex) ;

        if(code_=="appended") {
            if( Get_After_Keyword( line_, "offset=", '\"', offs_) ){
                convert_string( offs_, offset_ ) ;
                data_.setOffset(offset_) ;
            }
            else{
                success = false ;
            };
        }

        return success ;
    }

    else{
        return false ;
    };


};

/*!
 * Converts a Field information to string as requested by VTK format.
 * @param[in]  field_       Field information
 * @return string in VTK format
 */
std::string  VTKUtils::convertDataArrayToString( const VTKField &field_ ){

    std::stringstream os("") ;

    os << "        <DataArray "
        << "type=\"" << VTKUtils::convertEnumToString( field_.getType() ) << "\" "
        << "Name=\"" << field_.getName() << "\" "
        << "NumberOfComponents=\""<< unsigned(field_.getComponents()) << "\" "
        << "format=\"" << VTKUtils::convertEnumToString(field_.getCodification()) << "\" ";

    if( field_.getCodification() == VTKFormat::APPENDED ){
        os << "offset=\"" << field_.getOffset() << "\" " ;
    };

    os << ">" ;

    return( os.str() ) ;       


};

/*!
 * Converts a parallel field information to string as requested by VTK format.
 * @param[in]  field_       Field information
 * @return string in VTK format
 */
std::string  VTKUtils::convertPDataArrayToString( const VTKField &field_ ){

    std::stringstream  os("") ;

    os << "        <PDataArray "
        << "type=\"" << VTKUtils::convertEnumToString(field_.getType()) << "\" "
        << "Name=\"" << field_.getName() << "\" "
        << "NumberOfComponents=\""<< unsigned(field_.getComponents()) << "\" " 
        << ">" ;

    return( os.str() ) ;

};

/*!
 * Converts a VTKLocation into string for DataArray format
 * @param[in]  loc VTKLocation to be converted
 * @return string to be used in DataArray
 */
std::string VTKUtils::convertEnumToString( const VTKLocation &loc ){
   
    switch(loc){
        case VTKLocation::CELL :
            return("Cell");
        case VTKLocation::POINT :
            return("Point");
        case VTKLocation::UNDEFINED :
            return("Undefined") ;
        default:
            return("Undefined") ;
    };
};

/*!
 * Converts a VTKFormat into string for DataArray format
 * @param[in]  loc VTKFormat to be converted
 * @return string to be used in DataArray
 */
std::string VTKUtils::convertEnumToString( const VTKFormat &cod ){
    
    switch(cod){
        case VTKFormat::ASCII :
            return("ascii");
        case VTKFormat::APPENDED :
            return("appended");
        case VTKFormat::UNDEFINED :
            return("Undefined") ;
        default:
            return("Undefined") ;
    };
};

/*!
 * Converts a VTKDataType into string for DataArray format
 * @param[in]  type VTKDataType to be converted
 * @return string to be used in DataArray
 */
std::string VTKUtils::convertEnumToString( const VTKDataType &type ){
    
    switch(type){
        case VTKDataType::Int8 :
            return("Int8");
        case VTKDataType::Int16 :
            return("Int16");
        case VTKDataType::Int32 :
            return("Int32");
        case VTKDataType::Int64 :
            return("Int64");
        case VTKDataType::UInt8 :
            return("UInt8");
        case VTKDataType::UInt16 :
            return("UInt16");
        case VTKDataType::UInt32 :
            return("UInt32");
        case VTKDataType::UInt64 :
            return("UInt64");
        case VTKDataType::Float32 :
            return("Float32");
        case VTKDataType::Float64 :
            return("Float64");
        case VTKDataType::UNDEFINED :
            return("Undefined");
        default:
            return("Undefined") ;
    };
};

/*!
 * Converts a std::string as read in vtk file to VTKLocation
 * @param[in] str string read in DataArray header
 * @param[out]  loc VTKLocation
 * @return  if str contained expected value
 */
bool VTKUtils::convertStringToEnum(  const std::string &str, VTKLocation &loc ){
    

    if( str == "Cell"){
        loc = VTKLocation::CELL;
        return(true);

    } else if ( str == "Point" ){
        loc = VTKLocation::POINT;
        return(true);

    } else {
        loc = VTKLocation::UNDEFINED ;
        return(false);
    };

};

/*!
 * Converts a std::string as read in vtk file to VTKFormat  
 * @param[in] str string read in DataArray header
 * @param[out]  cod VTKFormat
 * @return  if str contained expected value
 */
bool VTKUtils::convertStringToEnum( const  std::string &str, VTKFormat &cod ){
    

    if( str == "ascii"){
        cod = VTKFormat::ASCII ;
        return(true);

    } else if ( str == "appended" ){
        cod = VTKFormat::APPENDED ;
        return(true);

    } else {
        cod = VTKFormat::UNDEFINED ;
        return(false);
    };

};

/*!
 * Converts a std::string as read in vtk file to VTKDataType 
 * @param[in] str string read in DataArray header
 * @param[out]  type VTKDataType to be converted
 * @return  if str contained expected value
 */
bool VTKUtils::convertStringToEnum( const std::string &str, VTKDataType &type ){

    if ( str =="Int8" ){
        type = VTKDataType::Int8;
        return(true);
    }

    else if ( str =="Int16" ){
        type = VTKDataType::Int16;
        return(true);
    }

    else if ( str =="Int32" ){
        type = VTKDataType::Int32;
        return(true);
    }

    else if ( str =="Int64" ){
        type = VTKDataType::Int64;
        return(true);
    }

    else if ( str =="UInt8" ){
        type = VTKDataType::UInt8;
        return(true);
    }

    else if ( str =="UInt16" ){
        type = VTKDataType::UInt16;
        return(true);
    }

    else if ( str =="UInt32" ){
        type = VTKDataType::UInt32;
        return(true);
    }

    else if ( str =="UInt64" ){
        type = VTKDataType::UInt64;
        return(true);
    }

    else if ( str =="Float32" ){
        type = VTKDataType::Float32;
        return(true);
    }

    else if ( str =="Float64" ){
        type = VTKDataType::Float64;
        return(true);
    }

    else {
        type = VTKDataType::UNDEFINED;
        return(false);
    }

};

/*!
 *  Determines the basic VTK type from argument.
 *  @param[in]  type       argument type
 *  @return     basic VTK type 
 *
 */
VTKDataType VTKUtils::whichType( const std::type_info & type ){

    if( type == typeid(int8_t) ){
        return VTKDataType::Int8 ;
    }

    else if( type == typeid(uint8_t) ){
        return VTKDataType::UInt8 ;
    }

    else if( type == typeid(int16_t) ){
        return VTKDataType::Int16 ;
    }

    else if( type == typeid(uint16_t) ){
        return VTKDataType::UInt16 ;
    }

    else if( type == typeid(int32_t) ){
        return VTKDataType::Int32 ;
    }

    else if( type == typeid(uint32_t) ){
        return VTKDataType::UInt32 ;
    }

    else if( type == typeid(uint64_t) ){
        return VTKDataType::Int64 ;
    }

    else if( type == typeid(uint64_t) ){
        return VTKDataType::UInt64 ;
    }

    else if( type == typeid(float) ){
        return VTKDataType::Float32 ;
    }

    else if( type == typeid(double) ){
        return VTKDataType::Float64 ;
    }


    else{
        return VTKDataType::UNDEFINED ;
    };


};

/*!
 * @}
 */
