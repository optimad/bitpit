#include "VTK.hpp"

/*!
 * @ingroup    VTKUtils
 * @{
 *
 * @brief Utility fuctions for VTK
 */

/*!
 * Calculates the size in bytes of basic types supported by VTK
 * @param[in]   type       basic types supported by VTK [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
 * \return      size of basic type
 */
uint8_t VTKUtils::sizeOfType( std::string type ){

  uint8_t nbytes ;

  if( type == "Int8" || type == "UInt8"){
    nbytes = sizeof(int8_t) ;
  }

  if( type == "Int16" || type == "UInt16"){
    nbytes = sizeof(int16_t) ;
  }

  else if( type == "Int32" || type == "UInt32"){
    nbytes = sizeof(int32_t) ;
  }

  else if( type == "Int64" || type == "UInt64"){
    nbytes = sizeof(int64_t) ;
  }

  else if( type == "Float32" ){
    nbytes = sizeof(float) ;
  }

  else if( type == "Float64" ){
    nbytes = sizeof(double) ;
  };

  return nbytes ;
};

/*!  
 *  Returns number of vertices for a given element type.
 *  Codification of element type according to http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *  @param[in]  t           element type
 *  \return     number of vertices of element type
 */
uint8_t VTKUtils::getNNodeInElement( uint8_t t){

  int e;

  switch (t){

    case 1: e= 1; break;
    case 3: e= 2; break;
    case 5: case 21: e =3; break;
    case 8: case 9:  case 10: e= 4; break;
    case 11: case 12: case 23: case 24: e=8; break;
    case 13: case 22: e= 6; break;
    case 14: e=5; break;
    case 25: e=20; break;
    default:
      std::cout << "Element type not supported: " << t << std::endl;

  };

  return e;

};


/*!
 * Converts a string into a Field information.
 * @param[in]   line_       string to be converted
 * @param[out]  data_       Field information
 * \return      true if successful
 */
bool VTKUtils::convertStringToDataArray( std::string &line_, VTKField &data_  ){

    std::string type_, name_, code_, comp_, offs_ ;
    int    components_(1), offset_ ;
    
    bool  success(true) ;
    
    
    if( Keyword_In_String( line_, "<DataArray ") ){  
        success = success && Get_After_Keyword( line_, "type=", '\"', type_) ;
        success = success && Get_After_Keyword( line_, "Name=", '\"', name_) ;
        success = success && Get_After_Keyword( line_, "format=", '\"', code_) ;
    
        if( Get_After_Keyword( line_, "NumberOfComponents=", '\"', comp_)  ){
            convert_string( comp_, components_ ) ;
        };
    
        data_.setType(type_) ;
        data_.setName(name_) ;
        data_.setComponents(components_) ;
        data_.setCodification(code_) ;
    
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
 * @param[out]   str        string in VTK format
 * @param[in]  field_       Field information
 * \return      true if successful
 */
void  VTKUtils::convertDataArrayToString( std::string &str, VTKField &field_ ){

  std::stringstream  os("") ;

  os << "        <DataArray "
       << "type=\"" << field_.getType() << "\" "
       << "Name=\"" << field_.getName() << "\" "
       << "NumberOfComponents=\""<< unsigned(field_.getComponents()) << "\" "
       << "format=\"" << field_.getCodification() << "\" ";

  if( field_.getCodification() == "appended"){
    os << "offset=\"" << field_.getOffset() << "\" " ;
  };
  
  os << ">" ;

  str = os.str() ;       
       
  return ;
  
};

/*!
 * Converts a parallel field information to string as requested by VTK format.
 * @param[out]   str        string in VTK format
 * @param[in]  field_       Field information
 * \return      true if successful
 */
void  VTKUtils::convertPDataArrayToString( std::string &str, VTKField &field_ ){

  std::stringstream  os("") ;

  os << "        <PDataArray "
      << "type=\"" << field_.getType() << "\" "
      << "Name=\"" << field_.getName() << "\" "
      << "NumberOfComponents=\""<< unsigned(field_.getComponents()) << "\" " 
      << ">" ;

  str = os.str() ;

  return ;
  
};

/*!
 * @}
 */
