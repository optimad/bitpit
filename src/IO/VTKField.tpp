
namespace bitpit{

/*!
 * @ingroup VisualizationToolKit
 * @{
 */

/*!
 * @class VTKFieldWithVector
 * @brief A template class for VTKFields with data stored is a std::vector<T>
 */

/*!
 *  Destructor. 
 *  @tparam   T  type of the data field
 */
template<class T>
VTKFieldWithVector<T>::~VTKFieldWithVector( ){

    m_ptr = NULL ;
};

/*!
 *  Default constructor. 
 *  @tparam   T  type of the data field
 */
template<class T>
VTKFieldWithVector<T>::VTKFieldWithVector( ): VTKField( ) {

    m_ptr = NULL ;
};

/*!
 *  Copy constructor from VTKField and assignement of data
 *  @tparam   T  type of the data field
 *  @param[in]   other  VTKField to be copied
 *  @param[in]   data_   vector with field data 
 */
template<class T>
VTKFieldWithVector<T>::VTKFieldWithVector( const VTKField &other, std::vector<T> &data_ ): VTKField( other ) {

    this->setData( data_ ) ;
};

/*!
 * Constructor. 
 * @tparam   T  type of the data field
 * @param[in]   name_   name of data field
 * @param[in]   data_   vector with field data 
 */
template<class T>
VTKFieldWithVector<T>::VTKFieldWithVector( std::string name_, std::vector<T> &data_ ): VTKField( name_) {

    this->setData( data_ ) ;

};

template<class T>
void VTKFieldWithVector<T>::setData( std::vector<T> &data_ ){

    derived = true ;
    m_ptr = &data_ ;
    this->setDataType( VTKTypes::whichType(data_) );

};

template<class T>
void VTKFieldWithVector<T>::flushData( std::fstream &str) const {

    if( this->getCodification() == VTKFormat::APPENDED){
        genericIO::flushBINARY(str, *m_ptr);

    } else {
        genericIO::flushASCII(str, *m_ptr);
    }
};

template< class T>
//template< class U, typename std::enable_if< bitpit::is_vector< U >::value  && std::is_same<T,U>::value>::type* >
void VTKFieldWithVector<T>::absorbData( std::fstream &str) const {

    int n;

    m_ptr->resize( this->getElements() );

    if( this->getFieldType() == VTKFieldType::SCALAR || this->getFieldType() == VTKFieldType::VECTOR ){
        n = static_cast<int>( this->getFieldType() );

    } else if ( this->getFieldType() == VTKFieldType::CONSTANT ){
        n = static_cast<int>( this->getComponents() );

    };


    for( auto & element : (*m_ptr) ){
        //element.resize(n);
        bitpit::vtk::allocate( element, n );
    }

    if( this->getCodification() == VTKFormat::APPENDED){
        genericIO::absorbBINARY(str, *m_ptr);
    } else {
        genericIO::absorbASCII(str, *m_ptr);
    }
};

//template< class T>
//template< class U, typename std::enable_if< ! bitpit::is_vector< U >::value  && std::is_same<T,U>::value >::type* >
//void VTKFieldWithVector<T>::absorbData( std::fstream &str) const {
//
//    m_ptr->resize( this->getElements() );
//
//    if( this->getCodification() == VTKFormat::APPENDED){
//        genericIO::absorbBINARY(str, *m_ptr);
//    } else {
//        genericIO::absorbASCII(str, *m_ptr);
//    }
//};

/*!
 * @}
 */

}
