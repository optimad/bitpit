#include <cassert>

namespace bitpit{

/*!
 *  Add user data for input or output. 
 *  Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 *  @tparam     T        type of data stored in std::vector
 *  @param[in]  name_    name of field
 *  @param[in]  data_    user data to be stored internally; 
 */
template<class T>
VTKField** VTK::addData( std::string name_, std::vector<T> &data_ ){

    VTKField**   ptr0 = this->addData( name_ ) ;
    VTKField*    ptr1 = new VTKFieldWithVector<T>( **ptr0, data_ );

    std::swap( *ptr0, ptr1 );

    delete ptr1 ;
    ptr1 = NULL ;

    return ptr0 ;

}

/*!
 *  Add user data for input or output. 
 *  Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 *  @tparam     T        type of data stored in std::vector
 *  @param[in]  name_    name of field
 *  @param[in]  ftype_   type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 *  @param[in]  loc_     location of data [VTKLocation::CELL/VTKLocation::POINT]
 *  @param[in]  data_    user data to be stored internally; 
 */
template<class T>
VTKField** VTK::addData( std::string name_, VTKFieldType ftype_,  VTKLocation loc_, std::vector<T> &data_ ){

    VTKField** ptr = this->addData( name_, data_) ;
    (*ptr)->setLocation( loc_ );
    (*ptr)->setFieldType( ftype_ );

    return ptr;
}

/*!  
 *  Constructor for mesh of constant element type with "point" and "connectivity" stored in vectors.
 *  sets input parameters and calls default constructor
 *  @tparam     T0    type of "point" information
 *  @tparam     T1    type of "connectivity" information
 *  @param[in]  dir_  Directory of vtk file with final "/"
 *  @param[in]  name_ Name of vtk file without suffix
 *  @param[in]  type_ Type of element in grid
 *  @param[in]  point "point" data; point.size() must be number of points in grid; std::vector<std::vector<T0>> and std::vector<std::array<T0,3>> are supported
 *  @param[in]  connectivity "connectivity" data; connectivity.size() must be number of cells; std::vector<std::vector<T1>> and std::vector<std::array<T1,n>> are supported
 */
template<class T0, class T1>
VTKUnstructuredGrid::VTKUnstructuredGrid( std::string dir_, std::string name_, VTKElementType type_, std::vector<T0> &point, std::vector<T1> &connectivity ):VTKUnstructuredGrid( dir_, name_, type_ ){

    this->setGeomData( point, connectivity ) ;

    return ;

};

/*!
 * Set the geometry data to be handeled internally by linking to external data.
 * @tparam     T0       type of "point" information
 * @param[in]  name     Name of the field ["Points","offsets","types","connectivity"]
 * @param[in]  data     std::vector<T0> with data
 */
template<class T0>
void VTKUnstructuredGrid::setGeomData( std::string name, std::vector<T0> &data ){

    uint8_t     index ;

    if( name == "Points" ){
        index = 0 ; 

    } else if ( name == "offsets" ){
        index = 1 ; 

    } else if ( name == "types" ){
        index = 2 ; 

    } else if ( name == "connectivity" ){
        index = 3 ; 

    } else {
        assert(false);

    };

    VTKField*  ptr = new VTKFieldWithVector<T0>( *geometry[index], data ); 
    std::swap( ptr, geometry[index] );

    delete ptr ;
    ptr = NULL ;

    return;

};

/*!  
 * Set the geometry data ("points" and "connectivity" ) to be handeled internally by linking to external data for homogenous meshes.
 * @tparam     T0    type of "point" information
 * @tparam     T1    type of "connectivity" information
 * @param[in]  point "point" data; point.size() must be number of points in grid; std::vector<std::vector<T0>> and std::vector<std::array<T0,3>> are supported
 * @param[in]  connectivity "connectivity" data; connectivity.size() must be number of cells; std::vector<std::vector<T1>> and std::vector<std::array<T1,n>> are supported
 */
template<class T0, class T1>
void VTKUnstructuredGrid::setGeomData( std::vector<T0> &point, std::vector<T1> &connectivity ){

    this->setGeomData( "Points", point );
    this->setGeomData( "connectivity", connectivity );

    int nNode = point.size() ;
    int nCell = connectivity.size() ;

    setDimensions( nCell, nNode ) ;

    return ;

};

}
