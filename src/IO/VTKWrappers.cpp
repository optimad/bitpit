
#include "VTKWrappers.hpp"

/*!
 * @ingroup     VisualizationToolKit
 * @{
 * @class       VTKUnstructuredVec
 * @brief       VTKUnstructuredVec implements an generic interface for VTK_UnstructuredGrid for std::vector<> containers
 *
 * VTKUnstrVec is a final class for reading and  writing unstructured grids when the grid and user data are stored in std::vector<> container.
 * The supported types are defined in BVector and BArray and merged in SAll.
 * The grid must be composed of uniform elements, e.g. only thetraeda.
 *
 */


/*!
 * Default constructor
 */
VTKUnstructuredVec::VTKUnstructuredVec():VTKUnstructuredGrid<VTKUnstructuredVec>(){
};

/*!
 * Constructor for reading VTU files. sets the given input parameters
 * @param[in]   dir_    directory of vtu file with final "/"
 * @param[in]   name_   name of vtu file withstd::out suffix
 * @param[in]   codex_  codex used in file ["appended"/"ascii"]
 * @param[in]   type_   element type of unstructured grid. See http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
VTKUnstructuredVec::VTKUnstructuredVec( std::string dir_, std::string name_, std::string codex_, uint8_t type_ ):VTKUnstructuredGrid<VTKUnstructuredVec>( ){

    setNames( dir_, name_ );
    setCodex( codex_ ) ;

    type = type_ ;

};

/*!
 * Destructor
 */
VTKUnstructuredVec::~VTKUnstructuredVec(){
};

/*!
 * Finds a field within class by its name.
 * @param[in]    name_   name of the field
 * @param[std::out]  the_field   pointer to the field if found, else unaltered
 * @return      true if found, else false   
 */
bool VTKUnstructuredVec::getFieldByName( const std::string &name_, VTKUnstructuredVec::ufield *&the_field ){


    for(  auto &field : adata ){
        if( field.name == name_ ){
            the_field = &field ;
            return true ;
        };
    };

    return false ;
};

/*!
 * Writes the vtu file. Checks if all necessary information has been provided.
 * Will write only data which has been provided by VTKUnstructuredVec::AddData() or VTKUnstructuredVec::LinkData().
 */
void VTKUnstructuredVec::write(  ) {

    bool                        CanWrite(true) ;
    ufield*                     FPtr ; 
    std::string                 name ;
    std::vector<std::string>    TBD;

    CanWrite = CanWrite && getFieldByName( "Points", FPtr ) ; 
    CanWrite = CanWrite && getFieldByName( "connectivity", FPtr ) ; 

    if( CanWrite){

        for( auto &field : data ){
            name = field.getName() ;
            if( !getFieldByName(name,FPtr) ){
                TBD.push_back(name) ;
            };
        };

        for( auto &field : TBD ){
            VTK::removeData( field ) ;
        };

        VTK::write() ;
    }

    else{
        std::cout << " VTKUnstructuredVec::write() cannot write file since \"Points\" or \"connectivity\" vectors are missing" << std::endl ;

    };

    return;
};

/*!  
 *  CRPT interface for writing data to stream. Uses static visitor pattern.
 *  @param[in]  str         stream to write to
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be written. Either user data or grid data
 */
void VTKUnstructuredVec::flush( std::fstream &str, std::string codex, std::string name ) {


    if( codex == "ascii" && name == "types"){
        for( uint64_t n=0; n<nr_cells-1; n++) {
            flush_ascii( str, type  ) ;
            str << std::endl ;
        };
        flush_ascii( str, type  ) ;
    }

    else if( codex == "binary" && name == "types"){
        for( uint64_t n=0; n<nr_cells; n++) flush_binary( str, type  ) ;
    }

    else if( codex == "ascii" && name == "offsets"){
        uint64_t off_(0), nT( VTKUtils::getNNodeInElement( type ) ) ;
        for(uint64_t  n=0; n<nr_cells-1; n++) {
            off_ += nT ;
            flush_ascii( str, off_  ) ;
            str << std::endl ;
        };
        off_ += nT ;
        flush_ascii( str, off_  ) ;
    }

    else if( codex == "binary" && name == "offsets"){
        uint64_t off_(0), nT( VTKUtils::getNNodeInElement( type ) ) ;
        for( uint64_t n=0; n<nr_cells; n++) {
            off_ += nT ; 
            flush_binary( str, off_  ) ;
        };
    }

    else{

        ufield  *f_ ;

        if( getFieldByName( name, f_)){
            stream_visitor  visitor;
            visitor.setStream( str ) ;
            visitor.setCodex( codex ) ;
            visitor.setName( name ) ;
            visitor.setTask( "write" ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };

    };

    return ;

};

/*!  
 *  CRPT interface for reading data from stream. Uses static visitor pattern.
 *  @param[in]  str         stream to write to
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be written. Either user data or grid data
 */
void VTKUnstructuredVec::absorb( std::fstream &str, std::string codex, std::string name ) {

    ufield  *f_ ;

    if( getFieldByName( name, f_)){

        VTKField*   FPtr ;
        if( VTK::getFieldByName(name, FPtr) ){

            stream_visitor  visitor;
            visitor.setStream( str ) ;
            visitor.setCodex( codex ) ;
            visitor.setName( name ) ;
            visitor.setTask( "read" ) ;

            visitor.setSize( FPtr->getElements() ) ;
            visitor.setComponents( FPtr->getComponents() ) ;

            if( name == "connectivity") visitor.setComponents( VTKUtils::getNNodeInElement(type) ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };
    }

    return ;

};

/*!  
 * @struct      VTKUnstructuredVec::ufield
 * @brief       ufield stores the name and a pointer to field data, both geometry and user data
 */

/*!  
 * Default constructor.
 */
VTKUnstructuredVec::ufield::ufield( ){
};

/*!  
 * sets the file stream
 * @param[in]   str_    file stream to be used.
 */
void VTKUnstructuredVec::stream_visitor::setStream( std::fstream& str_){
    str = &str_ ;
};

/*!  
 * sets the file codex
 * @param[in]   codex_    codex_ to be used ["ascii"/"appended"]
 */
void VTKUnstructuredVec::stream_visitor::setCodex( std::string codex_){
    codex = codex_ ;
};

/*!  
 * sets the task
 * @param[in]   task_   sets the task ["read"/"write"]
 */
void VTKUnstructuredVec::stream_visitor::setTask( std::string task_){
    task = task_ ;
};

/*!  
 * sets the name of the field which should be visited.
 * In case name="connectivity" special treatment is performed 
 * @param[in]   name_   name of the field 
 */
void VTKUnstructuredVec::stream_visitor::setName( std::string name_){
    name = name_ ;
};

/*!  
 * sets the size of the filed
 * @param[in]   size_    size of the entire field
 */
void VTKUnstructuredVec::stream_visitor::setSize( uint64_t size_){
    size = size_ ;
};

/*!  
 * sets the number of components of the field
 * @param[in]   com_    numer of components
 */
void VTKUnstructuredVec::stream_visitor::setComponents( uint8_t com_){
    components = com_ ;
};

/*!
 @}
*/
