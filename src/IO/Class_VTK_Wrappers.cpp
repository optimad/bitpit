
#include"Class_VTK_Wrappers.hpp"

/*! ========================================================================================
 * @ingroup     VisualizationToolKit
 * @{
 * @class       VtkUnstrVec
 * @brief       VtkUnstrVec implements an generic interface for VTK_UnstructuredGrid for std::vector<> containers
 *
 * VTKUnstrVec is a final class for reading and  writing unstructured grids when the grid and user data are stored in std::vector<> container.
 * The supported types are defined in BVector and BArray and merged in SAll.
 * The grid must be composed of uniform elements, e.g. only thetraeda.
 *
 */


/*!====================================================================================================
 * Default constructor
 */
VtkUnstrVec::VtkUnstrVec():VTK_UnstructuredGrid<VtkUnstrVec>(){
};

/*!====================================================================================================
 * Constructor for reading VTU files. Sets the given input parameters
 * @param[in]   dir_    directory of vtu file with final "/"
 * @param[in]   name_   name of vtu file withstd::out suffix
 * @param[in]   codex_  codex used in file ["appended"/"ascii"]
 * @param[in]   type_   element type of unstructured grid. See http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
VtkUnstrVec::VtkUnstrVec( std::string dir_, std::string name_, std::string codex_, uint8_t type_ ):VTK_UnstructuredGrid<VtkUnstrVec>( ){

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

};

/*!====================================================================================================
 * Destructor
 */
VtkUnstrVec::~VtkUnstrVec(){
};

/*!====================================================================================================
 * Finds a field within class by its name.
 * @param[in]    name_   name of the field
 * @param[std::out]  the_field   pointer to the field if found, else unaltered
 * @return      true if found, else false   
 */
bool VtkUnstrVec::GetFieldByName( const std::string &name_, VtkUnstrVec::ufield *&the_field ){

    std::vector<VtkUnstrVec::ufield>::iterator    it_ ;

    for( it_=adata.begin(); it_!=adata.end(); ++it_){
        if( (*it_).name == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    return false ;
};

/*!====================================================================================================
 * Writes the vtu file. Checks if all necessary information has been provided.
 * Will write only data which has been provided by VtkUnstrVec::AddData() or VtkUnstrVec::LinkData().
 */
void VtkUnstrVec::Write(  ) {

    bool                        CanWrite(true) ;
    ufield*                     FPtr ; 
    std::string                 name ;
    std::vector<std::string>    TBD;

    CanWrite = CanWrite && GetFieldByName( "Points", FPtr ) ; 
    CanWrite = CanWrite && GetFieldByName( "connectivity", FPtr ) ; 

    if( CanWrite){

        for( unsigned i=0; i<nr_data; ++i){
            name = data[i].GetName() ;
            if( !GetFieldByName(name,FPtr) ){
                TBD.push_back(name) ;
            };
        };

        for( unsigned i=0; i<TBD.size(); ++i){
            VTK::RemoveData( TBD[i] ) ;
        };

        VTK::Write() ;
    }

    else{
        std::cout << " VtkUnstrVec::Write cannot write file since \"Points\" or \"connectivity\" vectors are missing" << std::endl ;

    };

    return;
};

/*!  ===================================================================================
 *  CRPT interface for writing data to stream. Uses static visitor pattern.
 *  @param[in]  str         stream to write to
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be written. Either user data or grid data
 */
void VtkUnstrVec::Flush( std::fstream &str, std::string codex, std::string name ) {


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
        uint64_t off_(0) ;
        for(uint64_t  n=0; n<nr_cells-1; n++) {
            off_ += NumberOfElements( type ) ; 
            flush_ascii( str, off_  ) ;
            str << std::endl ;
        };
        off_ += NumberOfElements( type ) ;
        flush_ascii( str, off_  ) ;
    }

    else if( codex == "binary" && name == "offsets"){
        uint64_t off_(0), nT( NumberOfElements( type ) ) ;
        for( uint64_t n=0; n<nr_cells; n++) {
            off_ += nT ; 
            flush_binary( str, off_  ) ;
        };
    }

    else{

        ufield  *f_ ;

        if( GetFieldByName( name, f_)){
            stream_visitor  visitor;
            visitor.SetStream( str ) ;
            visitor.SetCodex( codex ) ;
            visitor.SetName( name ) ;
            visitor.SetTask( "write" ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };

    };

    return ;

};

/*!  ===================================================================================
 *  CRPT interface for reading data from stream. Uses static visitor pattern.
 *  @param[in]  str         stream to write to
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be written. Either user data or grid data
 */
void VtkUnstrVec::Absorb( std::fstream &str, std::string codex, std::string name ) {

    ufield  *f_ ;

    if( GetFieldByName( name, f_)){

        VTK::Field_C*   FPtr ;
        if( VTK::GetFieldByName(name, FPtr) ){

            stream_visitor  visitor;
            visitor.SetStream( str ) ;
            visitor.SetCodex( codex ) ;
            visitor.SetName( name ) ;
            visitor.SetTask( "read" ) ;

            visitor.SetSize( FPtr->GetElements() ) ;
            visitor.SetComponents( FPtr->GetComponents() ) ;

            if( name == "connectivity") visitor.SetComponents( NumberOfElements(type) ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };
    }

    return ;

};


/*!  ===================================================================================
 * @struct      VtkUnstrVec::ufield
 * @brief       ufield stores the name and a pointer to field data, both geometry and user data
 */

/*!  ===================================================================================
 * Default constructor.
 */
VtkUnstrVec::ufield::ufield( ){
};

/*!  ===================================================================================
 * Sets the file stream
 * @param[in]   str_    file stream to be used.
 */
void VtkUnstrVec::stream_visitor::SetStream( std::fstream& str_){
    str = &str_ ;
};

/*!  ===================================================================================
 * Sets the file codex
 * @param[in]   codex_    codex_ to be used ["ascii"/"appended"]
 */
void VtkUnstrVec::stream_visitor::SetCodex( std::string codex_){
    codex = codex_ ;
};

/*!  ===================================================================================
 * Sets the task
 * @param[in]   task_   sets the task ["read"/"write"]
 */
void VtkUnstrVec::stream_visitor::SetTask( std::string task_){
    task = task_ ;
};

/*!  ===================================================================================
 * Sets the name of the field which should be visited.
 * In case name="connectivity" special treatment is performed 
 * @param[in]   name_   name of the field 
 */
void VtkUnstrVec::stream_visitor::SetName( std::string name_){
    name = name_ ;
};

/*!  ===================================================================================
 * Sets the size of the filed
 * @param[in]   size_    size of the entire field
 */
void VtkUnstrVec::stream_visitor::SetSize( uint64_t size_){
    size = size_ ;
};

/*!  ===================================================================================
 * Sets the number of components of the field
 * @param[in]   com_    numer of components
 */
void VtkUnstrVec::stream_visitor::SetComponents( uint8_t com_){
    components = com_ ;
};

/*!
 @}
*/

