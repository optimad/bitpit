#include"Class_VTK_Wrappers.hpp"

/*! ====================================================================================================
 * Constructor suited for writing. Sets the provided input information and calls default constructor
 * @tparam      T0      class used to store grid nodes
 * @tparam      T1      class used to store grid connectivity
 * @param[in]   dir_    directory of VTU file with final "/"
 * @param[in]   name_   name of VTU file withstd::out suffix
 * @param[in]   codex_  codex used for writing
 * @param[in]   type_   element type of unstructured grid. See http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 * @param[in]   points_ext  grid nodes; the size of the vector is used to determine the number of nodes in grid;
 * @param[in]   connectivity_ext   grid cell-node connectivity; the size of the vector is used to determine the number of cells in grid;
 */
template< class T0, class T1>
VtkUnstrVec::VtkUnstrVec( std::string dir_, std::string name_, std::string codex_, uint8_t type_, std::vector<T0> &points_ext, std::vector<T1> &connectivity_ext ) :VTK_UnstructuredGrid<VtkUnstrVec>( ){


    int ncells_, npoints_, nconn_;

    T0  dum0 ;
    T1  dum1 ;

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

    ncells_ = connectivity_ext.size() ;
    npoints_ = points_ext.size() ;

    nconn_ = ncells_ * NumberOfElements( type ) ;

    SetGeomTypes( WhichType(dum0), "UInt64", "UInt8", WhichType(dum1)  ) ;
    SetDimensions( ncells_, npoints_, nconn_ ) ;

    adata.resize(2) ;

    adata[0].DPtr = &points_ext ;
    adata[0].name = "Points" ;

    adata[1].DPtr = &connectivity_ext ;
    adata[1].name = "connectivity" ;

};


/*! ====================================================================================================
 * Adds scalar field data for std::output stored in std::vector< T >
 * @tparam      T       POD type 
 * @param[in]   data_   data for std::output
 * @param[in]   name_   name of data file
 * @param[in]   loc_    location of the data
 */
template<class T>
void VtkUnstrVec::AddData( std::vector<T> &data_, std::string name_, std::string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 1, WhichType(dum_), loc_ ) ;

    return;
};

/*! ====================================================================================================
 * Adds vector field data for std::output stored in std::vector< std::array<T,3> >
 * @tparam      T       POD type 
 * @param[in]   data_   data for std::output
 * @param[in]   name_   name of data file
 * @param[in]   loc_    location of the data
 */
template<class T>
void VtkUnstrVec::AddData( std::vector< std::array<T,3> > &data_, std::string name_, std::string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;

    return;
};

/*! ====================================================================================================
 * Adds vector field data for std::output stored in std::vector< std::vector<T>>.
 * Assumes that the size of each internal vector is equal 3.
 * @tparam      T       POD type 
 * @param[in]   data_   data for std::output
 * @param[in]   name_   name of data file
 * @param[in]   loc_    location of the data
 */
template<class T>
void VtkUnstrVec::AddData( std::vector< std::vector<T> > &data_, std::string name_, std::string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;


    return;
};

/*! ====================================================================================================
 * Adds field data for input.
 * @tparam      T       POD type 
 * @param[in]   data_   reference to data structure 
 * @param[in]   name_   name of data file. Should correspond to the name in the VTU file.
 */
template<class T>
void VtkUnstrVec::LinkData( std::vector<T> &data_, std::string name_ ){

    adata.push_back( ufield(name_,data_) ) ;

    return;
};


/*! ====================================================================================================
 * Constructor. Sets the input variables.
 * @tparam      T       POD type
 * @param[in]   name_   name of data file. 
 * @param[in]   data_   reference to data structure 
 */
template<class T>
VtkUnstrVec::ufield::ufield( std::string name_, std::vector<T>& data_): name(name_), DPtr(&data_){
} ;

/*! ====================================================================================================
 * @struct      VtkUnstrVec::stream_visitor
 * @brief       static visitor pattern used for reading and writing different data types
 *  stream_visitor distinguishes different data types and codexes to be read or written and calls the correct low level function
 */

/*! ====================================================================================================
 * opertor() for std::vector<> and std::vector<str::array<>> used in boost::apply_visitor. 
 * Distinguishes task ["read"/"write"] and codex ["ascii"/"binary"].
 * If task=="read" the data structure is resized.
 * @tparam      T       container type
 */
template <typename T>
void VtkUnstrVec::stream_visitor::operator()( T* t) const{
    if( task=="write" && codex=="ascii")  flush_ascii( *str, 1, *t) ;
    if( task=="write" && codex=="binary") flush_binary( *str, *t) ;
    if( task=="read"  && codex=="ascii")  {
        (*t).resize(size); 
        absorb_ascii( *str, *t) ;
    };
    if( task=="read"  && codex=="binary") {
        (*t).resize(size);
        absorb_binary( *str, *t) ;
    };
};

/*! ====================================================================================================
 * opertor() for std::vector<str::vector<>> used in boost::apply_visitor. 
 * Distinguishes task ["read"/"write"] and codex ["ascii"/"binary"].
 * If task=="read" the data structure is resized.
 * @tparam      T       container type
 */
template <typename T>
void VtkUnstrVec::stream_visitor::operator()( std::vector< std::vector<T> >* t) const{

    if( task=="write" && codex=="ascii")  flush_ascii( *str, 1, *t) ;

    if( task=="write" && codex=="binary") flush_binary( *str, *t) ;

    if( task=="read"  && codex=="ascii" )  {
        if( name == "connectivity"){
            (*t).resize( size/components , std::vector<T>( components , 0 ) ) ;
        }

        else{
            (*t).resize( size , std::vector<T>( components , 0 ) ) ;
        };

        absorb_ascii( *str, *t) ;
    };

    if( task=="read"  && codex=="binary") {
        if( name == "connectivity"){
            (*t).resize( size/components , std::vector<T>( components , 0 ) ) ;
        }

        else{
            (*t).resize( size , std::vector<T>( components , 0 ) ) ;
        };
        absorb_binary( *str, *t) ;
    };
};


