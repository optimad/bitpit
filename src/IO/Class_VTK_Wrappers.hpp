#ifndef __CLASS_VTK_WRAP_HH__
#define __CLASS_VTK_WRAP_HH__

#include <boost/mpl/vector.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/variant.hpp>

#include"Class_VTK.hpp"


// integer std::vectors
typedef std::vector< int >                  ivector1D;
typedef std::vector< ivector1D >            ivector2D;
typedef std::vector< ivector2D >            ivector3D;
typedef std::vector< ivector3D >            ivector4D;

// double vectors
typedef std::vector< double >               dvector1D;
typedef std::vector< dvector1D >            dvector2D;
typedef std::vector< dvector2D >            dvector3D;
typedef std::vector< dvector3D >            dvector4D;

// double array
typedef std::array< double,3 >              darray3E;
typedef std::vector< darray3E >             dvecarr3E;

typedef boost::mpl::vector<
    std::vector<int8_t>*  ,
    std::vector<int16_t>* ,
    std::vector<int32_t>* ,
    std::vector<int64_t>* ,
    std::vector<uint8_t>* ,
    std::vector<uint16_t>*, 
    std::vector<uint32_t>*, 
    std::vector<uint64_t>*,

    std::vector<float>*   ,
    std::vector<double>*  ,

    std::vector< std::vector<int8_t> >*  , 
    std::vector< std::vector<int16_t> >*  , 
    std::vector< std::vector<int32_t> >*  , 
    std::vector< std::vector<int64_t> >*  , 
    std::vector< std::vector<uint8_t> >*  , 
    std::vector< std::vector<uint16_t> >* ,
    std::vector< std::vector<uint32_t> >* ,
    std::vector< std::vector<uint64_t> >* ,

    std::vector< std::vector<float> >*  , 
    std::vector< std::vector<double> >*
    >::type BVector;

typedef boost::mpl::vector<
    std::vector< std::array<int16_t,3> >*  , 
    std::vector< std::array<int32_t,3> >*  , 
    std::vector< std::array<int64_t,3> >*  , 
    std::vector< std::array<uint16_t,3> >* ,
    std::vector< std::array<uint32_t,3> >* ,
    std::vector< std::array<uint64_t,3> >* ,

    std::vector< std::array<int16_t,4> >*  ,
    std::vector< std::array<int32_t,4> >*  ,
    std::vector< std::array<int64_t,4> >*  ,
    std::vector< std::array<uint16_t,4> >* ,
    std::vector< std::array<uint32_t,4> >* ,
    std::vector< std::array<uint64_t,4> >* ,

    std::vector< std::array<int16_t,8> >*  ,
    std::vector< std::array<int32_t,8> >*  ,
    std::vector< std::array<int64_t,8> >*  ,
    std::vector< std::array<uint16_t,8> >* ,
    std::vector< std::array<uint32_t,8> >* ,
    std::vector< std::array<uint64_t,8> >* ,

    std::vector< std::array<float,3> >*    ,
    std::vector< std::array<double,3> >*   
    >::type BArray;

typedef boost::mpl::copy<
    BVector::type, 
    boost::mpl::back_inserter<BArray> 
    >::type SAll ;

typedef boost::make_variant_over< SAll >::type bv ;

struct stream_visitor : public boost::static_visitor<>{

    private:
    fstream*     str ;
    uint8_t      components ;
    uint64_t     size ;
    std::string       codex ;
    std::string       name ;
    std::string       task ;

    public:

    template <typename T>
    void operator()( T* t) const{
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

    template <typename T>
    void operator()( std::vector< std::vector<T> >* t) const{

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


    void SetStream( fstream& str_){
        str = &str_ ;
    };

    void SetCodex( std::string codex_){
        codex = codex_ ;
    };

    void SetTask( std::string task_){
        task = task_ ;
    };

    void SetName( std::string name_){
        name = name_ ;
    };

    void SetSize( uint64_t size_){
        size = size_ ;
    };

    void SetComponents( uint8_t com_){
        components = com_ ;
    };


};
 

class VtkUnstrVec : public VTK_UnstructuredGrid<VtkUnstrVec>{


    friend VTK_UnstructuredGrid<VtkUnstrVec> ;

    private:
    typedef  VTK_UnstructuredGrid<VtkUnstrVec> Base ;

    struct ufield{
        std::string          name;
        bv              DPtr ;

        template<class T>
        ufield( std::string name_, std::vector<T>& data_): name(name_), DPtr(&data_){} ;
        ufield(){} ;
    };

    uint8_t         type ;    
    std::vector<ufield>   adata;

    void Flush( fstream &str, std::string codex_, std::string name ) ;
    void Absorb( fstream &str, std::string codex_, std::string name ) ;
    
    public:
    VtkUnstrVec( ) ;

    template< class T0, class T1>
    VtkUnstrVec( std::string dir_, std::string name_, std::string codex_, uint8_t type_, std::vector<T0> &points_ext, std::vector<T1> &connectivity_external ) ;
    VtkUnstrVec( std::string dir_, std::string name_, std::string codex_, uint8_t type_ ) ;

   ~VtkUnstrVec( ) ;

    template< class T>
    void    AddData( std::vector<T> &data, std::string name_, std::string loc_ ) ; 
        
    template< class T>
    void    AddData( std::vector< std::array<T,3> > &data, std::string name_, std::string loc_ ) ; 

    template< class T>
    void    AddData( std::vector< std::vector<T> > &data, std::string name_, std::string loc_ ) ; 

    template< class T>
    void    LinkData( std::vector<T> &data, std::string name_ ) ; 

    void    Write() ;
 
    private:
    bool    GetFieldByName( const std::string &name_, ufield*& the_field ) ;


};

#include"Class_VTK_Wrappers.tpp"

#endif
