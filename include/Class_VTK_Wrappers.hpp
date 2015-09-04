#ifndef __CLASS_VTK_WRAP_HH__
#define __CLASS_VTK_WRAP_HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/variant.hpp>

#include"Class_VTK.hpp"

using namespace std;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// double array
typedef array< double,3 >              darray3E;
typedef vector< darray3E >             dvecarr3E;

typedef boost::mpl::vector<
    vector<int8_t>*  ,
    vector<int16_t>* ,
    vector<int32_t>* ,
    vector<int64_t>* ,
    vector<uint8_t>* ,
    vector<uint16_t>*, 
    vector<uint32_t>*, 
    vector<uint64_t>*
    > SInt ;

typedef boost::mpl::vector<
    vector<float>*   ,
    vector<double>*  
    >::type SFloat;

typedef boost::mpl::vector<
    vector< array<int16_t,3> >*  , 
    vector< array<int32_t,3> >*  , 
    vector< array<int64_t,3> >*  , 
    vector< array<uint16_t,3> >* ,
    vector< array<uint32_t,3> >* ,
    vector< array<uint64_t,3> >* ,

    vector< array<int16_t,4> >*  ,
    vector< array<int32_t,4> >*  ,
    vector< array<int64_t,4> >*  ,
    vector< array<uint16_t,4> >* ,
    vector< array<uint32_t,4> >* ,
    vector< array<uint64_t,4> >* ,

    vector< array<int16_t,8> >*  ,
    vector< array<int32_t,8> >*  ,
    vector< array<int64_t,8> >*  ,
    vector< array<uint16_t,8> >* ,
    vector< array<uint32_t,8> >* ,
    vector< array<uint64_t,8> >* 
    >::type SIArray;


typedef boost::mpl::vector<
    vector< array<float,3> >*    ,
    vector< array<double,3> >*   
    >::type SFArray;

typedef boost::mpl::vector<
    vector< vector<int16_t> >*  , 
    vector< vector<int32_t> >*  , 
    vector< vector<int64_t> >*  , 
    vector< vector<uint16_t> >* ,
    vector< vector<uint32_t> >* ,
    vector< vector<uint64_t> >* 
    >::type SIVector ;

typedef boost::mpl::vector<
    vector< vector<float> >*  , 
    vector< vector<double> >*
    >::type SFVector ;

typedef boost::mpl::copy<
    SInt::type, 
    boost::mpl::back_inserter<SFloat> 
    >::type SSingle ;

typedef boost::mpl::copy<
    SIArray::type, 
    boost::mpl::back_inserter<SFArray> 
    >::type SArray ;

typedef boost::mpl::copy<
    SIVector::type, 
    boost::mpl::back_inserter<SFVector> 
    >::type SVector ;

typedef boost::mpl::copy<
    SVector::type, 
    boost::mpl::back_inserter<SArray> 
    >::type SMulti ;

typedef boost::mpl::copy<
    SSingle::type, 
    boost::mpl::back_inserter<SMulti> 
    >::type SAll ;

typedef boost::make_variant_over< SAll >::type bv ;

struct stream_visitor : public boost::static_visitor<>{

    private:
    fstream*     str ;
    uint8_t      components ;
    uint64_t     size ;
    string       codex ;
    string       name ;
    string       task ;

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
    void operator()( vector< vector<T> >* t) const{

        if( task=="write" && codex=="ascii")  flush_ascii( *str, 1, *t) ;

        if( task=="write" && codex=="binary") flush_binary( *str, *t) ;

        if( task=="read"  && codex=="ascii" )  {
            if( name == "connectivity"){
                (*t).resize( size/components , vector<T>( components , 0 ) ) ;
            }

            else{
                (*t).resize( size , vector<T>( components , 0 ) ) ;
            };

            absorb_ascii( *str, *t) ;
        };

        if( task=="read"  && codex=="binary") {
            if( name == "connectivity"){
                (*t).resize( size/components , vector<T>( components , 0 ) ) ;
            }

            else{
                (*t).resize( size , vector<T>( components , 0 ) ) ;
            };
            absorb_binary( *str, *t) ;
        };
    };


    void SetStream( fstream& str_){
        str = &str_ ;
    };

    void SetCodex( string codex_){
        codex = codex_ ;
    };

    void SetTask( string task_){
        task = task_ ;
    };

    void SetName( string name_){
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
        VTK::Field_C*   FPtr;
        bv              DPtr ;
    };

    uint8_t         type ;    
    vector<ufield>   adata;

    void Flush( fstream &str, string codex_, string name ) ;
    void Absorb( fstream &str, string codex_, string name ) ;
    
    public:
    VtkUnstrVec( ) ;

    template< class T0, class T1>
    VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_, vector<T0> &points_ext, vector<T1> &connectivity_external ) ;
    VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_ ) ;


   ~VtkUnstrVec( ) ;


    template< class T>
    void    AddData( vector<T> &data, string name_, string loc_ ) ; 
        
    template< class T>
    void    AddData( vector< array<T,3> > &data, string name_, string loc_ ) ; 

    template< class T>
    void    AddData( vector< vector<T> > &data, string name_, string loc_ ) ; 

    private:
    bool    GetFieldByName( const string &name_, ufield*& the_field ) ;



   

};


#include"Class_VTK_Wrappers.tpp"

#endif
