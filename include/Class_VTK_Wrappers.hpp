#ifndef __CLASS_VTK_WRAP_HH__
#define __CLASS_VTK_WRAP_HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

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

class VtkUnstrVec : public VTK_UnstructuredGrid<VtkUnstrVec>{


    friend VTK_UnstructuredGrid<VtkUnstrVec> ;

    private:
    typedef  VTK_UnstructuredGrid<VtkUnstrVec> Base ;

    struct sfield{
        string      name ;     
        string      location ;     
        dvector1D   *data ;
    };

    struct vfield{
        string      name ;     
        string      location ;     
        dvecarr3E   *data ;
    };
    
    string          codex ;
    int             type ;    
    dvecarr3E       *points ;
    ivector2D       *connectivity ;

    vector<sfield>   scalar_data; 
    vector<vfield>   vector_data;

    void Flush( fstream &str, string codex_, string name ) ;
    void Absorb( fstream &str, string codex_, string name ) ;
    
    public:
    VtkUnstrVec( ) ;
    VtkUnstrVec( string dir_, string name_, string codex_, int type_, dvecarr3E &points_ext, ivector2D &connectivity_external ) ;

   ~VtkUnstrVec( ) ;

    void    AddData( dvector1D &data, string name_, string loc_ ) ; 
    void    AddData( dvecarr3E &data, string name_, string loc_ ) ; 
   

};


#include"Class_VTK_Wrappers.tpp"

#endif
