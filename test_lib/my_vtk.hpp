#ifndef __CLASS_my_VTK__HH__
#define __CLASS_my_VTK__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include"Class_VTK.hpp"

using namespace std; 

typedef array<double,3>  vec1;
typedef array<int,8>     node_v;

class my_vtk_unstr : public VTK_UnstructuredGrid<my_vtk_unstr>{

  private:
    typedef  VTK_UnstructuredGrid<my_vtk_unstr> Base ;

    vector<vec1>    *points ;
    vector<node_v>  *connectivity ;

  public:
    my_vtk_unstr( ) ;
    my_vtk_unstr( string dir_, string name_, string codex_, int ncells_, int npoints_, int nconn_ ) ;
   ~my_vtk_unstr( ) ;

    void Link_Data( vector<vec1> &points_ext, vector<node_v> &connectivity_external ) ;

    void Flush( fstream &str, string codex_, string name ) ;
    void Absorb( fstream &str, string codex_, string name ) ;

    int Get_Connectivity_Size() ;

};

#include"my_vtk.tpp"

#endif
