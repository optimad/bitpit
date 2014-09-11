#include <cmath>
#include <vector>
#include <array>
#include <iostream>

//#include "Class_VTK.hpp"
#include "my_vtk.hpp"
#include "STL_IOFunct.hpp"
#include "DGF_IOFunct.hpp"


typedef array<double,3>  vec1;
typedef array<int,8>     node_v;

int main(){

  vector<vec1>    points ;
  vector<node_v>  connectivity ;


  points.resize(8) ;
  connectivity.resize(1) ;

  points[0][0]  =  0.  ;
  points[0][1]  =  0.  ;
  points[0][2]  =  0.  ;

  points[1][0]  =  1.  ;
  points[1][1]  =  0.  ;
  points[1][2]  =  0.  ;

  points[2][0]  =  0.  ;
  points[2][1]  =  1.  ;
  points[2][2]  =  0.  ;

  points[3][0]  =  1.  ;
  points[3][1]  =  1.  ;
  points[3][2]  =  0.  ;

  points[4][0]  =  0.  ;
  points[4][1]  =  0.  ;
  points[4][2]  =  1.  ;

  points[5][0]  =  1.  ;
  points[5][1]  =  0.  ;
  points[5][2]  =  1.  ;

  points[6][0]  =  0.  ;
  points[6][1]  =  1.  ;
  points[6][2]  =  1.  ;

  points[7][0]  =  1.  ;
  points[7][1]  =  1.  ;
  points[7][2]  =  1.  ;


  for( int i=0; i<8; i++)  connectivity[0][i] = i ; 



  {  //TESTING VTK

    int     ncells, nvertices, nconnectivity, k ;

    ncells = 1;
    nvertices = 8 ;
    nconnectivity = 8 ;
 
    my_vtk_unstr   vtk(".", "ustr1", "appended", ncells, nvertices, nconnectivity );

    vtk.Set_Parallel(7,0) ;
    vtk.Set_Counter(3) ;

    vtk.Link_Data( points, connectivity ) ;

    vtk.Write() ;
    vtk.Write() ;
    vtk.Write() ;
    vtk.Write() ;

//    my_vtk_unstr   vtk2;
//   
//    vtk2.Set_Parallel(7,0) ;
//    vtk2.Set_Names( ".", "ustr1" ) ;
//    vtk2.Link_Data( points, connectivity ) ;
//    vtk2.Read() ;

  };


  {  //TESTING STL_IO
  
  };


  {  //TESTING DGF_IO

  };

  return 0;
};
