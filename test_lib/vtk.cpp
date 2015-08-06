#include "Class_VTK_Wrappers.hpp"


void Demo_BasicVTK(){

  dvecarr3E    points, points2, points3 ;
  ivector2D    connectivity, connectivity2, connectivity3 ;

  dvector1D    pressure2, pressure3 ;
  dvecarr3E    velocity2, velocity3 ;


  points.resize(8) ;
  connectivity.resize(1) ;
  connectivity[0].resize(8) ;

  pressure2.resize(8) ;
  velocity2.resize(1) ;

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

  for( int i=0; i<8; i++)  pressure2[i] = (double) i ; 
  velocity2[0][0] = 1. ;
  velocity2[0][1] = 2. ;
  velocity2[0][2] = 3. ;


  { //Write only grid to VTK
  VtkUnstrVec   vtk(".", "ustr1", "appended", 11, points, connectivity );
  vtk.Write() ;
  }
 
  { //Read grid, add data and rexport
  VtkUnstrVec   vtk2(".", "ustr1", "appended", 11, points2, connectivity2 );
  vtk2.Read() ;
   
  vtk2.AddData( pressure2, "press", "Point") ;
  vtk2.AddData( velocity2, "vel", "Cell") ;
  
  vtk2.SetNames("./", "ustr2") ;
  vtk2.Write() ;
  }

  { //Read grid and data and rexport
  VtkUnstrVec   vtk3("./", "ustr2", "appended", 11, points3, connectivity3 );
  

  vtk3.ReadMetaData() ;
  
  vtk3.AddData( velocity3, "vel", "Cell") ;
  vtk3.AddData( pressure3, "press", "Point") ;

  vtk3.ReadData() ;

  pressure3 = 2.* pressure3 ;

  vtk3.SetNames("./", "ustr3") ;
  vtk3.Write() ;
  }



  return ;
};
