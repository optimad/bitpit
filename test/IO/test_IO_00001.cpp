
#include <iostream>

#include "VTKWrappers.hpp"

using namespace std;

int main()
{

    vector<array<double,3>>     points ;
    vector<vector<int>>          connectivity ;

    vector<double>              pressure ;
    vector<array<double,3>>     velocity ;


    points.resize(8) ;
    connectivity.resize(1) ;
    connectivity[0].resize(8) ;

    pressure.resize(8) ;
    velocity.resize(1) ;

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

    for( int i=0; i<8; i++)  pressure[i] = (double) i ;
    velocity[0][0] = 1. ;
    velocity[0][1] = 2. ;
    velocity[0][2] = 3. ;

    { //Write only grid to VTK in ascii format
        cout << "Write only grid to VTK in ascii format" << endl;

        VTKUnstructuredVec   vtk(".", "ustr1", VTKFormat::ASCII, VTKElementType::VOXEL, points, connectivity );
        vtk.write() ;

    }

    { //Write grid and data to VTK in appended mode
        cout << "Write grid and data to VTK in appended mode" << endl;

        VTKUnstructuredVec   vtk(".", "ustr2", VTKFormat::APPENDED, VTKElementType::VOXEL, points, connectivity );
        vtk.addData( pressure, "press", VTKLocation::POINT ) ;
        vtk.addData( velocity, "vel", VTKLocation::CELL ) ;
        vtk.write() ;
    }

    { //Read grid and data, rename and rexport
        cout << "Read grid and data, rename and rexport" << endl;

        vector<array<double,3>> Ipoints ;
        vector<vector<int>>     Iconnectivity ;

        vector<double>    Ipressure ;
        vector<array<double,3>>    Ivelocity ;


        VTKUnstructuredVec   vtk(".", "ustr2", VTKFormat::APPENDED, VTKElementType::VOXEL );

        vtk.linkData( Ipoints, "Points") ;
        vtk.linkData( Iconnectivity, "connectivity") ;
        vtk.linkData( Ipressure, "press") ;
        vtk.linkData( Ivelocity, "vel") ;

        vtk.read() ;

        Ipressure = 2. * Ipressure ;

        vtk.setNames("./", "ustr3") ;
        vtk.write() ;
    }


    { //Read grid and data from Paraview-generated file, rename and rexport
        cout << "Read grid and data from Paraview-generated file, rename and rexport" << endl;

        vector< array<float,3>  >    Ipoints ;
        vector< vector<int64_t> >    Iconnectivity ;

        vector<float>   label ;
        vector<int64_t> ids ;

        VTKUnstructuredVec   vtk("./data", "selection", VTKFormat::APPENDED, VTKElementType::TRIANGLE );

        vtk.linkData( Ipoints, "Points") ;
        vtk.linkData( Iconnectivity, "connectivity") ;
        vtk.linkData( label, "STLSolidLabeling") ;
        vtk.linkData( ids, "vtkOriginalCellIds") ;

        vtk.read() ;

        vtk.setNames("./", "my_selection") ;
        vtk.write() ;
    }

}
