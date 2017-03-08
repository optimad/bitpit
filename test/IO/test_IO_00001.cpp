/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/


#include <iostream>

#include "VTK.hpp"

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

        bitpit::VTKUnstructuredGrid  vtk(".", "ustr1", bitpit::VTKElementType::VOXEL );

        vtk.setDimensions(1,8) ;
        vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
        vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;

        vtk.write() ;

    }

    { //Write grid and data to VTK in appended mode
        cout << "Write grid and data to VTK in appended mode" << endl;

        bitpit::VTKUnstructuredGrid  vtk(".", "ustr2", bitpit::VTKElementType::VOXEL );
        vtk.setDimensions(1,8) ;

        vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
        vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
        vtk.addData( "press", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, pressure) ;
        vtk.addData( "vel", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, velocity) ;

        vtk.write() ;
    }

    { //Read grid and data, rename and rexport
        cout << "Read grid and data, rename and rexport" << endl;

        vector<array<double,3>> Ipoints ;
        vector<vector<int>>     Iconnectivity ;

        vector<double>    Ipressure ;
        vector<array<double,3>>    Ivelocity ;


        bitpit::VTKUnstructuredGrid  vtk(".", "ustr2", bitpit::VTKElementType::VOXEL );

        vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
        vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
        vtk.addData( "press", Ipressure) ;

        vtk.read() ;

        vtk.removeData( "vel" ) ;

        Ipressure *= 2. ;

        vtk.setNames("./", "ustr3") ;
        vtk.write() ;
    }


    { //Read grid and data from Paraview-generated file, rename and rexport
        cout << "Read grid and data from Paraview-generated file, rename and rexport" << endl;

        vector< array<float,3>  >    Ipoints ;
        vector< vector<int64_t> >    Iconnectivity ;

        vector<float>   label ;
        vector<int64_t> cids, pids ;

        bitpit::VTKUnstructuredGrid  vtk("./data", "selection", bitpit::VTKElementType::TRIANGLE );
        vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
        vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
        vtk.addData( "STLSolidLabeling", label) ;
        vtk.addData( "vtkOriginalCellIds", cids) ;
        vtk.addData( "vtkOriginalPointIds", pids) ;

        vtk.read() ;

        vtk.setNames("./", "mySelection") ;
        vtk.setCounter( 0 ) ;
        vtk.write() ;
        vtk.write() ;
        vtk.disableData( "STLSolidLabeling" ) ;
        vtk.write( "mySelection_noLabel" ) ;
        vtk.enableData( "STLSolidLabeling" ) ;
        vtk.setName( "otherSelection") ;
        vtk.write(bitpit::VTKWriteMode::NO_SERIES) ;
        vtk.write(bitpit::VTKWriteMode::NO_INCREMENT) ;
    }

}
