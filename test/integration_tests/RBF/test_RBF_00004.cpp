/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include <vector>
#include <array>
#include <algorithm>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_RBF.hpp"

using namespace bitpit;

//Description. Directly Interpolate a scalar (z-displacement) field

/*!
* Subtest 001
*
* Testing field interpolation.
*/
int subtest_001()
{
    // Variables
    std::vector< std::array<double,3> >         controlNodes;
    std::vector< std::array<double,3> >         points;
	std::vector< double >         				zDispl ;
	
    // fill control points
    int nCP = 4;
    controlNodes.resize(nCP);
    controlNodes[0] = std::array<double, 3>({0., 0., 0.});
    controlNodes[1] = std::array<double, 3>({0., 2., 0.});
    controlNodes[2] = std::array<double, 3>({2., 2., 0.});
    controlNodes[3] = std::array<double, 3>({2., 0., 0.});

    // fill displacements of control points
	zDispl.resize(nCP);
    zDispl[0] = 1.;
    zDispl[1] = 1.;
    zDispl[2] = 2.;
    zDispl[3] = 1.;

    
    // RBF used to directly interpolate zDispl field
    std::vector<double> supportRadii;
    supportRadii.resize(nCP);
    supportRadii[0] = 1.;
    supportRadii[1] = 0.75;
    supportRadii[2] = 1.25;
    supportRadii[3] = 0.75;

    bitpit::RBF paraMorph;
    RBFBasisFunction funct = RBFBasisFunction::WENDLANDC2;
	paraMorph.setFunction(funct);
    paraMorph.setSupportRadius(supportRadii);
    paraMorph.enablePolynomial();

    std::vector<int> nIndex = paraMorph.addNode(controlNodes);
	paraMorph.addData(zDispl);
    int err = paraMorph.solve() ;
	if(err > 0 ) 	return 1;
	
    // create moving points
    points.resize(nCP+1);
    points[0] = std::array<double, 3>({0., 0., 0.});
    points[1] = std::array<double, 3>({0., 2., 0.});
    points[2] = std::array<double, 3>({2., 2., 0.});
    points[3] = std::array<double, 3>({2., 0., 0.});
    points[4] = std::array<double, 3>({1., 1., 0.});

    std::vector<double> disp;
    for( auto & point : points){
        std::cout << point << " -----> ";
        disp = paraMorph.evalRBF(point);
        point[2] += disp[0] ;
        std::cout << point << std::endl;
    };

    double val = points[4][2];
    bool check = !bitpit::utils::DoubleFloatingEqual()(val, 1.25);


    // Use a different function
    paraMorph.removeAllData();
    funct = RBFBasisFunction::THINPLATE;
    paraMorph.setFunction(funct);
    paraMorph.setSupportRadius(1.);
    paraMorph.enablePolynomial();

    paraMorph.addData(zDispl);
    err = paraMorph.solve();
    if (err > 0)
        return 1;

    // reset moving points
    points[0] = std::array<double, 3>({0., 0., 0.});
    points[1] = std::array<double, 3>({0., 2., 0.});
    points[2] = std::array<double, 3>({2., 2., 0.});
    points[3] = std::array<double, 3>({2., 0., 0.});
    points[4] = std::array<double, 3>({1., 1., 0.});

    // move points
    for (auto &point : points) {
        std::cout << point << " -----> ";
        disp = paraMorph.evalRBF(point);
        point[2] += disp[0];
        std::cout << point << std::endl;
    };

    val = points[4][2];
    check = check || !bitpit::utils::DoubleFloatingEqual()(val, 1.25);
    return check;
}

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        bitpit::log::cout() << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
