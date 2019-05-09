/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_IO.hpp"
#include "bitpit_RBF.hpp"

using namespace std;

/*!
* Subtest 001
*
* Testing basis evaluation.
*/
int subtest_001()
{
    bitpit::RBF     myRBF;

	bool check;
	check = (myRBF.evalBasis(0.9) > 0.0) && (myRBF.evalBasis(1.1) == 0);
	std::cout<<"value basis RBF in 0.9  "<<myRBF.evalBasis(0.9)<<std::endl;
	std::cout<<"value basis RBF in 1.1  "<<myRBF.evalBasis(1.1)<<std::endl;	
	int err = (int)(!check);

	return err;
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
