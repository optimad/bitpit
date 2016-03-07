/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include "ParaTree.hpp"
#include "PabloUniform.hpp"

using namespace std;
using namespace bitpit;

// =================================================================================== //

void test001() {

	/**<Instantation and setup of a default (named bitpit) logfile.*/
	int errorFlag;
	int nproc;
	int	rank;
#if BITPIT_ENABLE_MPI==1
	MPI_Comm comm = MPI_COMM_WORLD;
	errorFlag = MPI_Comm_size(comm,&nproc);
	errorFlag = MPI_Comm_rank(comm,&rank);
#else
	nproc = 1;
	rank = 0;
#endif
	log::cout().setParallel(nproc, rank);
	log::cout() << fileVerbosity(log::NORMAL);
	log::cout() << consoleVerbosity(log::QUIET);

    /**<Instantation of a 2D para_tree object with default constructor.*/
    ParaTree ptreedefault;
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo001_default");


	/**<Instantation and setup of a custom (named custom) logfile.*/
	log::cout("custom").setParallel(nproc, rank);
	log::cout("custom") << fileVerbosity(log::NORMAL);
	log::cout("custom") << consoleVerbosity(log::QUIET);

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    int dim, maxlevel;
    dim = 2;
    maxlevel = 20;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    PabloUniform ptreecustom(X,Y,Z,L,dim,maxlevel,"custom");
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo001_custom");

    return ;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if BITPIT_ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        test001() ;

#if BITPIT_ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
