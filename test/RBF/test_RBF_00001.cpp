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


#include <iostream>

#include "bitpit_RBF.hpp"

using namespace std;

int main()
{

    bitpit::RBF     myRBF ;

	bool check;
	check =  (myRBF.evalBasis(0.9) > 0.0) && (myRBF.evalBasis(1.1) == 0);
	std::cout<<"value basis RBF in 0.9  "<<myRBF.evalBasis(0.9)<<std::endl;
	std::cout<<"value basis RBF in 1.1  "<<myRBF.evalBasis(1.1)<<std::endl;	
	int err = (int)(!check);
	return err;
}
