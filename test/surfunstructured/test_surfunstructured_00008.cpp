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

// ========================================================================== //
//           ** BitPit mesh ** Test 008 for class SurfUnstructured **         //
//                                                                            //
// Test multi-solid ASCI STL export                                           //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <iostream>

// BitPit
# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_surfunstructured.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// SUBTEST #001 Test multi-solid ASCI STL export ruotines                     //
// ========================================================================== //
int subtest_001(
    void
) {

    // Create the mesh
    SurfUnstructured *mesh = new SurfUnstructured(0);
    mesh->setExpert(true);

    std::vector<array<double,3>> verts(5);
    verts[0] = {{0,0,0}};
    verts[1] = {{1,0,0}};
    verts[2] = {{1,1,0}};
    verts[3] = {{0,1,0}};
    verts[4] = {{0.5,0.5,0}};

    std::vector<vector<long>> conns(4, vector<long>(3));
    conns[0] = {{3,0,4}};
    conns[1] = {{0,1,4}};
    conns[2] = {{1,2,4}};
    conns[3] = {{2,3,4}};

    for (int i = 0; i < 5; ++i) {
        mesh->addVertex(verts[i], i);
    }

    for (int j = 0; j < 4; ++j) {
        auto it = mesh->addCell(ElementInfo::TRIANGLE, true, conns[j], j);
        it->setPID(j);
    }

    bool isBinary = false;
    bool isMulti = true;
    mesh->exportSTL("multiPIDSquare.stl", isBinary, isMulti, false);

    return 0;
}

// ========================================================================== //
// MAIN FOR TEST #00006                                                       //
// ========================================================================== //
int main(
    void
) {

int                             err = 0;

err = subtest_001();
if (err > 0) return(10 + err);

return err;

}
