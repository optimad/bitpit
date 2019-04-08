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

#include <array>

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_operators.hpp"
#include "bitpit_volunstructured.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing basic features of a 2D patch.
*/
int subtest_001()
{
    std::array<double, 3> minPoint;
    std::array<double, 3> maxPoint;

    log::manager().initialize(log::COMBINED);
    log::cout() << "Testing unstructured patches\n";

    log::cout() << "\n\n:: 2D unstructured patch ::\n";

    VolUnstructured *patch_2D = new VolUnstructured(0, 2);
    patch_2D->getVTK().setName("unstructured_uniform_patch_2D");

    patch_2D->addVertex({{0.00000000, 0.00000000, 0.00000000}},  1);
    patch_2D->addVertex({{0.00000000, 1.00000000, 0.00000000}},  2);
    patch_2D->addVertex({{1.00000000, 1.00000000, 0.00000000}},  3);
    patch_2D->addVertex({{1.00000000, 0.00000000, 0.00000000}},  4);
    patch_2D->addVertex({{1.00000000, 0.50000000, 0.00000000}},  5);
    patch_2D->addVertex({{0.25992107, 1.00000000, 0.00000000}},  6);
    patch_2D->addVertex({{0.58740113, 1.00000000, 0.00000000}},  7);
    patch_2D->addVertex({{0.00000000, 0.75000000, 0.00000000}},  8);
    patch_2D->addVertex({{0.00000000, 0.50000000, 0.00000000}},  9);
    patch_2D->addVertex({{0.00000000, 0.25000000, 0.00000000}}, 10);
    patch_2D->addVertex({{0.25992107, 0.00000000, 0.00000000}}, 11);
    patch_2D->addVertex({{0.58740113, 0.00000000, 0.00000000}}, 12);
    patch_2D->addVertex({{0.42807699, 0.41426491, 0.00000000}}, 13);
    patch_2D->addVertex({{0.30507278, 0.69963441, 0.00000000}}, 14);
    patch_2D->addVertex({{0.64032722, 0.68239464, 0.00000000}}, 15);
    patch_2D->addVertex({{0.24229808, 0.24179558, 0.00000000}}, 16);
    patch_2D->addVertex({{0.67991107, 0.28835559, 0.00000000}}, 17);
    patch_2D->addVertex({{0.22034760, 0.48841527, 0.00000000}}, 18);
    patch_2D->addVertex({{0.43952167, 0.18888322, 0.00000000}}, 19);

    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 3,  7, 15}}));
    patch_2D->addCell(ElementType::QUAD,     std::vector<long>({{ 1, 11, 16, 10}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 8,  9, 18}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 8, 18, 14}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 3, 15,  5}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 9, 10, 18}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{10, 16, 18}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 4, 17, 12}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 4,  5, 17}}));
    patch_2D->addCell(ElementType::QUAD,     std::vector<long>({{13, 17, 15, 14}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{11, 12, 19}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 19, 17}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{12, 17, 19}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 14, 18}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{ 5, 15, 17}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{11, 19, 16}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 16, 19}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{13, 18, 16}}));
    patch_2D->addCell(ElementType::TRIANGLE, std::vector<long>({{14, 15,  7}}));
    patch_2D->addCell(ElementType::POLYGON,  std::vector<long>({{ 5,  2,  8, 14, 7, 6}}));

    patch_2D->buildAdjacencies();
    patch_2D->buildInterfaces();

    patch_2D->update();
    patch_2D->write();

    log::cout() << std::endl;
    log::cout() << ">> 2D volume\n";
    log::cout() << std::endl;

    double volume_2D = 0.;
    for (const Cell &cell : patch_2D->getCells()) {
        volume_2D += patch_2D->evalCellVolume(cell.getId());
    }

    log::cout() << " - Total volume : " << volume_2D << std::endl;

    double volume_expected_2D = 1.0;
    if (std::abs(volume_2D - volume_expected_2D) > 1e-12) {
        throw std::runtime_error("Volume of the 2D patch doesn't match the expected value");
    }

    log::cout() << std::endl;
    log::cout() << ">> 2D surface area\n";
    log::cout() << std::endl;

    double surfaceArea_2D = 0.;
    for (const Interface &interface : patch_2D->getInterfaces()) {
        if (!interface.isBorder()) {
           continue;
        }

        surfaceArea_2D += patch_2D->evalInterfaceArea(interface.getId());
    }

    log::cout() << " - Surface area : " << surfaceArea_2D << std::endl;

    double surfaceArea_expected_2D = 4.0;
    if (std::abs(surfaceArea_2D - surfaceArea_expected_2D) > 1e-12) {
        throw std::runtime_error("Surface area of the 2D patch doesn't match the expected value");
    }

    log::cout() << std::endl;
    log::cout() << ">> 2D bounding box\n";
    log::cout() << std::endl;

    patch_2D->getBoundingBox(minPoint, maxPoint);

    log::cout() << "  x : (" << minPoint[0] << ", " << maxPoint[0] << ")\n";
    log::cout() << "  y : (" << minPoint[1] << ", " << maxPoint[1] << ")\n";
    log::cout() << "  z : (" << minPoint[2] << ", " << maxPoint[2] << ")\n";

    log::cout() << std::endl;
    log::cout() << ">> 2D neighbour test\n";

    std::vector<long> neighs_2D;

    long cellId_2D = 19;
    log::cout() << std::endl;
    log::cout() << "Cell id: " << cellId_2D << std::endl << std::endl;

    log::cout() << "Face neighbours (complete list):\n";
    neighs_2D = patch_2D->findCellFaceNeighs(cellId_2D);
    for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
        log::cout() << " - " << neighs_2D[i] << std::endl;
    }

    log::cout() << "\nVertex neighbours (complete list):\n";
    neighs_2D = patch_2D->findCellVertexNeighs(cellId_2D, true);
    for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
        log::cout() << " - " << neighs_2D[i] << std::endl;
    }

    log::cout() << "\nVertex neighbours (excuding face neighbours):\n";
    neighs_2D = patch_2D->findCellVertexNeighs(cellId_2D, false);
    for (unsigned int i = 0; i < neighs_2D.size(); ++i) {
        log::cout() << " - " << neighs_2D[i] << std::endl;
    }

    log::cout() << std::endl;

    delete patch_2D;

    return 0;
}

/*!
* Subtest 002
*
* Testing basic features of a 3D patch.
*/
int subtest_002()
{
    std::array<double, 3> minPoint;
    std::array<double, 3> maxPoint;

    log::cout() << "\n\n:: 3D unstructured mesh ::\n";

    VolUnstructured *patch_3D = new VolUnstructured(0, 3);
    patch_3D->getVTK().setName("unstructured_uniform_patch_3D");

    patch_3D->addVertex({{0.00000000, 0.00000000,  0.00000000}},  1);
    patch_3D->addVertex({{1.00000000, 0.00000000,  0.00000000}},  2);
    patch_3D->addVertex({{1.00000000, 1.00000000,  0.00000000}},  3);
    patch_3D->addVertex({{0.00000000, 1.00000000, -0.75000000}},  4);
    patch_3D->addVertex({{0.00000000, 0.00000000,  1.00000000}},  5);
    patch_3D->addVertex({{1.00000000, 0.00000000,  1.00000000}},  6);
    patch_3D->addVertex({{1.00000000, 1.00000000,  1.00000000}},  7);
    patch_3D->addVertex({{0.00000000, 1.00000000,  1.00000000}},  8);
    patch_3D->addVertex({{1.00000000, 0.00000000,  0.54678323}},  9);
    patch_3D->addVertex({{0.50000000, 0.00000000,  1.00000000}}, 10);
    patch_3D->addVertex({{0.00000000, 0.00000000,  0.54678323}}, 11);
    patch_3D->addVertex({{1.00000000, 1.00000000,  0.54678323}}, 12);
    patch_3D->addVertex({{0.50000000, 1.00000000,  1.00000000}}, 13);
    patch_3D->addVertex({{0.00000000, 1.00000000,  0.54678323}}, 14);
    patch_3D->addVertex({{1.00000000, 0.50000000,  1.00000000}}, 15);
    patch_3D->addVertex({{0.00000000, 0.50000000,  1.00000000}}, 16);
    patch_3D->addVertex({{0.51053620, 0.00000000,  0.34680184}}, 17);
    patch_3D->addVertex({{0.36278402, 0.00000000,  0.68603230}}, 18);
    patch_3D->addVertex({{0.69618860, 0.00000000,  0.73234294}}, 19);
    patch_3D->addVertex({{0.51053620, 1.00000000,  0.34680184}}, 20);
    patch_3D->addVertex({{0.36278402, 1.00000000,  0.68603230}}, 21);
    patch_3D->addVertex({{0.69618860, 1.00000000,  0.73234294}}, 22);
    patch_3D->addVertex({{1.00000000, 0.51053620,  0.34680184}}, 23);
    patch_3D->addVertex({{1.00000000, 0.36278402,  0.68603230}}, 24);
    patch_3D->addVertex({{1.00000000, 0.69618860,  0.73234294}}, 25);
    patch_3D->addVertex({{0.00000000, 0.51053620,  0.34680184}}, 26);
    patch_3D->addVertex({{0.00000000, 0.36278402,  0.68603230}}, 27);
    patch_3D->addVertex({{0.00000000, 0.69618860,  0.73234294}}, 28);
    patch_3D->addVertex({{0.50000000, 0.50000000,  1.00000000}}, 29);
    patch_3D->addVertex({{0.75000000, 0.25000000,  1.00000000}}, 30);
    patch_3D->addVertex({{0.25000000, 0.25000000,  1.00000000}}, 31);
    patch_3D->addVertex({{0.00000000, 0.00000000, -1.50000000}}, 32);
    patch_3D->addVertex({{1.00000000, 0.00000000, -1.00000000}}, 33);
    patch_3D->addVertex({{1.00000000, 1.00000000, -0.50000000}}, 34);
    patch_3D->addVertex({{0.00000000, 1.00000000, -1.00000000}}, 35);
    patch_3D->addVertex({{0.00000000, 0.00000000, -2.25000000}}, 36);
    patch_3D->addVertex({{1.00000000, 0.00000000, -1.66666666}}, 37);
    patch_3D->addVertex({{1.00000000, 1.00000000, -2.12500000}}, 38);
    patch_3D->addVertex({{0.00000000, 1.00000000, -2.00000000}}, 39);
    patch_3D->addVertex({{0.00000000, 0.00000000, -3.00000000}}, 40);
    patch_3D->addVertex({{1.00000000, 0.00000000, -3.00000000}}, 41);
    patch_3D->addVertex({{1.00000000, 1.00000000, -3.00000000}}, 42);
    patch_3D->addVertex({{0.00000000, 1.00000000, -3.00000000}}, 43);
    patch_3D->addVertex({{0.00000000, 0.00000000, -4.00000000}}, 44);
    patch_3D->addVertex({{1.00000000, 0.00000000, -4.00000000}}, 45);
    patch_3D->addVertex({{1.00000000, 1.00000000, -4.00000000}}, 46);
    patch_3D->addVertex({{0.00000000, 1.00000000, -4.00000000}}, 47);
    patch_3D->addVertex({{0.50000000, 0.00000000, -4.00000000}}, 48);
    patch_3D->addVertex({{1.00000000, 0.50000000, -4.00000000}}, 49);
    patch_3D->addVertex({{0.50000000, 1.00000000, -4.00000000}}, 50);
    patch_3D->addVertex({{0.33333333, 0.66666666, -4.00000000}}, 51);
    patch_3D->addVertex({{0.00000000, 0.50000000, -4.00000000}}, 52);

    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 22, 25, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 27, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 28, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26, 23, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 25, 23}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 21, 22, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 28, 27, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23, 26, 29, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 23, 29, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 29, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 29, 23, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 29, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{29, 25, 23, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 23,  3, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{31, 18, 27, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{22, 12, 25, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 26, 11, 27}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 19, 24, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 21, 14, 28}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 22,  7, 25}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 29, 22, 25}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{25, 12, 23, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23,  4,  3,  2}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26, 18, 11, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 9, 24, 23, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21, 26, 14, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{28, 21, 13, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 18, 19, 10}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 29, 18}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 31, 18, 10}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21,  8, 28, 13}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 29, 28}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 13, 28,  8}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 22,  7}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{13, 15, 29, 22}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{15, 24, 29, 25}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{21, 13, 29, 22}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{28, 16, 27, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{31, 18,  5, 27}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 30, 15, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 31, 27, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 11,  5, 27}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 19, 24,  6}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 1, 26, 11, 17}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{19,  9, 24,  6}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{14, 21,  8, 28}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 31,  5, 10}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23, 12,  3, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{19, 30, 10,  6}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 7, 22, 12, 25}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 30,  6, 15}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{16, 31,  5, 27}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{23,  9, 17,  2}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{26,  4, 14, 20}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 17, 24, 19}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{24, 17, 18, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{18, 24, 30, 19}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{30, 24, 18, 29}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 17, 26,  1}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 2, 26, 17, 23}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4,  2, 26,  1}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{ 4, 26,  2, 23}}));
    patch_3D->addCell(ElementType::WEDGE,      std::vector<long>({{ 2, 4, 1, 33, 35, 32}}));
    patch_3D->addCell(ElementType::WEDGE,      std::vector<long>({{ 4, 2, 3, 35, 33, 34}}));
    patch_3D->addCell(ElementType::PYRAMID,    std::vector<long>({{36, 37, 38, 39, 33}}));
    patch_3D->addCell(ElementType::PYRAMID,    std::vector<long>({{39, 38, 34, 35, 33}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{35, 36, 39, 33}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{32, 35, 33, 36}}));
    patch_3D->addCell(ElementType::HEXAHEDRON, std::vector<long>({{42, 43, 40, 41, 38, 39, 36, 37}}));
    patch_3D->addCell(ElementType::POLYHEDRON, std::vector<long>({{11,
                                                                         4, 42, 43, 40, 41,
                                                                         5, 52, 51, 50, 49, 48,
                                                                         3, 41, 40, 48,
                                                                         3, 42, 41, 49,
                                                                         3, 43, 42, 50,
                                                                         3, 40, 43, 52,
                                                                         3, 43, 51, 52,
                                                                         3, 43, 50, 51,
                                                                         3, 42, 49, 50,
                                                                         3, 41, 48, 49,
                                                                         3, 40, 52, 48
                                                                       }}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{52, 51, 47, 43}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{51, 50, 47, 43}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{44, 48, 52, 40}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{48, 45, 49, 41}}));
    patch_3D->addCell(ElementType::TETRA,      std::vector<long>({{50, 49, 46, 42}}));

    patch_3D->buildAdjacencies();
    patch_3D->buildInterfaces();

    patch_3D->update();
    patch_3D->write();

    log::cout() << std::endl;
    log::cout() << ">> 3D volume\n";
    log::cout() << std::endl;

    double volume_3D = 0.;
    for (const Cell &cell : patch_3D->getCells()) {
        volume_3D += patch_3D->evalCellVolume(cell.getId());
    }

    log::cout() << " - Total volume : " << volume_3D << std::endl;

    double volume_expected_3D = 5.0;
    if (std::abs(volume_3D - volume_expected_3D) > 1e-12) {
        throw std::runtime_error("Volume of the 3D patch doesn't match the expected value");
    }

    log::cout() << std::endl;
    log::cout() << ">> 3D surface area\n";
    log::cout() << std::endl;

    double surfaceArea_3D = 0.;
    for (const Interface &interface : patch_3D->getInterfaces()) {
        if (!interface.isBorder()) {
           continue;
        }

        surfaceArea_3D += patch_3D->evalInterfaceArea(interface.getId());
    }

    log::cout() << "  Surface area : " << surfaceArea_3D << std::endl;

    double surfaceArea_expected_3D = 22.0;
    if (std::abs(surfaceArea_3D - surfaceArea_expected_3D) > 1e-12) {
        throw std::runtime_error("Surface area of the 3D patch doesn't match the expected value");
    }

    log::cout() << std::endl;
    log::cout() << ">> 3D divergence\n";
    log::cout() << std::endl;

    std::array<double, 3> field = {{1., 2., 3.}};

    double divergence_3D = 0.;
    for (const Cell &cell : patch_3D->getCells()) {
        long cellId = cell.getId();
        int nCellInterfaces = cell.getInterfaceCount();
        const long *interfaces = cell.getInterfaces();

        double cellDivergence_3D = 0.;
        for (int k = 0; k < nCellInterfaces; ++k) {
            long interfaceId = interfaces[k];
            const Interface &interface = patch_3D->getInterface(interfaceId);

            double area = patch_3D->evalInterfaceArea(interfaceId);
            std::array<double, 3> normal = patch_3D->evalInterfaceNormal(interfaceId);

            int sign = 1;
            if (cellId != interface.getOwner()) {
                sign *= -1;
            }

            cellDivergence_3D += sign * area * dotProduct(field, normal);
        }
        cellDivergence_3D /= patch_3D->evalCellVolume(cell.getId());

        divergence_3D += std::abs(cellDivergence_3D);
    }

    log::cout() << "  Field      : " << field << std::endl;
    log::cout() << "  Divergence : " << divergence_3D << std::endl;

    double divergence_expected_3D = 0.0;
    if (std::abs(divergence_3D - divergence_expected_3D) > 1e-12) {
        throw std::runtime_error("Divergence of the 3D patch doesn't match the expected value");
    }

    log::cout() << std::endl;
    log::cout() << ">> 3D bounding box\n";
    log::cout() << std::endl;

    patch_3D->getBoundingBox(minPoint, maxPoint);

    log::cout() << "  x : (" << minPoint[0] << ", " << maxPoint[0] << ")\n";
    log::cout() << "  y : (" << minPoint[1] << ", " << maxPoint[1] << ")\n";
    log::cout() << "  z : (" << minPoint[2] << ", " << maxPoint[2] << ")\n";

    log::cout() << std::endl;
    log::cout() << ">> 3D neighbour test\n";

    std::vector<long> neighs_3D;

    long cellId_3D = 69;
    log::cout() << std::endl;
    log::cout() << "Cell id: " << cellId_3D << std::endl << std::endl;

    log::cout() << "Face neighbours (complete list):\n";
    neighs_3D = patch_3D->findCellFaceNeighs(cellId_3D);
    for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
        log::cout() << " - " << neighs_3D[i] << std::endl;
    }

    log::cout() << "\nEdge neighbours (complete list):\n";
    neighs_3D = patch_3D->findCellEdgeNeighs(cellId_3D, true);
    for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
        log::cout() << " - " << neighs_3D[i] << std::endl;
    }

    log::cout() << "\nEdge neighbours (excuding face neighbours):\n";
    neighs_3D = patch_3D->findCellEdgeNeighs(cellId_3D, false);
    for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
        log::cout() << " - " << neighs_3D[i] << std::endl;
    }

    log::cout() << "\nVertex neighbours (complete list):\n";
    neighs_3D = patch_3D->findCellVertexNeighs(cellId_3D, true);
    for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
        log::cout() << " - " << neighs_3D[i] << std::endl;
    }

    log::cout() << "\nVertex neighbours (excuding face and edge neighbours):\n";
    neighs_3D = patch_3D->findCellVertexNeighs(cellId_3D, false);
    for (unsigned int i = 0; i < neighs_3D.size(); ++i) {
        log::cout() << " - " << neighs_3D[i] << std::endl;
    }

    log::cout() << std::endl;

    delete patch_3D;

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Testing basic features of volunstructured patches" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }

        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
