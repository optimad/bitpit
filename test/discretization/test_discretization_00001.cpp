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

#include <array>
#if BITPIT_ENABLE_MPI==1
#   include <mpi.h>
#endif

#include "bitpit_IO.hpp"
#include "bitpit_discretization.hpp"
#include "bitpit_volcartesian.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing computation of reconstruction polynomial in cartesian 2D configuration
* using cell values only
*/
int subtest_001()
{
    int dimensions = 2;

    // Generate patch
    std::array<double, 3> patchOrigin = {{-2.5, -2.5, 0.}};
    double length = 5.;
    double h = 1;

    VolCartesian patch(dimensions, patchOrigin, length, h);
    patch.switchMemoryMode(VolCartesian::MEMORY_NORMAL);
    patch.write("subtest_001");

    // Test costant reconstruction
    {
        log::cout() << "Testing 2D cartesian constant" << std::endl;

        long cellId;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};
        std::vector<std::array<double, 3>> vertexCoords(ReferenceElementInfo::getInfo(ElementType::PIXEL).nVertices);

        Reconstruction reconstruction(0, dimensions);

        // Base cell
        cellId = 12;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());
    }

    // Test linear reconstruction
    {
        log::cout() << "Testing 2D cartesian linear, using a cross-type stencil " << std::endl;

        long cellId;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};
        std::vector<std::array<double, 3>> vertexCoords(ReferenceElementInfo::getInfo(ElementType::PIXEL).nVertices);

        Reconstruction reconstruction(1, dimensions);

        // Base cell
        cellId = 12;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // West neighbour
        cellId = 11;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // East neighbour
        cellId = 13;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // South neighbour
        cellId = 7;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // North neighbour
        cellId = 17;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());

        // Compute gradient stencil
        int stencilSize = reconstruction.getEquationCount();

        std::array<double, 3> point = {{0.,0.,0.}};

        StencilVector pointGradientStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointGradientStencil.appendItem(k, {{0.,0.,0.}});
        }
        reconstruction.computeGradientWeights(reconstructionOrigin, point, pointGradientStencil.weightData());
        pointGradientStencil.display(log::cout());
    }

    // Test linear reconstruction
    {
        log::cout() << "Testing 2D cartesian linear, using S, W and SW neighbours" << std::endl;

        long cellId;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};
        std::vector<std::array<double, 3>> vertexCoords(ReferenceElementInfo::getInfo(ElementType::PIXEL).nVertices);

        Reconstruction reconstruction(1, dimensions);

        // Base cell
        cellId = 12;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // West neighbour
        cellId = 11;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // South neighbour
        cellId = 7;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // South-west neighbour
        cellId = 6;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());
    }

    // Test quadratic reconstruction
    {
        log::cout() << "Testing 2D cartesian quadratic, using a block-type stencil " << std::endl;

        long cellId;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};
        std::vector<std::array<double, 3>> vertexCoords(ReferenceElementInfo::getInfo(ElementType::PIXEL).nVertices);

        Reconstruction reconstruction(2, dimensions);

        // Base cell
        cellId = 12;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // Neighbours
        for (int i = -1; i < 2; ++i) {
            for (int j = -1; j < 2; ++j) {
                if (i == 0 && j == 0) {
                    continue;
                }

                cellId = patch.getCellLinearId(2 + i, 2 + j, 0);
                patch.getCellVertexCoordinates(cellId, vertexCoords.data());
                reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
                reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());
            }
        }

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());
    }

    return 0;
}

/*!
* Subtest 002
*
* Testing computation of reconstruction polynomials in cartesian 2D configuration
* using cell values and derivatives
*/
int subtest_002()
{
    int dimensions = 2;

    // Generate patch
    std::array<double, 3> origin = {{-2.5, -2.5, 0.}};
    double length = 5.;
    double h = 1;

    VolCartesian patch(dimensions, origin, length, h);
    patch.switchMemoryMode(VolCartesian::MEMORY_NORMAL);
    patch.write("subtest_002");

    // Generate reconstruction
    {
        log::cout() << "Testing 2D cartesian linear, with cell values and derivatives" << std::endl;

        long cellId;
        std::array<double, 3> point;
        std::array<double, 3> direction;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};
        std::vector<std::array<double, 3>> vertexCoords(ReferenceElementInfo::getInfo(ElementType::PIXEL).nVertices);

        Reconstruction reconstruction(1, dimensions);

        // Base cell
        cellId = 12;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // West neighbour
        cellId = 11;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // South neighbour
        cellId = 7;
        patch.getCellVertexCoordinates(cellId, vertexCoords.data());
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addCellAverageEquation(reconstructionType, patch.getCell(cellId), reconstructionOrigin, vertexCoords.data());

        // Derivative equation
        point = {{0.5 * h, 0., 0.}};
        direction = {{1., 0., 0.}};
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointDerivativeEquation(reconstructionType, reconstructionOrigin, point, direction);

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());
    }

    return 0;
}

/*!
* Subtest 003
*
* Testing computation of reconstruction polynomials in 1D configuration
*/
int subtest_003()
{
    int dimensions = 1;

    double h = 1.;

    // Generate reconstruction
    {
        log::cout() << "Testing 1D stretched mesh, quadratic, with cell values " << std::endl;

        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};

        Reconstruction reconstruction(2, dimensions);
        int nCoeffs = reconstruction.getCoefficientCount();

        // Add equations
        std::array<double, 3> point_0 = {{0., 0., 0.}};
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_0);

        std::array<double, 3> point_1 = {{-0.75 * h, 0., 0.}};
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_1);

        std::array<double, 3> point_2 = {{1.5 * h, 0., 0.}};
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_2);

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());

        // Apply reconstruction
        std::vector<double> expectedCoeffs(nCoeffs);
        std::vector<double> pointValues(nCoeffs);
        std::array<double, 3> dist;

        expectedCoeffs[0] = 1.1;
        expectedCoeffs[1] = 2.2;
        expectedCoeffs[2] = 3.3;

        dist = point_0;
        pointValues[0] = expectedCoeffs[0] + expectedCoeffs[1] * dist[0] + 0.5  *  expectedCoeffs[2] * dist[0] * dist[0];

        dist = point_1;
        pointValues[1] = expectedCoeffs[0] + expectedCoeffs[1] * dist[0] + 0.5  *  expectedCoeffs[2] * dist[0] * dist[0];

        dist = point_2;
        pointValues[2] = expectedCoeffs[0] + expectedCoeffs[1] * dist[0] + 0.5  *  expectedCoeffs[2] * dist[0] * dist[0];

        ReconstructionPolynomial polynomial;
        reconstruction.assemblePolynomial(reconstructionOrigin, pointValues.data(), &polynomial);
        const double *evaluatedCoeffs = polynomial.getCoefficients();

        log::cout() << " Expected polynomial coefficients:  " << expectedCoeffs << std::endl;

        log::cout() << " Evaluated polynomial coefficients:";
        for (int k = 0; k < nCoeffs; ++k) {
            log::cout() << " " << evaluatedCoeffs[k];
        }
        log::cout() << std::endl;

        for (int k = 0; k < nCoeffs; ++k) {
            if (!utils::DoubleFloatingEqual()(expectedCoeffs[k], evaluatedCoeffs[k], 10)) {
                log::cout() << " Error in evaluation of polynomial coefficients" << std::endl;
                log::cout() << " Difference between expected and evaluated value is:  " << (expectedCoeffs[k] - evaluatedCoeffs[k]) << std::endl;
                return 1;
            }
        }
    }

    return 0;
}

/*!
* Subtest 004
*
* Testing computation of reconstruction polynomials in 1D configuration
*/
int subtest_004()
{
    int dimensions = 1;

    // Generate reconstruction stencil
    {
        log::cout() << "Testing 1D uniform mesh, linear, with cell values " << std::endl;

        double h = 1.;
        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};

        Reconstruction reconstruction(1, dimensions);

        // Add equations
        std::array<double, 3> point_0 = {{0., 0., 0.}};
        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_0);

        std::array<double, 3> point_1 = {{-h, 0., 0.}};
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_1);

        std::array<double, 3> point_2 = {{h, 0., 0.}};
        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_2);

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());

        // Apply reconstruction
        int stencilSize = reconstruction.getEquationCount();

        std::array<double, 3> point = {{0.5 * h, 0., 0.}};
        std::array<double, 3> direction = {{1., 0., 0.}};

        StencilScalar pointValueStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointValueStencil.appendItem(k, 0.);
        }
        reconstruction.computeValueWeights(reconstructionOrigin, point, pointValueStencil.weightData());
        pointValueStencil.display(log::cout());

        StencilScalar pointDerivativeStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointDerivativeStencil.appendItem(k, 0.);
        }

        reconstruction.computeDerivativeWeights(reconstructionOrigin, point, direction, pointDerivativeStencil.weightData());
        pointDerivativeStencil.display(log::cout());
    }

    return 0;
}

/*!
* Subtest 005
*
* Testing computation of reconstruction polynomials in 1D configuration
*/
int subtest_005()
{
    int dimensions = 2;

    // Generate reconstruction stencil
    {
        log::cout() << "Testing 2D uniform, under-determined, linear, with cell values " << std::endl;

        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};

        Reconstruction reconstruction(1, dimensions);

        // Add equations
        std::array<double, 3> point_0 = {{ 0., 0., 0.}};
        std::array<double, 3> point_1 = {{-1., 0., 0.}};
        std::array<double, 3> point_2 = {{ 1., 0., 0.}};

        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_0);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_1);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_2);

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());

        // Apply reconstruction
        int stencilSize = reconstruction.getEquationCount();

        std::array<double, 3> point = point_0;

        StencilScalar pointValueStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointValueStencil.appendItem(k, 0.);
        }
        reconstruction.computeValueWeights(reconstructionOrigin, point, pointValueStencil.weightData());
        pointValueStencil.display(log::cout());

        StencilVector pointGradientStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointGradientStencil.appendItem(k, {{0.,0.,0.}});
        }

        reconstruction.computeGradientWeights(reconstructionOrigin, point, pointGradientStencil.weightData());
        pointGradientStencil.display(log::cout());
    }

    // Generate reconstruction stencil
    {
        log::cout() << "Testing 2D uniform, under-determined, quadratic, with cell values " << std::endl;

        Reconstruction::ReconstructionType reconstructionType;
        std::array<double, 3> reconstructionOrigin = {{0., 0., 0.}};

        Reconstruction reconstruction(2, dimensions);

        // Add equations
        std::array<double, 3> point_0 = {{ 0., 0., 0.}};
        std::array<double, 3> point_1 = {{-1., 0., 0.}};
        std::array<double, 3> point_2 = {{ 1., 0., 0.}};
        std::array<double, 3> point_3 = {{ 0.,-1., 0.}};
        std::array<double, 3> point_4 = {{ 0., 1., 0.}};

        reconstructionType = Reconstruction::TYPE_CONSTRAINT;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_0);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_1);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_2);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_3);

        reconstructionType = Reconstruction::TYPE_LEAST_SQUARE;
        reconstruction.addPointValueEquation(reconstructionType, reconstructionOrigin, point_4);

        // Assembly reconstruction
        reconstruction.assemble();
        reconstruction.display(log::cout());

        // Apply reconstruction
        int stencilSize = reconstruction.getEquationCount();

        std::array<double, 3> point = point_0;

        StencilScalar pointValueStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointValueStencil.appendItem(k, 0.);
        }
        reconstruction.computeValueWeights(reconstructionOrigin, point, pointValueStencil.weightData());
        pointValueStencil.display(log::cout());

        StencilVector pointGradientStencil;
        for (int k = 0; k < stencilSize; ++k) {
            pointGradientStencil.appendItem(k, {{0.,0.,0.}});
        }

        reconstruction.computeGradientWeights(reconstructionOrigin, point, pointGradientStencil.weightData());
        pointGradientStencil.display(log::cout());
    }

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
    log::cout() << "Testing calculation of reconstruction stencils" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_003();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_004();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_005();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
