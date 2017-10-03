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
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>
# include <iostream>
#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// BitPit
# include "bitpit_common.hpp"
# include "bitpit_IO.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;

// ========================================================================== //
// GENERATE A TEST NON-MANIFOLD SURFACE TRIANGULATION FOR TESTS.              //
// ========================================================================== //
void generateTestTriangulation(
    SurfUnstructured                &mesh
) {

	// ========================================================================== //
	// void generateTestTriangulation(                                            //
	//     SurfUnstructured            &mesh)                                     //
	//                                                                            //
	// Generate a non-manifold surface triangulation for tests.                   //
	// ========================================================================== //
	// INPUT                                                                      //
	// ========================================================================== //
	// - mesh    : SurfUnstructured, surface mesh patch                           //
	// ========================================================================== //
	// OUTPUT                                                                     //
	// ========================================================================== //
	// - none                                                                     //
	// ========================================================================== //

	// ========================================================================== //
	// VARIABLES DECLARATION                                                      //
	// ========================================================================== //

	// Local variables
	long                    nV = 27;
	long                    nS = 33;

	// Counters
	int                     i;

	// ========================================================================== //
	// INITIALIZE TRIANGULATION                                                   //
	// ========================================================================== //
	{
		// Scope variables ------------------------------------------------------ //

		// Reserve memory for vertex & cell storage ----------------------------- //
		mesh.reserveVertices(nV);
		mesh.reserveCells(nS);
	}

	// ========================================================================== //
	// GENERATE VERTEX LIST                                                       //
	// ========================================================================== //
	{
		// Scope variables ------------------------------------------------------ //
		double                      off = -0.5;
		array<double, 3>            vertex;

		// 0-row ---------------------------------------------------------------- //
		vertex.fill(0.0);
		for (i = 0; i < 8; ++i) {
			vertex[0] = double(i);
			mesh.addVertex(vertex);
		} //next i

		// 1-row ---------------------------------------------------------------- //
		vertex.fill(0.0);
		vertex[1] = 1.0;
		for (i = 0; i < 9; ++i) {
			vertex[0] = double(i) + off;
			mesh.addVertex(vertex);
		} //next i

		// 2-row ---------------------------------------------------------------- //
		vertex.fill(0.0);
		vertex[1] = 2.0;
		for (i = 0; i < 8; ++i) {
			vertex[0] = double(i);
			mesh.addVertex(vertex);
		} //next i

		// Orthogonal element(s) ------------------------------------------------ //
		vertex[0] = 0.5*(mesh.getVertex(3)[0] + mesh.getVertex(12)[0]);
		vertex[1] = 0.5*(mesh.getVertex(3)[1] + mesh.getVertex(12)[1]);
		vertex[2] = 0.5 * sqrt(3.0);
		mesh.addVertex(vertex);
		vertex[0] = 0.5*(mesh.getVertex(12)[0] + mesh.getVertex(21)[0]);
		vertex[1] = 0.5*(mesh.getVertex(12)[1] + mesh.getVertex(21)[1]);
		vertex[2] = 0.5 * sqrt(3.0);
		mesh.addVertex(vertex);
	}

	// ========================================================================== //
	// GENERATE CONNECTIVITY                                                      //
	// ========================================================================== //
	{
		// Scope variables ------------------------------------------------------ //
		long int                    off;
		vector<long>                connectivity(3);

		// 0-row ---------------------------------------------------------------- //
		off = 8;
		for (i = 0; i < 7; ++i) {
			connectivity[0] = i;
			connectivity[1] = i + 1 + off;
			connectivity[2] = i + off;
			mesh.addCell(ElementType::TRIANGLE, true, connectivity);
			connectivity[0] = i;
			connectivity[1] = i + 1;
			connectivity[2] = i + 1 + off;
			mesh.addCell(ElementType::TRIANGLE, true, connectivity);
		} //next i
		connectivity[0] = i;
		connectivity[1] = i + 1 + off;
		connectivity[2] = i + off;
		mesh.addCell(ElementType::TRIANGLE, true, connectivity);

		// 1-row ---------------------------------------------------------------- //
		off = 9;
		for (i = 8; i < 15; ++i) {
			connectivity[0] = i;
			connectivity[1] = i + 1;
			connectivity[2] = i + off;
			mesh.addCell(ElementType::TRIANGLE, true, connectivity);
			connectivity[0] = i + 1;
			connectivity[1] = i + 1 + off;
			connectivity[2] = i + off;
			mesh.addCell(ElementType::TRIANGLE, true, connectivity);
		} //next i
		connectivity[0] = i;
		connectivity[1] = i + 1;
		connectivity[2] = i + off;
		mesh.addCell(ElementType::TRIANGLE, true, connectivity);

		// Orthogonal element --------------------------------------------------- //
		connectivity[0] = 3;
		connectivity[1] = 12;
		connectivity[2] = 25;
		mesh.addCell(ElementType::TRIANGLE, true, connectivity);
		connectivity[0] = 12;
		connectivity[1] = 26;
		connectivity[2] = 25;
		mesh.addCell(ElementType::TRIANGLE, true, connectivity);
		connectivity[0] = 12;
		connectivity[1] = 21;
		connectivity[2] = 26;
		mesh.addCell(ElementType::TRIANGLE, true, connectivity);
	}

}

// ========================================================================== //
// SUBTEST #001 Creating a mesh a renumber its ids.		                      //
//	returning 0 if successfull, 1 if failed.
// ========================================================================== //
int subtest_001()
{
	long counter;

	// Local variables
	SurfUnstructured                        mesh(2, 3);
	mesh.setExpert(true);

	// Generate a Dummy Triangulation
	generateTestTriangulation(mesh);
	mesh.buildAdjacencies();
	mesh.buildInterfaces();

	// Mess around with the numbering
	mesh.getVertices().updateId(0, 108);
	mesh.getVertices().updateId(1, 0);
	mesh.getVertex(0).setId(0);
	mesh.getVertices().updateId(108, 1);
	mesh.getVertex(1).setId(1);

	mesh.getCells().updateId(3, 108);
	mesh.getCells().updateId(4, 3);
	mesh.getCell(3).setId(3);
	mesh.getCells().updateId(108, 4);
	mesh.getCell(4).setId(4);

	// Get information on the triangulation
	size_t nVertices   = mesh.getVertexCount();
	size_t nCells      = mesh.getCellCount();
	size_t nInterfaces = mesh.getInterfaceCount();

	// Set the offsets
	long vertexOffset    = 0;
	long cellOffset      = 5;
	long interfaceOffset = 1012;

	// Creating expected id lists
	std::vector<long> expectedVertexIds(nVertices);
	counter = vertexOffset;
	for(auto & val : expectedVertexIds){
		val = counter;
		++counter;
	}

	std::vector<long> expectedCellIds(nCells);
	counter = cellOffset;
	for(auto & val : expectedCellIds){
		val = counter;
		++counter;
	}

	std::vector<long> expectedInterfaceIds(nInterfaces);
	counter = interfaceOffset;
	for(auto & val : expectedInterfaceIds){
		val = counter;
		++counter;
	}
	
	// Renumbering mesh
	mesh.consecutiveRenumber(vertexOffset,cellOffset,interfaceOffset);
	mesh.write("renumberedMesh");
	
	// Check renumbering of vertices
	if (mesh.getVertices().size() != nVertices) {
		return 1;
	}

	counter=0;
	for (const Vertex &vertex : mesh.getVertices()) {
		if (vertex.getId() != expectedVertexIds[counter]) {
			return 1;
		}
		++counter;
	}

	log::cout() << " - Vertices renumbered successfully" << std::endl;
	
	// Check renumbering of cells
	if (mesh.getCells().size() != nCells) {
		return 1;
	}

	counter=0;
	for (const Cell &cell : mesh.getCells()) {
		if (cell.getId() != expectedCellIds[counter]) {
			return 1;
		}

		++counter;
	}
	
	log::cout() << " - Cells renumbered successfully" << std::endl;

	// Check renumbering of interfaces
	if (mesh.getInterfaces().size() != nInterfaces) {
		return 1;
	}

	counter=0;
	for (const Interface &interface : mesh.getInterfaces()) {
		if (interface.getId() != expectedInterfaceIds[counter]) {
			return 1;
		}

		++counter;
	}

	log::cout() << " - Interfaces renumbered successfully" << std::endl;

	return 0;
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
        log::cout() << exception.what();
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
