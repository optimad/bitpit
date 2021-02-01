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
#include <mpi.h>
#endif

#include "bitpit_IO.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/*!
	\example voloctree_adaptation_example_00001.cpp

	\brief 3D mesh adaptation using voloctree

	This example creates a 3D octree mesh on the square domain [0,1]x[0,1].
	On this domain two fields, a scalar and a vector, are defined.

	Mesh is first refined and then coarsened, during each mesh adaption step
	fields are remapped on the updated mesh.

	<b>To run</b>: ./voloctree_adaptation_example_00001.cpp \n
*/

// Auxiliary class to write fields in VTK format.
class FieldStreamer : public VTKBaseStreamer {

public:
    FieldStreamer(const PatchKernel &patch, const PiercedStorage<double, long> &scalarField, const PiercedStorage<std::array<double, 3>, long> &vectorField)
        : m_patch(patch), m_scalarField(scalarField), m_vectorField(vectorField)
    {
    };

    void flushData(std::fstream &stream, const std::string &name, VTKFormat format)
    {
        assert(format == VTKFormat::APPENDED);
        BITPIT_UNUSED(format);

        if (name == "scalarField") {
            for (const Cell &cell : m_patch.getVTKCellWriteRange()) {
                long id = cell.getId();
                genericIO::flushBINARY(stream, m_scalarField.at(id));
            }
        } else if (name == "vectorField") {
            for (const Cell &cell : m_patch.getVTKCellWriteRange()) {
                long id = cell.getId();
                genericIO::flushBINARY(stream, m_vectorField.at(id));
            }
        }
    };

private:
    const PatchKernel &m_patch;
    const PiercedStorage<double, long> &m_scalarField;
    const PiercedStorage<std::array<double, 3>, long> &m_vectorField;

};

// Main program
int main(int argc, char *argv[])
{
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);


#if ENABLE_MPI == 1
    // Initialize MPI
    MPI_Init(&argc, &argv);
#endif

    //
    // Domain initialization
    //

    // Initializing mesh
    log::cout() << " Initializing mesh..." << std::endl;

    double length = 1.0;
    std::array<double, 3> center  = {{0.0,0.0,0.0}};
    std::array<double, 3> minimum = center  - length / 2.;
    double dh = length / 16.;

#if BITPIT_ENABLE_MPI
    VolOctree mesh(3, minimum, length, dh, MPI_COMM_NULL);
#else
    VolOctree mesh(3, minimum, length, dh);
#endif
    mesh.initializeAdjacencies();
    mesh.update();

    // Initialize fields
    log::cout() << "Initializing fields..." << std::endl;

    PiercedStorage<double, long> scalarField;
    PiercedStorage<std::array<double, 3>, long> vectorField;
    scalarField.setDynamicKernel(&(mesh.getCells()), PiercedVector<Cell>::SYNC_MODE_JOURNALED);
    vectorField.setDynamicKernel(&(mesh.getCells()), PiercedVector<Cell>::SYNC_MODE_JOURNALED);
    for (const Cell &cell : mesh.getCells()) {
        const long cellId = cell.getId();
        std::array<double, 3> cellCentroid = mesh.evalCellCentroid(cellId);
        double r = std::sqrt(cellCentroid[0] * cellCentroid[0] + cellCentroid[1] * cellCentroid[1] + cellCentroid[2] * cellCentroid[2]);

        scalarField.set(cellId, r);
        vectorField.set(cellId, {{cellCentroid[0] / r, cellCentroid[1] / r, cellCentroid[2] / r}});
    }

    // Initialize mesh output
    FieldStreamer dataStreamer(mesh, scalarField, vectorField);

    mesh.getVTK().setName("voloctree_adaptation_example_00001");
    mesh.getVTK().setCounter();
    mesh.getVTK().addData<double>("scalarField", VTKFieldType::SCALAR, VTKLocation::CELL, &dataStreamer);
    mesh.getVTK().addData<std::array<double, 3>>("vectorField", VTKFieldType::VECTOR, VTKLocation::CELL, &dataStreamer);

    // Write mesh
    mesh.write();

    //
    // Mesh adaptation
    //
    bool trackAdaptation = true;
    bool squeeshPatchStorage = false;
    std::vector<adaption::Info> adaptionData;

    // Refinement
    log::cout() << "Performing refinement..." << std::endl;

    int nRefinements = 2;
    for (int i = 0; i < nRefinements; ++i) {
        // Mark cells that need refinement
        //
        // All cells inside a circle centered in the origin and with a radius
        // equal to 0.2 will be refined.
        log::cout() << "  Marking cells that need refinement..." << std::endl;

        for (const Cell &cell : mesh.getCells()) {
            long cellId = cell.getId();
            std::array<double, 3> cellCentroid = mesh.evalCellCentroid(cellId);

            double r = std::sqrt(cellCentroid[0] * cellCentroid[0] + cellCentroid[1] * cellCentroid[1] + cellCentroid[2] * cellCentroid[2]);
            if (r < 0.2) {
                mesh.markCellForRefinement(cellId);
            }
        }

        // Update the mesh
        //
        // Adaption data structure contains the information needed for mapping
        // the fields from the previous mesh to the updated mesh.
        adaptionData = mesh.adaptionPrepare(trackAdaptation);

        // Save values needed for remapping
        //
        // Data of cells no longer in the mesh are moved in a temporary
        // container to be ued later for remapping.
        PiercedVector<double, long> previousScalarField;
        PiercedVector<std::array<double, 3>, long> previousVectorField;
        for (const adaption::Info &adaptionInfo : adaptionData) {
            // Consider only cell refinements
            if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_REFINEMENT) {
                continue;
            }

            // Save parent data
            for (long previousId : adaptionInfo.previous) {
                previousScalarField.insert(previousId, std::move(scalarField.at(previousId)));
                previousVectorField.insert(previousId, std::move(vectorField.at(previousId)));
            }

        }

        adaptionData = mesh.adaptionAlter(trackAdaptation, squeeshPatchStorage);

        // Map the fields on newly created cells
        //
        // On refinement one parent cell has been divided in many children
        // cells. At every children it will be assigned the value of their
        // parent.
        for (const adaption::Info &adaptionInfo : adaptionData) {
            // Consider only cell refinements
            if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_REFINEMENT) {
                continue;
            }

            // Assign data to children
            long parentId = adaptionInfo.previous.front();
            for (long currentId : adaptionInfo.current) {
                scalarField.set(currentId, previousScalarField.at(parentId));
                vectorField.set(currentId, previousVectorField.at(parentId));
            }
        }

        mesh.adaptionCleanup();

        // Write mesh
        mesh.write();
    }

    // Coarsening
    int nofCoarsening = 1;
    for (int i = 0; i < nofCoarsening; ++i) {
        // Mark cells that need coarsening
        //
        // All cells inside a circle centered in the origin and with a radius
        // equal to 0.15 will be coarsened.
        for (const Cell &cell : mesh.getCells()) {
            long cellId = cell.getId();
            std::array<double, 3> cellCentroid = mesh.evalCellCentroid(cellId);

            double r = std::sqrt(cellCentroid[0] * cellCentroid[0] + cellCentroid[1] * cellCentroid[1] + cellCentroid[2] * cellCentroid[2]);
            if (r < 0.15) {
                mesh.markCellForCoarsening(cellId);
            }
        }

        adaptionData = mesh.adaptionPrepare(trackAdaptation);

        // Save values needed for remapping
        //
        // Data of cells no longer in the mesh are moved in a temporary
        // container to be ued later for remapping.
        PiercedVector<double, long> previousScalarField;
        PiercedVector<std::array<double, 3>, long> previousVectorField;
        for (const adaption::Info &adaptionInfo : adaptionData) {
            // Consider only cell coarsenings
            if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_COARSENING) {
                continue;
            }

            // Save parent data
            for (long previousId : adaptionInfo.previous) {
                previousScalarField.insert(previousId, std::move(scalarField.at(previousId)));
                previousVectorField.insert(previousId, std::move(vectorField.at(previousId)));
            }
        }

        // Update the mesh
        //
        // Adaption data structure contains the information needed for mapping
        // the fields from the previous mesh to the updated mesh.
        adaptionData = mesh.adaptionAlter(trackAdaptation, squeeshPatchStorage);

        // Map the fields on newly created cells
        //
        // On coarsening many parent cell have been merged together to create
        // one child cell. The child fields will be set equal to the average
        // of the parent fields.
        for (const adaption::Info &adaptionInfo : adaptionData) {
            // Consider only cell coarsening
            if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_COARSENING) {
                continue;
            }

            // Evaluate parent average
            int nParents = adaptionInfo.previous.size();
            double scalarParentAverage = 0.;
            std::array<double, 3> vectorParentAverage = {{0., 0., 0.}};
            for (long parentId : adaptionInfo.previous) {
                scalarParentAverage += previousScalarField.at(parentId);
                vectorParentAverage += previousVectorField.at(parentId);
            }
            scalarParentAverage /= nParents;
            vectorParentAverage /= (double) nParents;

            // Assign data to the child
            long childId = adaptionInfo.current.front();
            scalarField.set(childId, scalarParentAverage);
            vectorField.set(childId, vectorParentAverage);
        }

        mesh.adaptionCleanup();

        // Write mesh
        mesh.write();
    }

#if ENABLE_MPI == 1
    // Finalize MPI
    MPI_Finalize();
#endif
}
