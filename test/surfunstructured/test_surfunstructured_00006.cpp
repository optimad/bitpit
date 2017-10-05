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
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_surfunstructured.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing dump/restore of a 2D octree patch.
*
* \param patch_2D is the patch that will be created by the test
* \param patch_2D_restored is the patch that will be restored by the test
*/
int subtest_001(SurfUnstructured *patch_2D, SurfUnstructured *patch_2D_restored)
{
	int archiveVersion = 1;

	log::cout() << "  >> 2D unstructured surface patch" << std::endl;

	// Create the patch
	log::cout() << "Creating 2D patch..." << std::endl;

	const std::string fielname_2D = "./data/cube.stl";

	patch_2D = new SurfUnstructured(2, 2);
	patch_2D->importSTL(fielname_2D);
	patch_2D->getVTK().setName("surfunstructured_patch_2D");

	log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
	log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

	patch_2D->write();

	// Dump the patch
	log::cout() << "Dumping 2D patch..." << std::endl;

	std::string header2D = "2D unstructured surface patch";
	OBinaryArchive binaryWriter2D("surfunstructured_patch_2D.dat", archiveVersion, header2D);
	patch_2D->dump(binaryWriter2D.getStream());
	binaryWriter2D.close();

	// Delete the patch
	log::cout() << "Deleting 2D patch..." << std::endl;

	// Restore the patch
	log::cout() << "Restoring 2D patch..." << std::endl;

	patch_2D_restored = new SurfUnstructured();
	IBinaryArchive binaryReader2D("surfunstructured_patch_2D.dat");
	patch_2D_restored->restore(binaryReader2D.getStream());
	binaryReader2D.close();

	log::cout() << "Restored cell count:   " << patch_2D_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

	patch_2D_restored->getVTK().setName("surfunstructured_patch_2D_restored");
	patch_2D_restored->write();

    return 0;
}

/*!
* Subtest 002
*
* Testing dump/restore through the patch manager.
*/
int subtest_002()
{
    int archiveVersion = 1;

    // Dump all the patches
	log::cout() << "Dumping patch manager..." << std::endl;

	std::string headerPM = "2D unstructured surface patch";
	OBinaryArchive binaryWriterPM("surfunstructured_patch_PM.dat", archiveVersion, headerPM);
	patch::manager().dumpAll(binaryWriterPM.getStream());
	binaryWriterPM.close();

	// Restore all the patches
	log::cout() << "Restoring patches through patch manager..." << std::endl;

	SurfUnstructured *patch_2D_PM_restored = static_cast<SurfUnstructured *>(patch::manager().get(0));
	patch_2D_PM_restored->reset();

	IBinaryArchive binaryReaderPM("surfunstructured_patch_PM.dat");
	patch::manager().restoreAll(binaryReaderPM.getStream());
	binaryReaderPM.close();

	log::cout() << "Restored cell count (2D):   " << patch_2D_PM_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count (2D): " << patch_2D_PM_restored->getVertexCount() << std::endl;

	patch_2D_PM_restored->getVTK().setName("surfunstructured_patch_2D_restored_PM");
	patch_2D_PM_restored->write();

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
    log::cout() << "Testing dump/restore of unstructured surface patches" << std::endl;

    int status;
    try {
        SurfUnstructured *patch_2D = nullptr;
        SurfUnstructured *patch_2D_restored = nullptr;

        status = subtest_001(patch_2D, patch_2D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_002();
        if (status != 0) {
            return status;
        }

        delete patch_2D;
        delete patch_2D_restored;
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
