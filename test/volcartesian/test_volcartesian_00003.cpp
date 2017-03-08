/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
#include "bitpit_volcartesian.hpp"

using namespace bitpit;

int main(int argc, char *argv[]) {

#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

	log::manager().initialize(log::COMBINED);
	log::cout() << "Testing dump/restore of a Cartesian patch" << std::endl;

	std::array<double, 3> origin = {{-10., -10., -10.}};
	double length = 20;
	double dh = 0.5;

	int archiveVersion = 1;

	//
	// 2D Test
	//
	log::cout() << "  >> 2D Cartesian patch" << std::endl;

	// Create the patch
	log::cout() << "Creating 2D patch..." << std::endl;

	VolCartesian *patch_2D = new VolCartesian(0, 2, origin, length, dh);
	patch_2D->getVTK().setName("cartesian_uniform_patch_2D");
	patch_2D->setMemoryMode(VolCartesian::MEMORY_NORMAL);

	log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
	log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

	patch_2D->write();

	// Dump the patch
	log::cout() << "Dumping 2D patch..." << std::endl;

	std::string header2D = "2D Cartesian patch";
	OBinaryArchive binaryWriter2D("cartesian_uniform_patch_2D", archiveVersion, header2D);
	patch_2D->dump(binaryWriter2D.getStream());
	binaryWriter2D.close();

	// Delete the patch
	log::cout() << "Deleting 2D patch..." << std::endl;

	delete patch_2D;

	// Restore the patch
	log::cout() << "Restoring 2D patch..." << std::endl;

	VolCartesian *patch_2D_restored = new VolCartesian();
	IBinaryArchive binaryReader2D("cartesian_uniform_patch_2D");
	patch_2D_restored->restore(binaryReader2D.getStream());
	binaryReader2D.close();

	log::cout() << "Restored cell count:   " << patch_2D_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

	patch_2D_restored->getVTK().setName("cartesian_uniform_patch_2D_restored");
	patch_2D_restored->write();

	//
	// 3D Test
	//
	log::cout() << "  >> 3D Cartesian patch" << std::endl;

	// Create the patch
	log::cout() << "Creating 3D patch..." << std::endl;

	VolCartesian *patch_3D = new VolCartesian(1, 3, origin, length, dh);
	patch_3D->getVTK().setName("cartesian_uniform_patch_3D");
	patch_3D->setMemoryMode(VolCartesian::MEMORY_NORMAL);

	log::cout() << "Cell count:   " << patch_3D->getCellCount() << std::endl;
	log::cout() << "Vertex count: " << patch_3D->getVertexCount() << std::endl;

	patch_3D->write();

	// Dump the patch
	log::cout() << "Dumping 3D patch..." << std::endl;

	std::string header3D = "3D Cartesian patch";
	OBinaryArchive binaryWriter3D("cartesian_uniform_patch_3D", archiveVersion, header3D);
	patch_3D->dump(binaryWriter3D.getStream());
	binaryWriter3D.close();

	// Delete the patch
	log::cout() << "Deleting 3D patch..." << std::endl;

	delete patch_3D;

	// Restore the patch
	log::cout() << "Restoring 3D patch..." << std::endl;

	VolCartesian *patch_3D_restored = new VolCartesian();
	IBinaryArchive binaryReader3D("cartesian_uniform_patch_3D");
	patch_3D_restored->restore(binaryReader3D.getStream());
	binaryReader3D.close();

	log::cout() << "Restored cell count:   " << patch_3D_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count: " << patch_3D_restored->getVertexCount() << std::endl;

	patch_3D_restored->getVTK().setName("cartesian_uniform_patch_3D_restored");
	patch_3D_restored->write();

	//
	// Patch Manager test
	//

	// Dump all the patches
	log::cout() << "Dumping patch manager..." << std::endl;

	std::string headerPM = "2D and 3D Cartesian patch";
	OBinaryArchive binaryWriterPM("cartesian_uniform_patch_PM", archiveVersion, headerPM);
	patch::manager().dumpAll(binaryWriterPM.getStream());
	binaryWriterPM.close();

	// Delete old patches
	log::cout() << "Deleting existing patches..." << std::endl;

	delete patch_2D_restored;
	delete patch_3D_restored;

	// Restore all the patches
	log::cout() << "Restoring patches through patch manager..." << std::endl;

	VolCartesian *patch_2D_PM_restored = new VolCartesian();
	VolCartesian *patch_3D_PM_restored = new VolCartesian();
	IBinaryArchive binaryReaderPM("cartesian_uniform_patch_PM");
	patch::manager().restoreAll(binaryReaderPM.getStream());
	binaryReaderPM.close();

	log::cout() << "Restored cell count (2D):   " << patch_2D_PM_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count (2D): " << patch_2D_PM_restored->getVertexCount() << std::endl;
	log::cout() << "Restored cell count (3D):   " << patch_3D_PM_restored->getCellCount() << std::endl;
	log::cout() << "Restored vertex count (3D): " << patch_3D_PM_restored->getVertexCount() << std::endl;

	patch_2D_PM_restored->getVTK().setName("cartesian_uniform_patch_2D_restored_PM");
	patch_2D_PM_restored->write();

	patch_3D_PM_restored->getVTK().setName("cartesian_uniform_patch_3D_restored_PM");
	patch_3D_PM_restored->write();

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif

}
