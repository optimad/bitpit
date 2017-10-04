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
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing dump/restore of a 2D octree patch.
*
* \param patch_2D is the patch that will be created by the test
* \param patch_2D_restored is the patch that will be restored by the test
*/
int subtest_001(VolOctree *patch_2D, VolOctree *patch_2D_restored)
{
    std::array<double, 3> origin = {{-8., -8., -8.}};
    double length = 16;
    double dh = 1;

    int archiveVersion = 1;

    log::cout() << "  >> 2D octree patch" << std::endl;

    // Create the patch
    log::cout() << "Creating 2D patch..." << std::endl;

    patch_2D = new VolOctree(2, origin, length, dh);
    patch_2D->getVTK().setName("octree_uniform_patch_2D");
    patch_2D->update();

    // Refine the patch
    std::vector<uint32_t> refineList2D;
    refineList2D.push_back(  7);
    refineList2D.push_back( 13);
    refineList2D.push_back( 15);
    refineList2D.push_back( 26);
    refineList2D.push_back( 27);
    refineList2D.push_back( 31);
    refineList2D.push_back( 37);
    refineList2D.push_back( 39);
    refineList2D.push_back( 49);
    refineList2D.push_back( 50);
    refineList2D.push_back( 51);
    refineList2D.push_back( 53);
    refineList2D.push_back( 55);
    refineList2D.push_back( 63);
    refineList2D.push_back( 78);
    refineList2D.push_back(100);
    refineList2D.push_back(102);
    refineList2D.push_back(105);
    refineList2D.push_back(108);
    refineList2D.push_back(109);
    refineList2D.push_back(110);
    refineList2D.push_back(135);
    refineList2D.push_back(141);
    refineList2D.push_back(143);
    refineList2D.push_back(146);
    refineList2D.push_back(147);
    refineList2D.push_back(151);
    refineList2D.push_back(153);
    refineList2D.push_back(154);
    refineList2D.push_back(155);
    refineList2D.push_back(157);
    refineList2D.push_back(159);
    refineList2D.push_back(165);
    refineList2D.push_back(167);
    refineList2D.push_back(183);
    refineList2D.push_back(198);
    refineList2D.push_back(204);
    refineList2D.push_back(206);
    refineList2D.push_back(225);
    refineList2D.push_back(228);
    refineList2D.push_back(229);
    refineList2D.push_back(230);

    for (long id : refineList2D) {
        patch_2D->markCellForRefinement(id);
    }

    for (int k = 0; k < 2; ++k) {
        std::vector<adaption::Info> adaptionData = patch_2D->update(true);
        for (auto &adaptionInfo : adaptionData) {
            if (adaptionInfo.entity != adaption::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_REFINEMENT) {
                continue;
            }

            for (long id : adaptionInfo.current) {
                patch_2D->markCellForRefinement(id);
            }
        }
    }
    patch_2D->update();

    // Show patch info
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->write();

    // Dump the patch
    log::cout() << "Dumping 2D patch..." << std::endl;

    std::string header2D = "2D octree patch";
    OBinaryArchive binaryWriter2D("octree_uniform_patch_2D", archiveVersion, header2D);
    patch_2D->dump(binaryWriter2D.getStream());
    binaryWriter2D.close();

    // Delete the patch
    log::cout() << "Deleting 2D patch..." << std::endl;

    delete patch_2D;

    // Restore the patch
    log::cout() << "Restoring 2D patch..." << std::endl;

    patch_2D_restored = new VolOctree();
    IBinaryArchive binaryReader2D("octree_uniform_patch_2D");
    patch_2D_restored->restore(binaryReader2D.getStream());
    binaryReader2D.close();

    log::cout() << "Restored cell count:   " << patch_2D_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count: " << patch_2D_restored->getVertexCount() << std::endl;

    patch_2D_restored->getVTK().setName("octree_uniform_patch_2D_restored");
    patch_2D_restored->write();

    return 0;
}

/*!
* Subtest 002
*
* Testing dump/restore of a 3D octree patch.
*
* \param patch_3D is the patch that will be created by the test
* \param patch_3D_restored is the patch that will be restored by the test
*/
int subtest_002(VolOctree *patch_3D, VolOctree *patch_3D_restored)
{
    std::array<double, 3> origin = {{-8., -8., -8.}};
    double length = 16;
    double dh = 1;

    int archiveVersion = 1;

    log::cout() << "  >> 3D octree patch" << std::endl;

    // Create the patch
    log::cout() << "Creating 3D patch..." << std::endl;

    patch_3D = new VolOctree(3, origin, length, dh);
    patch_3D->getVTK().setName("octree_uniform_patch_3D");
    patch_3D->update();

    // Refine the patch
    std::vector<uint32_t> refineList3D;
    refineList3D.push_back(2351);
    refineList3D.push_back(2365);
    refineList3D.push_back(2367);
    refineList3D.push_back(2422);
    refineList3D.push_back(2423);
    refineList3D.push_back(2431);
    refineList3D.push_back(2477);
    refineList3D.push_back(2479);
    refineList3D.push_back(2533);
    refineList3D.push_back(2534);
    refineList3D.push_back(2535);
    refineList3D.push_back(2541);
    refineList3D.push_back(2543);
    refineList3D.push_back(2559);
    refineList3D.push_back(2878);
    refineList3D.push_back(2988);
    refineList3D.push_back(2990);
    refineList3D.push_back(2997);
    refineList3D.push_back(3004);
    refineList3D.push_back(3005);
    refineList3D.push_back(3006);
    refineList3D.push_back(3375);
    refineList3D.push_back(3389);
    refineList3D.push_back(3391);
    refineList3D.push_back(3430);
    refineList3D.push_back(3431);
    refineList3D.push_back(3439);
    refineList3D.push_back(3445);
    refineList3D.push_back(3446);
    refineList3D.push_back(3447);
    refineList3D.push_back(3453);
    refineList3D.push_back(3455);
    refineList3D.push_back(3501);
    refineList3D.push_back(3503);
    refineList3D.push_back(3567);
    refineList3D.push_back(3886);
    refineList3D.push_back(3900);
    refineList3D.push_back(3902);
    refineList3D.push_back(4005);
    refineList3D.push_back(4012);
    refineList3D.push_back(4013);
    refineList3D.push_back(4014);

    for (long id : refineList3D) {
        patch_3D->markCellForRefinement(id);
    }

    for (int k = 0; k < 2; ++k) {
        std::vector<adaption::Info> adaptionData = patch_3D->update(true);
        for (auto &adaptionInfo : adaptionData) {
            if (adaptionInfo.entity != adaption::ENTITY_CELL) {
                continue;
            } else if (adaptionInfo.type != adaption::TYPE_REFINEMENT) {
                continue;
            }

            for (long id : adaptionInfo.current) {
                patch_3D->markCellForRefinement(id);
            }
        }
    }
    patch_3D->update();

    // Show patch info
    log::cout() << "Cell count:   " << patch_3D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_3D->getVertexCount() << std::endl;

    patch_3D->write();

    // Dump the patch
    log::cout() << "Dumping 3D patch..." << std::endl;

    std::string header3D = "3D octree patch";
    OBinaryArchive binaryWriter3D("octree_uniform_patch_3D", archiveVersion, header3D);
    patch_3D->dump(binaryWriter3D.getStream());
    binaryWriter3D.close();

    // Delete the patch
    log::cout() << "Deleting 3D patch..." << std::endl;

    delete patch_3D;

    // Restore the patch
    log::cout() << "Restoring 3D patch..." << std::endl;

    patch_3D_restored = new VolOctree();
    IBinaryArchive binaryReader3D("octree_uniform_patch_3D");
    patch_3D_restored->restore(binaryReader3D.getStream());
    binaryReader3D.close();

    log::cout() << "Restored cell count:   " << patch_3D_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count: " << patch_3D_restored->getVertexCount() << std::endl;

    patch_3D_restored->getVTK().setName("octree_uniform_patch_3D_restored");
    patch_3D_restored->write();

    return 0;
}

/*!
* Subtest 003
*
* Testing dump/restore through the patch manager.
*/
int subtest_003()
{
    int archiveVersion = 1;

    // Dump all the patches
    log::cout() << "Dumping patch manager..." << std::endl;

    std::string headerPM = "2D and 3D octree patch";
    OBinaryArchive binaryWriterPM("octree_uniform_patch_PM", archiveVersion, headerPM);
    patch::manager().dumpAll(binaryWriterPM.getStream());
    binaryWriterPM.close();

    // Restore all the patches
    log::cout() << "Restoring patches through patch manager..." << std::endl;

    VolOctree *patch_2D_PM_restored = new VolOctree();
    VolOctree *patch_3D_PM_restored = new VolOctree();
    IBinaryArchive binaryReaderPM("octree_uniform_patch_PM");
    patch::manager().restoreAll(binaryReaderPM.getStream());
    binaryReaderPM.close();

    log::cout() << "Restored cell count (2D):   " << patch_2D_PM_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count (2D): " << patch_2D_PM_restored->getVertexCount() << std::endl;
    log::cout() << "Restored cell count (3D):   " << patch_3D_PM_restored->getCellCount() << std::endl;
    log::cout() << "Restored vertex count (3D): " << patch_3D_PM_restored->getVertexCount() << std::endl;

    patch_2D_PM_restored->getVTK().setName("octree_uniform_patch_2D_restored_PM");
    patch_2D_PM_restored->write();

    patch_3D_PM_restored->getVTK().setName("octree_uniform_patch_3D_restored_PM");
    patch_3D_PM_restored->write();

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
    log::cout() << "Testing dump/restore of octree patches" << std::endl;

    int status;
    try {
        VolOctree *patch_2D = nullptr;
        VolOctree *patch_2D_restored = nullptr;
        VolOctree *patch_3D = nullptr;
        VolOctree *patch_3D_restored = nullptr;

        status = subtest_001(patch_2D, patch_2D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_002(patch_3D, patch_3D_restored);
        if (status != 0) {
            return status;
        }

        status = subtest_003();
        if (status != 0) {
            return status;
        }

        delete patch_2D;
        delete patch_2D_restored;
        delete patch_3D_restored;
        delete patch_3D;
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
