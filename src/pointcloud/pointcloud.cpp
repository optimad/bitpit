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

#include "bitpit_common.hpp"

#include "pointcloud.hpp"

namespace bitpit {

/*!
    \class PointCloud
    \ingroup pointpatches

    \brief The PointCloud class defines a point cloud.

    PointCloud defines a point cloud.
*/

#if BITPIT_ENABLE_MPI==1
/*!
    Creates an uninitialized partitioned patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
PointCloud::PointCloud(MPI_Comm communicator)
    : PointKernel(communicator, true)
#else
/*!
    Creates an uninitialized serial patch.
*/
PointCloud::PointCloud()
    : PointKernel(true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param dimension is the dimension of the patch
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
PointCloud::PointCloud(int dimension, MPI_Comm communicator)
    : PointKernel(PatchManager::AUTOMATIC_ID, dimension, communicator, true)
#else
/*!
    Creates a patch.

    \param dimension is the dimension of the patch
*/
PointCloud::PointCloud(int dimension)
    : PointKernel(PatchManager::AUTOMATIC_ID, dimension, true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch.

    If a null comunicator is provided, a serial patch will be created, this
    means that each processor will be unaware of the existence of the other
    processes.

    \param id is the id of the patch
    \param dimension is the dimension of the patch
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
PointCloud::PointCloud(int id, int dimension, MPI_Comm communicator)
    : PointKernel(id, dimension, communicator, true)
#else
/*!
    Creates a patch.

    \param id is the id of the patch
    \param dimension is the dimension of the patch
*/
PointCloud::PointCloud(int id, int dimension)
    : PointKernel(id, dimension, true)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
    Creates a patch restoring the patch saved in the specified stream.

    The number of processes in the communicator should be equal to the number
    of processes of the communicator used when dumping the patch.

    \param stream is the stream to read from
    \param communicator is the communicator to be used for exchanging data
    among the processes
*/
PointCloud::PointCloud(std::istream &stream, MPI_Comm communicator)
    : PointKernel(communicator, false)
#else
/*!
    Creates a patch restoring the patch saved in the specified stream.

    \param stream is the stream to read from
*/
PointCloud::PointCloud(std::istream &stream)
    : PointKernel(false)
#endif
{
    // Restore the patch
    restore(stream);
}

/*!
    Creates a clone of the pach.

    \result A clone of the pach.
*/
std::unique_ptr<PatchKernel> PointCloud::clone() const
{
    return std::unique_ptr<PointCloud>(new PointCloud(*this));
}

/*!
 * Enables or disables expert mode.
 *
 * When expert mode is enabled, it will be possible to change the
 * patch using low level functions (e.g., it will be possible to
 * add individual cells, add vertices, delete cells, ...).
 *
 * \param expert if true, the expert mode will be enabled
 */
void PointCloud::setExpert(bool expert)
{
    PointKernel::setExpert(expert);
}

/*!
 *  Get the version associated to the binary dumps.
 *
 *  \result The version associated to the binary dumps.
 */
int PointCloud::_getDumpVersion() const
{
    const int DUMP_VERSION = 1;

    return DUMP_VERSION;
}

/*!
 *  Write the patch to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void PointCloud::_dump(std::ostream &stream) const
{
#if BITPIT_ENABLE_MPI==1
    // Dump works only for serial calculations
    if (getProcessorCount() != 1) {
        throw std::runtime_error ("Dump of lineunstructured is implemented only for serial calculations.");
    }
#endif

    // Save the vertices
    dumpVertices(stream);

    // Save the cells
    dumpCells(stream);

    // Save the interfaces
    dumpInterfaces(stream);
}

/*!
 *  Restore the patch from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void PointCloud::_restore(std::istream &stream)
{
#if BITPIT_ENABLE_MPI==1
    // Restore works only for serial calculations
    if (getProcessorCount() != 1) {
        throw std::runtime_error ("Restore of lineunstructured is implemented only for serial calculations.");
    }
#endif

    // Restore the vertices
    restoreVertices(stream);

    // Restore the cells
    restoreCells(stream);

    // Restore the interfaces
    restoreInterfaces(stream);
}

/*!
 * Locates the cell the contains the point.
 *
 * If the point is not inside the patch, the function returns the id of the
 * null element.
 *
 * NOTE: this function is not implemented yet.
 *
 * \param[in] point is the point to be checked
 * \result Returns the linear id of the cell the contains the point. If the
 * point is not inside the patch, the function returns the id of the null
 * element.
 */
long PointCloud::locatePoint(const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

    throw std::runtime_error ("The function 'locatePoint' is not implemented yet");

    return false;
}

}
