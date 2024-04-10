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

#include "volume_kernel.hpp"

namespace bitpit {

/*!
	\class VolumeKernel
	\ingroup volumepatches

	\brief The VolumeKernel class provides an interface for defining
	volume patches.

	VolumeKernel is the base class for defining voulme patches.
*/

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(AdaptionMode adaptionMode)
	: PatchKernel(adaptionMode)
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
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(dimension, communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(int dimension, AdaptionMode adaptionMode)
	: PatchKernel(dimension, adaptionMode)
#endif
{
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a patch.

	If a null comunicator is provided, a serial patch will be created, this
	means that each processor will be unaware of the existence of the other
	processes.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param communicator is the communicator to be used for exchanging data
	among the processes. If a null comunicator is provided, a serial patch
	will be created
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, AdaptionMode adaptionMode)
	: PatchKernel(id, dimension, communicator, haloSize, adaptionMode)
#else
/*!
	Creates a patch.

	\param id is the id that will be assigned to the patch
	\param dimension is the dimension of the patch
	\param adaptionMode is the adaption mode that will be used for the patch
*/
VolumeKernel::VolumeKernel(int id, int dimension, AdaptionMode adaptionMode)
	: PatchKernel(id, dimension, adaptionMode)
#endif
{
}

/*!
	Get the codimension of the patch in the volume space.

	\result The codimension of the patch in the volume space.
*/
int VolumeKernel::getVolumeCodimension() const
{
	return 0;
}

/*!
	Get the codimension of the patch in the surface space.

	\result The codimension of the patch in the surface space.
*/
int VolumeKernel::getSurfaceCodimension() const
{
	return -1;
}

/*!
	Get the codimension of the patch in the line space.

	\result The codimension of the patch in the line space.
*/
int VolumeKernel::getLineCodimension() const
{
	return -2;
}

/*!
	Get the codimension of the patch in the point space.

	\result The codimension of the patch in the point space.
*/
int VolumeKernel::getPointCodimension() const
{
	return -3;
}

/*!
	Extracts the external envelope and appends it to the given patch.

	The external envelope is composed by all the free faces of the patch.

	\param[in,out] envelope is the patch to which the external envelope
	will be appended
*/
void VolumeKernel::extractEnvelope(SurfaceKernel &envelope) const
{
	PatchKernel::extractEnvelope(envelope);
}

/*!
	Checks if the specified point is inside the patch.

	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the patch, false otherwise.
 */
bool VolumeKernel::isPointInside(double x, double y, double z) const
{
	return isPointInside({{x, y, z}});
}

/*!
	Checks if the specified point is inside a cell.

	\param[in] id is the index of the cells
	\param[in] x is the x coordinate of the point
	\param[in] y is the y coordinate of the point
	\param[in] z is the z coordinate of the point
	\result Returns true if the point is inside the cell, false otherwise.
 */
bool VolumeKernel::isPointInside(long id, double x, double y, double z) const
{
	return isPointInside(id, {{x, y, z}});
}

/*!
	Get the counter-clockwise ordered list of vertex for the specified face.

	\param cell is the cell
	\param face is the face
	\result The counter-clockwise ordered list of vertex for the specified face.
*/
ConstProxyVector<long> VolumeKernel::getFaceOrderedVertexIds(const Cell &cell, int face) const
{
    ConstProxyVector<long> faceVertexIds = cell.getFaceVertexIds(face);
    if (areFaceVerticesOrdered(cell, face)) {
        return faceVertexIds;
    }

    std::size_t nFaceVertices = faceVertexIds.size();
    ConstProxyVector<long> orderedFaceVertexIds(ConstProxyVector<long>::INTERNAL_STORAGE, nFaceVertices);
    ConstProxyVector<long>::storage_pointer orderedFaceVertexIdsStorage = orderedFaceVertexIds.storedData();
    for (std::size_t k = 0; k < nFaceVertices; ++k) {
        int vertex = getFaceOrderedLocalVertex(cell, face, k);
        orderedFaceVertexIdsStorage[k] = faceVertexIds[vertex];
    }

    return orderedFaceVertexIds;
}

/*!
	Check if the vertices of the specified face are counter-clockwise ordered.

	\param cell is the cell
	\param face is the face
	\result Return true if the vertices of the specified face are counter-clockwise
	ordered, false otherwise.
*/
bool VolumeKernel::areFaceVerticesOrdered(const Cell &cell, int face) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 2) {
        return true;
    }

    // Check if vertices are ordered
    ElementType faceType = cell.getFaceType(face);
    switch (faceType) {

    case (ElementType::POLYGON):
    {
        return true;
    }

    default:
    {
        assert(faceType != ElementType::UNDEFINED);
        const Reference2DElementInfo &faceInfo = static_cast<const Reference2DElementInfo &>(ReferenceElementInfo::getInfo(faceType));
        return faceInfo.areVerticesCCWOrdered();
    }

    }
}

/*!
	Get the local index of the vertex occupying the n-th position in the counter-clockwise
	ordered list of vertex ids.

	\param cell is the cell
	\param face is the face
	\param[in] n is the requested position
	\result The local index of the vertex occupying the n-th position in the counter-clockwise
	ordered list of vertex ids.
*/
int VolumeKernel::getFaceOrderedLocalVertex(const Cell &cell, int face, std::size_t n) const
{
    // Early return for low-dimension elements
    if (getDimension() <= 2) {
        return n;
    }

    // Get the request index
    ElementType faceType = cell.getFaceType(face);
    switch (faceType) {

    case (ElementType::POLYGON):
    {
        return n;
    }

    default:
    {
        assert(faceType != ElementType::UNDEFINED);
        const Reference2DElementInfo &faceInfo = static_cast<const Reference2DElementInfo &>(ReferenceElementInfo::getInfo(faceType));
        return faceInfo.getCCWOrderedVertex(n);
    }

    }
}

}
