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

#ifndef __BITPIT_POD_COMMON_HPP__
#define __BITPIT_POD_COMMON_HPP__

#include "bitpit_containers.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

/*!
    \ingroup POD
    \brief The namespace 'pod' contains structures for working with the POD class.
 */
namespace pod {
typedef PiercedStorage<double> ScalarStorage;                   /**< Scalar Storage type definition. */
typedef PiercedStorage<std::array<double, 3>> VectorStorage;    /**< Vector Storage type definition. */

/**
 * \ingroup POD
 * \brief The SnapFile structure is used to store the file names inside POD classes.
 */
struct SnapshotFile
{
    std::string directory;   /**< Name of the directory. */
    std::string name;        /**< Name of the file. */

    /**
     * Creates a new object.#
     */
    SnapshotFile(const std::string &_directory = std::string(), const std::string &_name = std::string())
    : directory(_directory), name(_name)
    {
    }

};

/**
 * \ingroup POD
 *
 * \brief The PODfield structure is used to store the fields inside POD classes.
 */
struct PODField
{
    std::shared_ptr<VolumeKernel>           meshStorage; /**< Mesh owned by the field */
    VolumeKernel*                           mesh;        /**< Mesh of the field.*/
    std::unique_ptr<PiercedStorage<bool>>   mask;        /**< Mask: true where the values of the field are defined.*/
    std::unique_ptr<pod::ScalarStorage>     scalar;      /**< Scalar fields.*/
    std::unique_ptr<pod::VectorStorage>     vector;      /**< Vector fields.*/

    /**
     * Creates a new empty pod field.
     */
    PODField()
        : mesh(nullptr)
    {
    }

    /**
     * Copy a pod field.
     */
    PODField(const PODField & field) :
        mask(new PiercedStorage<bool>(*field.mask)),
        scalar(new PiercedStorage<double>(*field.scalar)),
        vector(new PiercedStorage<std::array<double,3>>(*field.vector))
    {
        meshStorage = field.meshStorage;
        mesh        = field.mesh;
    }

    /**
     * Move a field.
     */
    PODField(PODField && field) = default;

    /**
     * Creates a new empty pod field with fixed number of scalar and vector fields and optional linked kernel.
     *
     * \param[in] nsf is the number of scalar fields.
     * \param[in] nvf is the number of vector fields.
     * \param[in] lmesh Linked mesh.
     * \param[in] lkernel Linked kernel.
     */

    PODField(int nsf, int nvf, VolumeKernel *lmesh = nullptr, const PiercedKernel<long> *lkernel = nullptr)
    {
        mesh   = lmesh;
        mask   = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, lkernel));
        scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(nsf, lkernel));
        vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(nvf, lkernel));
    }

    /**
     * Assignement operator using copy and swap
     *
     * \param[in] other another pod field whose content is copied in this
     * pod field
     * \return the assigneed field
     */
    PODField & operator=(PODField other)
    {
        this->swap(other);

        return *this;
    }

    /**
     * Exchanges the content of the pod field by the content of other pod
     * field received in input.
     *
     * \param[in] other another pod field whose content is swapped with that of
     * this pod field.
     */
    void swap(PODField &other)
    {
        std::swap(meshStorage, other.meshStorage);
        std::swap(mesh, other.mesh);
        std::swap(mask, other.mask);
        std::swap(scalar, other.scalar);
        std::swap(vector, other.vector);
    }

    /**
     * Set the linked static kernel to a pod field.
     *
     * \param[in] lkernel Linked kernel
     */
    void setStaticKernel(const PiercedKernel<long> *lkernel)
    {
        mask->unsetKernel();
        mask->setStaticKernel(lkernel);
        scalar->unsetKernel();
        scalar->setStaticKernel(lkernel);
        vector->unsetKernel();
        vector->setStaticKernel(lkernel);
    }

    /**
     * Set the linked dynamic kernel to a pod field.
     *
     * \param[in] lkernel Linked kernel
     */
    void setDynamicKernel(PiercedKernel<long> *lkernel)
    {
        mask->setDynamicKernel(lkernel, PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
        scalar->setDynamicKernel(lkernel, PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
        vector->setDynamicKernel(lkernel, PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    }

    /**
     * Set the linked mesh to a pod field.
     *
     * \param[in] lmesh Linked mesh
     */
    void setMesh(VolumeKernel *lmesh)
    {
        meshStorage.reset();
        mesh = lmesh;
    }

    /**
     * Set the linked mesh to a pod field.
     *
     * \param[in] lmesh Linked mesh
     */
    void setMesh(const std::shared_ptr<VolumeKernel> &lmesh)
    {
        meshStorage = lmesh;
        mesh        = meshStorage.get();
    }

    /**
     *Clear the pod field.
     * If the field is set to be mesh owner the mesh is destroyed.
     */
    void clear()
    {
        meshStorage.reset();

        mesh = nullptr;
        mask.reset();
        scalar.reset();
        vector.reset();
    }

};

/**
 * \ingroup POD
 *
 * \brief The PODMode structure is used to store the modes inside pod classes.
 */
struct PODMode
{
    std::unique_ptr<pod::ScalarStorage> scalar; /**< Scalar fields.*/
    std::unique_ptr<pod::VectorStorage> vector; /**< Vector fields.*/

    /**
     * Creates a new empty pod mode.
     */
    PODMode()
    {
        scalar = nullptr;
        vector = nullptr;
    }

    /**
     * Creates a new empty pod mode with fixed number of scalar and vector fields and optional linked kernel.
     *
     * \param[in] nsf # of scalar fields.
     * \param[in] nvf # of vector fields
     * \param[in] lkernel Linked kernel.
     */
    PODMode(int nsf, int nvf, const PiercedKernel<long>* lkernel = nullptr)
    {
        scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(nsf, lkernel));
        vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(nvf, lkernel));
    }

    /**
     * Clear a pod mode.
     */
    void clear()
    {
        scalar.reset();
        vector.reset();
    }

};

}

}

#endif
