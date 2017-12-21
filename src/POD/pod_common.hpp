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

#include "bitpit_containers.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

/*!
    \ingroup POD
    \brief The namespace 'pod' contains structures for working with the POD class.
*/
namespace pod {
typedef PiercedStorage<double> ScalarStorage;
typedef PiercedStorage<std::array<double, 3>> VectorStorage;

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
    SnapshotFile(const std::string _directory = std::string(), const std::string _name = std::string())
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
    bool                                    meshOwner;  /**<True if the mesh is released during clear. */
    VolumeKernel*                           mesh;   /**< Mesh of the field.*/
    std::unique_ptr<PiercedStorage<bool>>   mask;   /**< Mask: true where the values of the field are defined.*/
    std::unique_ptr<pod::ScalarStorage>     scalar; /**< Scalar fields.*/
    std::unique_ptr<pod::VectorStorage>     vector; /**< Vector fields.*/

    /**
     * Creates a new empty pod field.
     */
    PODField()
    {
        meshOwner = false;
        mesh   = nullptr;
        mask   = nullptr;
        scalar = nullptr;
        vector = nullptr;
    }

    /**
     * Copy a pod field.
     */
    PODField(const PODField & field) :
        mask(new PiercedStorage<bool>(*field.mask)),
        scalar(new PiercedStorage<double>(*field.scalar)),
        vector(new PiercedStorage<std::array<double,3>>(*field.vector))
    {
        meshOwner = field.meshOwner;
        mesh   = field.mesh;
    }

    /**
     * Destroy a pod field.
     */
    ~PODField()
    {
        clear();
    }

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
     * Set the linked kernel to a pod field.
     *
     * \param[in] lkernel Linked kernel
     */
    void setStaticKernel(const PiercedKernel<long> *lkernel)
    {
        mask->setStaticKernel(lkernel);
        scalar->setStaticKernel(lkernel);
        vector->setStaticKernel(lkernel);
    }

    /**
     * Set the linked mesh to a pod field.
     *
     * \param[in] lmesh Linked mesh
     */
    void setMesh(VolumeKernel *lmesh)
    {
        mesh = lmesh;
    }

    /**
     * Set the uniqueness of the mesh of a pod field.
     * If set the mesh will be destroyed with the pod field.
     *
     * \param[in] owner Owner of the mesh
     */
    void setMeshOwner(bool owner = true)
    {
        meshOwner = owner;
    }

    void clear()
    {
        if (meshOwner){
            delete mesh;
        }
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

    void clear()
    {
        scalar.reset();
        vector.reset();
    }

};

}

}
