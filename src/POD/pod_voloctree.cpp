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

#include "pod_voloctree.hpp"

namespace bitpit {

/**
 * \class PODVolOctree
 * \ingroup POD
 *
 * \brief The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 *
 * The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 */

/**
 * Creates a new PODVolOctree object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the pod object. MPI_COMM_WORLD is the default value.
 */
PODVolOctree::PODVolOctree(MPI_Comm comm):PODKernel(comm)
# else
PODVolOctree::PODVolOctree():PODKernel()
# endif
{
}

/**
 * Destructor of PODVolOctree
 */
PODVolOctree::~PODVolOctree()
{
}

/**
 * Create the mesh for POD.
 *
 * \return Pointer to instantiated mesh.
 */
VolumeKernel* PODVolOctree::createMesh()
{
    VolumeKernel *mesh = new VolOctree();

#if BITPIT_ENABLE_MPI
    mesh->setCommunicator(getCommunicator());
#endif

    return mesh;
}

/**
 * Map a POD Field TO the POD mesh.
 *
 * \param[in] field Input POD Field
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \return Mapped POD field on target cells of POD mesh.
 */
pod::PODField PODVolOctree::mapPODFieldToPOD(const pod::PODField & field, const std::unordered_set<long> * targetCells)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields on pod mesh
    std::size_t nsf = field.scalar->getFieldCount();
    std::size_t nvf = field.vector->getFieldCount();

    const PiercedStorage<mapping::Info> & m_mapper = getMeshMapper().getMapping();

    pod::PODField mappedField(nsf, nvf, getMesh(), &getMesh()->getCells());

    for (long id : *targetCells){
        std::size_t rawIndex = m_meshPOD->getCells().getRawIndex(id);

        double *datamappedS = mappedField.scalar->data(id);
        std::array<double,3> *datamappedV = mappedField.vector->data(id);

        if (m_mapper[id].type == adaption::Type::TYPE_RENUMBERING){
            bool dataB = field.mask->at(m_mapper[id].previous[0]);
            mappedField.mask->set(id, dataB);

            double *dataS = field.scalar->data(m_mapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *dataSi = dataS + i;
                double *datamappedSi = datamappedS + i;
                (*datamappedSi) = (*dataSi);
            }
            std::array<double,3> *dataV = field.vector->data(m_mapper[id].previous[0]);
            for (std::size_t i = 0; i < nvf; i++) {
                std::array<double,3> *dataVi = dataV + i;
                std::array<double,3> *datamappedVi = datamappedV + i;
                (*datamappedVi) = (*dataVi);
            }
        }
        else if (m_mapper[id].type == adaption::Type::TYPE_COARSENING){
            mappedField.mask->set(id, true);

            for (std::size_t i = 0; i < nsf; i++){
                double *datamappedSi = datamappedS + i;
                (*datamappedSi) = 0.0;
            }
            for (std::size_t i = 0; i < nvf; i++){
                std::array<double,3> *datamappedVi = datamappedV + i;
                (*datamappedVi) = std::array<double,3>{0.0, 0.0, 0.0};
            }

            bool dataB, dataMappedB = true;
            double *dataS;
            std::array<double,3> *dataV;
            double volmapped = getRawCellVolume(rawIndex);
            for (long idd : m_mapper[id].previous){
                dataB = field.mask->at(idd);
                dataMappedB &= dataB;

                dataS = field.scalar->data(idd);
                double vol = field.mesh->evalCellVolume(idd);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) += (*dataSi) * vol / volmapped;
                }
                dataV = field.vector->data(idd);
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) += (*dataVi) * vol / volmapped;
                }
            }
            mappedField.mask->set(id, dataMappedB);

        }
        else if (m_mapper[id].type == adaption::Type::TYPE_REFINEMENT){
            bool dataB = field.mask->at(m_mapper[id].previous[0]);
            double *dataS = field.scalar->data(m_mapper[id].previous[0]);
            std::array<double,3> *dataV = field.vector->data(m_mapper[id].previous[0]);

            mappedField.mask->set(id, dataB);

            for (std::size_t i = 0; i < nsf; i++){
                double *dataSi = dataS + i;
                double *datamappedSi = datamappedS + i;
                (*datamappedSi) = (*dataSi);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                std::array<double,3> *dataVi = dataV + i;
                std::array<double,3> *datamappedVi = datamappedV + i;
                (*datamappedVi) = (*dataVi);
            }
        }
    }

    return mappedField;

}

/**
 * Map a POD Field FROM the POD mesh.
 *
 * \param[out] field Output POD Field
 * \param[in] targetCells Pointer to set of target cells of output field mesh (cells with mapped values as result)
 * \param[in] mappedField POD field on target cells of output field mesh.
 */
void PODVolOctree::mapPODFieldFromPOD(pod::PODField & field, const std::unordered_set<long> * targetCells,
        const pod::PODField & mappedField)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : field.mesh->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields from pod mesh
    std::size_t nsf = field.scalar->getFieldCount();
    std::size_t nvf = field.vector->getFieldCount();

    const PiercedStorage<mapping::Info> & m_invmapper = getMeshMapper().getInverseMapping();

    for (long id : *targetCells){
        double *dataS = field.scalar->data(id);
        std::array<double,3> *dataV = field.vector->data(id);

        if (m_invmapper[id].type == adaption::Type::TYPE_RENUMBERING){
            bool datamappedB = mappedField.mask->at(m_invmapper[id].previous[0]);
            field.mask->set(id, datamappedB);
            double *datamappedS = mappedField.scalar->data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *dataSi = dataS + i;
                double *datamappedSi = datamappedS + i;
                (*dataSi) = (*datamappedSi);
            }
            std::array<double,3>  *datamappedV = mappedField.vector->data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nvf; i++) {
                std::array<double,3> *dataVi = dataV + i;
                std::array<double,3> *datamappedVi = datamappedV + i;
                (*dataVi) = (*datamappedVi);
            }
        }
        else if (m_invmapper[id].type == adaption::Type::TYPE_COARSENING){
            bool dataB = true;
            for (std::size_t i = 0; i < nsf; i++){
                double *dataSi = dataS + i;
                (*dataSi) = 0.0;
            }
            for (std::size_t i = 0; i < nvf; i++){
                std::array<double,3> *dataVi = dataV + i;
                (*dataVi) = std::array<double,3>{0.0, 0.0, 0.0};
            }
            bool datamappedB;
            double *datamappedS;
            std::array<double,3> *datamappedV;
            double vol = field.mesh->evalCellVolume(id);
            for (long idd : m_invmapper[id].previous){
                std::size_t rawIndexIdd = m_meshPOD->getCells().getRawIndex(idd);
                datamappedB = mappedField.mask->at(idd);
                dataB &= datamappedB;
                datamappedS = mappedField.scalar->data(idd);
                double volmapped = getRawCellVolume(rawIndexIdd);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*dataSi) += (*datamappedSi) * volmapped / vol;
                }
                datamappedV = mappedField.vector->data(idd);
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*dataVi) += (*datamappedVi) * volmapped / vol;
                }
            }
            field.mask->set(id, dataB);
        }
        else if (m_invmapper[id].type == adaption::Type::TYPE_REFINEMENT){
            bool datamappedB = mappedField.mask->at(m_invmapper[id].previous[0]);
            field.mask->set(id, datamappedB);
            double *datamappedS = mappedField.scalar->data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *dataSi = dataS + i;
                double *datamappedSi = datamappedS + i;
                (*dataSi) = (*datamappedSi);
            }
            std::array<double,3> *datamappedV = mappedField.vector->data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nvf; i++) {
                std::array<double,3> *dataVi = dataV + i;
                std::array<double,3> *datamappedVi = datamappedV + i;
                (*dataVi) = (*datamappedVi);
            }
        }
    }
}

/**
 * Map a storage of fields TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 * \return Fields storage on target cells of POD mesh.
 */
PiercedStorage<double> PODVolOctree::mapFieldsToPOD(const PiercedStorage<double> & fields, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells,
        const std::vector<std::size_t> &scalarIds,
        const std::vector<std::array<std::size_t, 3>> &vectorIds)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields on pod mesh
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();

    const PiercedStorage<mapping::Info> & m_mapper = getMeshMapper().getMapping();

    PiercedStorage<double> mappedFields(fields.getFieldCount(), &getMesh()->getCells());

    for (long id : *targetCells){
        std::size_t rawIndex = m_meshPOD->getCells().getRawIndex(id);

        double *datamapped = mappedFields.data(id);
        if (m_mapper[id].type == adaption::Type::TYPE_RENUMBERING){
            const double *data = fields.data(m_mapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                const double *datai = data + scalarIds[i];
                double *datamappedi = datamapped + scalarIds[i];
                (*datamappedi) = (*datai);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    const double *datai = data + vectorIds[i][j];
                    double *datamappedi = datamapped + vectorIds[i][j];
                    (*datamappedi) = (*datai);
                }
            }
        }
        else if (m_mapper[id].type == adaption::Type::TYPE_COARSENING){
            for (std::size_t i = 0; i < nsf; i++){
                double *datamappedi = datamapped + scalarIds[i];
                (*datamappedi) = 0.0;
            }
            for (std::size_t i = 0; i < nvf; i++){
                for (std::size_t j = 0; j < 3; j++) {
                    double *datamappedi = datamapped + vectorIds[i][j];
                    (*datamappedi) = 0.0;
                }
            }
            double volmapped = getRawCellVolume(rawIndex);
            for (long idd : m_mapper[id].previous){
                const double *data = fields.data(idd);
                double vol = mesh->evalCellVolume(idd);

                for (std::size_t i = 0; i < nsf; i++){
                    const double *datai = data + scalarIds[i];
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) += (*datai) * vol / volmapped;
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        const double *datai = data + vectorIds[i][j];
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) += (*datai) * vol / volmapped;
                    }
                }
            }
        }
        else if (m_mapper[id].type == adaption::Type::TYPE_REFINEMENT){
            const double *data = fields.data(m_mapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                const double *datai = data + scalarIds[i];
                double *datamappedi = datamapped + scalarIds[i];
                (*datamappedi) = (*datai);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    const double *datai = data + vectorIds[i][j];
                    double *datamappedi = datamapped + vectorIds[i][j];
                    (*datamappedi) = (*datai);
                }
            }
        }
    }

    return mappedFields;

}

/**
 * Map a storage of fields FROM the POD mesh.
 *
 * \param[out] field Output Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of input field mesh mesh (cells with mapped values as result)
 * \param[in]mappedFields  Fields storage on POD mesh to be mapped.
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 */
void PODVolOctree::mapFieldsFromPOD(PiercedStorage<double> & fields, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells, const PiercedStorage<double> & mappedFields,
        const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : mesh->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields from pod mesh
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();

    const PiercedStorage<mapping::Info> & m_invmapper = getMeshMapper().getInverseMapping();

    for (long id : *targetCells){
        double *data = fields.data(id);
        if (m_invmapper[id].type == adaption::Type::TYPE_RENUMBERING){
            const double *datamapped = mappedFields.data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                const double *datamappedi = datamapped + scalarIds[i];
                (*datai) = (*datamappedi);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    const double *datamappedi = datamapped + vectorIds[i][j];
                    (*datai) = (*datamappedi);
                }
            }
        }
        else if (m_invmapper[id].type == adaption::Type::TYPE_COARSENING){
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                (*datai) = 0.0;
            }
            for (std::size_t i = 0; i < nvf; i++){
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    (*datai) = 0.0;
                }
            }
            double vol = mesh->evalCellVolume(id);
            for (long idd : m_invmapper[id].previous){
                std::size_t rawIndexIdd = m_meshPOD->getCells().getRawIndex(idd);
                const double * datamapped = mappedFields.data(idd);
                double volmapped = getRawCellVolume(rawIndexIdd);

                for (std::size_t i = 0; i < nsf; i++){
                    double *datai = data + scalarIds[i];
                    const double *datamappedi = datamapped + scalarIds[i];
                    (*datai) += (*datamappedi) * volmapped / vol;
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datai = data + vectorIds[i][j];
                        const double *datamappedi = datamapped + vectorIds[i][j];
                        (*datai) += (*datamappedi) * volmapped / vol;
                    }
                }
            }
        }
        else if (m_invmapper[id].type == adaption::Type::TYPE_REFINEMENT){
            const double *datamapped = mappedFields.data(m_invmapper[id].previous[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                const double *datamappedi = datamapped + scalarIds[i];
                (*datai) = (*datamappedi);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    const double *datamappedi = datamapped + vectorIds[i][j];
                    (*datai) = (*datamappedi);
                }
            }
        }
    }

}

/**
 * Map a storage of boolean field TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \return Fields storage on target cells of POD mesh.
 */
PiercedStorage<bool> PODVolOctree::mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells)
{

    // Map bool field (number of fields = 1) on pod mesh
    PiercedStorage<bool> mappedField(1, &getMesh()->getCells());
    mapBoolFieldToPOD(field, mesh, targetCells, mappedField);
    return mappedField;

}

/**
 * Map a storage of boolean field TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \param[out] Fields storage on target cells of POD mesh.
 */
void PODVolOctree::mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells, PiercedStorage<bool> & mappedField)
{
    BITPIT_UNUSED(mesh);

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map bool field (number of fields = 1) on pod mesh
    const PiercedStorage<mapping::Info> & m_mapper = getMeshMapper().getMapping();

    mappedField.setStaticKernel(&getMesh()->getCells());
    mappedField.fill(false);

    for (long id : *targetCells){
        if (m_mapper[id].type == adaption::Type::TYPE_RENUMBERING){
            bool dataB = field.at(m_mapper[id].previous[0]);
            mappedField.set(id, dataB);
        }
        else if (m_mapper[id].type == adaption::Type::TYPE_COARSENING){
            mappedField.set(id, false);
            bool dataB, dataMappedB = true;
            for (long idd : m_mapper[id].previous){
                dataB = field.at(idd);
                dataMappedB &= dataB;
            }
            mappedField.set(id, dataMappedB);
        }
        else if (m_mapper[id].type == adaption::Type::TYPE_REFINEMENT){
            bool dataB = field.at(m_mapper[id].previous[0]);
            mappedField.set(id, dataB);
        }
    }
}

/**
 * Map a set of cells TO the POD mesh by using a pre-computed mapper.
 *
 * \param[in] targetCells Pointer to input cells set
 * return Target cells set of POD mesh.
 */
std::unordered_set<long> PODVolOctree::mapCellsToPOD(const std::unordered_set<long> * targetCells)
{
    std::unordered_set<long> mappedCells;
    const PiercedStorage<mapping::Info> & m_invmapper = getMeshMapper().getInverseMapping();
    for (const long & id : *targetCells){
        if (m_invmapper[id].type == adaption::Type::TYPE_RENUMBERING){
            mappedCells.insert(m_invmapper[id].previous[0]);
        }
        if (m_invmapper[id].type == adaption::Type::TYPE_REFINEMENT){
            mappedCells.insert(m_invmapper[id].previous[0]);
        }
        if (m_invmapper[id].type == adaption::Type::TYPE_COARSENING){
            for (long idd : m_invmapper[id].previous)
                mappedCells.insert(idd);
        }
    }

    return mappedCells;

}

void PODVolOctree::adaptMeshToMesh(VolumeKernel* meshToAdapt, VolumeKernel * meshReference)
{
    BITPIT_UNUSED(meshToAdapt);

    computeMapper(meshReference);

    bool adapt = true;

    while(adapt){

        adapt = false;

        const PiercedStorage<mapping::Info> & m_mapper = getMeshMapper().getMapping();

        for (Cell & cell : getMesh()->getCells()){
            long id = cell.getId();
            if (m_mapper[id].type == adaption::Type::TYPE_COARSENING){
                getMesh()->markCellForRefinement(id);
                adapt = true;
            }
        }

#if BITPIT_ENABLE_MPI==1
        MPI_Allreduce(MPI_IN_PLACE, &adapt, 1, MPI_C_BOOL, MPI_LOR, m_communicator);
#endif

        std::vector<adaption::Info> infoAdapt = getMesh()->adaptionPrepare(true);

        getMeshMapper().mappingAdaptionPreparare(infoAdapt, true);

        infoAdapt = getMesh()->adaptionAlter(true);

        getMeshMapper().mappingAdaptionUpdate(infoAdapt, true, true);

        getMesh()->adaptionCleanup();

    }

}

}
