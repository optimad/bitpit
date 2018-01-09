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

#include <cassert>
#include <sstream>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>

#if BITPIT_ENABLE_MPI
#   include <mpi.h>
#endif

#include <lapacke.h>
#include "bitpit_voloctree.hpp"

#include "pod.hpp"
#include "pod_kernel.hpp"

namespace bitpit {

/**
 * \class POD
 * \ingroup POD
 * \brief The POD class provides an interface for defining POD object.
 *
 * POD is the base class for defining POD object. \n
 * This class provides the interface methods, as the base get/set methods. \n
 * To execute a pod work-flow of a POD object call the
 * method run().
 *
 */

/**
 * Creates a new POD object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the pod object. MPI_COMM_WORLD is the default value.
 */
POD::POD(MPI_Comm comm)
# else
POD::POD()
# endif
: m_filter(1), m_sensorMask(1)
{

#if BITPIT_ENABLE_MPI
    m_communicator = MPI_COMM_NULL;
#endif

    m_podkernel = nullptr;
    m_meshType = MeshType::UNDEFINED;
    m_staticMesh = true;    //default static mesh TODO change it
    m_toUpdate = false;       //default, switched to true during set<...>
    m_nSnapshots = 0;
    m_nModes = std::numeric_limits<std::size_t>::max();
    m_energyLevel = 100;
    m_nScalarFields = 0;
    m_nVectorFields = 0;
    m_nReconstructionSnapshots = 0;

# if BITPIT_ENABLE_MPI
    initializeCommunicator(comm);
    MPI_Comm_size(m_communicator, &m_nProcs);
    MPI_Comm_rank(m_communicator, &m_rank);
#else
    m_rank = 0;
    m_nProcs = 1;
#endif

    m_memoryMode = MemoryMode::MEMORY_NORMAL;
    m_runMode = RunMode::COMPUTE;
    m_writeMode = WriteMode::DUMP;
    m_reconstructionMode = ReconstructionMode::MINIMIZATION;

    m_directory = ".";
    m_name = "pod";

}

/**
 * Destructor of POD
 */
POD::~POD()
{
    clear();
}

/**
 * Cleaning POD
 */
void POD::clear()
{

    m_meshType = MeshType::UNDEFINED;
    m_staticMesh = true;

    m_correlationMatrices.clear();
    m_modes.clear();
    m_minimizationMatrices.clear();

    m_database.clear();
    m_databaseReconstruction.clear();
    m_toUpdate = false;
    m_nSnapshots = 0;
    m_nScalarFields = 0;
    m_nVectorFields = 0;
    m_nModes = 0;
    m_nReconstructionSnapshots = 0;

# if BITPIT_ENABLE_MPI
    freeCommunicator();
# endif
}

/**
 * Set the name of the pod object.
 *
 * \param[in] name is the name that will be set.
 */
void POD::setName(const std::string &name)
{
    m_name = name;
}

/**
 * Get the name of the pod object.
 *
 * \return The name of the POD object.
 */
const std::string & POD::getName()
{
    return m_name;
}

/**
 * Set the input/output directory.
 *
 * \param[in] directory is the name of the directory.
 */
void POD::setDirectory(const std::string &directory)
{
    m_directory = directory + "/";
}

/**
 * Get the directory to write/read POD results.
 *
 * \return The name of the input/output object.
 */
const std::string & POD::getDirectory()
{
    return m_directory;
}

/**
 * Add a snapshot to the POD database.
 *
 * \param[in] directory is the directory that contains the snapshot.
 * \param[in] name is the name of the snapshot.
 */
void POD::addSnapshot(const std::string &directory, const std::string &name)
{
    pod::SnapshotFile file(directory, name);
    addSnapshot(file);
}

/**
 * Add a snapshot to the POD database.
 *
 * \param[in] name is the name of the snapshot.
 */
void POD::addSnapshot(const pod::SnapshotFile &file)
{
    m_database.push_back(file);
    m_nSnapshots++;
}

/**
 * Set the snapshots of the POD database.
 *
 * \param[in] database Vector of snapshot file names.
 */
void POD::setSnapshots(const std::vector<pod::SnapshotFile> &database)
{
    m_database = database;
    m_nSnapshots = database.size();
    m_toUpdate = true;
}

/**
 * Add a snapshot to the reconstruction database.
 *
 * \param[in] directory is the directory that contains the snapshot.
 * \param[in] name is the name of the snapshot.
 */
void POD::addReconstructionSnapshot(const std::string &directory, const std::string &name)
{
    pod::SnapshotFile file(directory, name);
    addReconstructionSnapshot(file);
}

/**
 * Add a snapshot to the reconstruction database.
 *
 * \param[in] name is the name of the snapshot.
 */
void POD::addReconstructionSnapshot(const pod::SnapshotFile &file)
{
    m_databaseReconstruction.push_back(file);
    m_nReconstructionSnapshots++;
}

/**
 * Set the number of retained POD modes during POD evaluation.
 *
 * \param[in] nmodes Number of retained modes.
 */
void POD::setModeCount(std::size_t nmodes)
{
    m_nModes = nmodes;
    m_toUpdate = true;
}

/**
 * Get the number of POD modes.
 *
 * \return The number of POD modes.
 */
std::size_t POD::getModeCount()
{
    return m_nModes;
}

/**
 * Set the energy percentage level (e%) used to retain a limited number of POD modes
 * during POD evaluation. m_nModes such that sum(lambda[1:m_nModes])/etot >= e% and sum(lambda[1:m_nModes-1])/etot < e%.
 *
 * \param[in] energy is the target energy percentage level
 */
void POD::setEnergyLevel(double energy)
{
    m_energyLevel = energy;
}

/**
 * Get the percentage energy level of retained POD modes.
 *
 * \return The energy level.
 */
double POD::getEnergyLevel()
{
    return m_energyLevel;
}

/**
 * Set the type of the POD mesh, where the fields and modes are defined.
 *
 * \param[in] type Type of the POD mesh.
 */
void POD::setMeshType(POD::MeshType type)
{

    if (m_meshType != MeshType::UNDEFINED)
        throw std::runtime_error ("POD mesh type already set. Change not allowed.");

    m_meshType = type;

    switch (type)
    {
    case POD::MeshType::VOLOCTREE:
#if BITPIT_ENABLE_MPI
        m_podkernel = std::unique_ptr<PODKernel>(new PODVolOctree(m_communicator));
#endif
        m_podkernel = std::unique_ptr<PODKernel>(new PODVolOctree());
        break;

    default:
        break;
    }
}

/**
 * Get the type of the POD mesh.
 *
 * \return Type of the POD mesh.
 */
POD::MeshType POD::getMeshType()
{
    return m_meshType;
}

/**
 * Set if the mesh, where the fields are defined,
 * is a Static mesh, fixed for the whole database, or it is
 * an AMR mesh (differently adapted for each field in the database).
 *
 * \param[in] flag if set to true the mesh is considered static.
 */
void POD::setStaticMesh(bool flag)
{
    m_staticMesh = flag;
}

/**
 * Set the memory mode of the POD object.
 *
 * \param[in] mode Memory mode. If set to MEMORY_LIGHT the POD modes if present
 * are erased and then they are read from file if needed. If set to MEMORY_NORMAL,
 * the POD modes are kept in memory when computed or read for the first time.
 */
void POD::setMemoryMode(POD::MemoryMode mode)
{
    if (m_memoryMode == mode) {
        return;
    }

    switch (mode) {

    case MemoryMode::MEMORY_NORMAL:

        m_memoryMode = MemoryMode::MEMORY_NORMAL;
        break;

    case MemoryMode::MEMORY_LIGHT:

        m_memoryMode = MemoryMode::MEMORY_LIGHT;
        m_modes.clear();
        break;

    }
}

/**
 * Get the memory mode of the POD object.
 */
POD::MemoryMode POD::getMemoryMode()
{
    return m_memoryMode;
}

/**
 * Set the run mode of the POD object.
 *
 * \param[in] mode is the run mode. If set to COMPUTE POD modes, mean field
 * and POD mesh are computed during execute command,
 * if set to RESTORE they are restored from file.
 */
void POD::setRunMode(POD::RunMode mode)
{
    if (m_runMode == mode) {
        return;
    }

    switch (mode) {

    case RunMode::COMPUTE:

        m_runMode = RunMode::COMPUTE;
        break;

    case RunMode::RESTORE:

        m_runMode = RunMode::RESTORE;
        break;

    }
}

/**
 * Get the execution mode of the POD object.
 */
POD::RunMode POD::getRunMode()
{
    return m_runMode;
}

/**
 * Set the write mode of the POD object.
 *
 * \param[in] mode is the write mode.
 */
void POD::setWriteMode(POD::WriteMode mode)
{
    if (m_writeMode == mode)
        return;

    m_writeMode = mode;
}

/**
 * Get the write mode of the POD object.
 */
POD::WriteMode POD::getWriteMode()
{
    return m_writeMode;
}

/**
 * Set the reconstruction mode of the POD object.
 *
 * \param[in] mode is the reconstruction mode. If set to MINIMIZATION solve the least-squares problem,
 * if set to PROJECTION uses orthogonal projection.
 */
void POD::setReconstructionMode(POD::ReconstructionMode mode)
{
    if (m_reconstructionMode == mode)
        return;

    switch (mode) {

    case ReconstructionMode::MINIMIZATION:

        m_reconstructionMode = ReconstructionMode::MINIMIZATION;
        break;

    case ReconstructionMode::PROJECTION:

        m_reconstructionMode = ReconstructionMode::PROJECTION;
        break;

    }
}

/**
 * Get the reconstruction mode of the POD object.
 *
 * \return The reconstruction mode.
 */
POD::ReconstructionMode POD::getReconstructionMode()
{
    return m_reconstructionMode;
}

/**
 * Set the boolean field used as mask for the POD reconstruction.
 * \param[in] mask Reference to boolean field used as mask to chose the active cells.
 * If set to true, the corresponding cell is used during the minimization [or projection],
 * otherwise it is discarded.
 */
void POD::setSensorMask(const PiercedStorage<bool> &mask)
{

    m_sensorMask.setStaticKernel(&(m_podkernel->getMesh()->getCells()));
    for (auto &cell : m_podkernel->getMesh()->getCells()) {
        long id = cell.getId();
        m_sensorMask[id] = mask[id];
    }

    fillListID(m_sensorMask);
    m_toUpdate = true;

}

/**
 * Get the number of snapshots in the database.
 *
 * \return The number of snapshots.
 */
std::size_t POD::getSnapshotCount()
{
    return m_nSnapshots;
}

/**
 * Get the names of the scalar fields involved in the POD.
 *
 * \return The names of scalar fields.
 */
std::vector<std::string> POD::getScalarNames()
{
    return m_nameScalarFields;
}

/**
 * Get the names of the vector fields (three components) involved in the POD.
 *
 * \return The names of vector fields.
 */
std::vector<std::array<std::string,3>> POD::getVectorNames()
{
    return m_nameVectorFields;
}

/**
 * Get the names of the fields involved in the POD in one list (scalar fields and components of vector fields).
 *
 * \return The names of fields.
 */
std::vector<std::string> POD::getFieldsNames()
{
    std::vector<std::string> names;
    for (std::string ss : m_nameScalarFields)
        names.push_back(ss);
    for (std::array<std::string,3> ssa : m_nameVectorFields){
        for (std::string ss : ssa)
            names.push_back(ss);
    }

    return names;
}

/**
 * Get a pointer to the POD mesh.
 *
 * \return The const pointer to the POD mesh.
 */
const VolumeKernel* POD::getMesh()
{
    if (m_podkernel == nullptr)
        throw std::runtime_error ("POD mesh not built");

    return m_podkernel->getMesh();
}

/**
 * Get a reference to the POD mean field.
 *
 * \return The const reference to the POD mean field.
 */
const pod::PODMode & POD::getMean()
{
    return m_mean;
}

/**
 * Get a reference to the POD modes.
 *
 * \return The const reference to the POD modes.
 */
const std::vector<pod::PODMode> & POD::getModes()
{
    return m_modes;
}

/**
 * Get the POD coefficients of the last reconstructed field.
 *
 * \return The POD modal coefficients of the last reconstructed field.
 */
std::vector<std::vector<double> > POD::getReconstructionCoeffs()
{
    return m_reconstructionCoeffs;
}

/**
 * Get the IDs of active cells.
 */
const std::vector<long> & POD::getListID()
{
    return m_listActiveIDs;
}

/**
 * Get the number of internal cells in the list of IDs of active cells.
 */
std::size_t POD::getListIDInternalCount()
{
    return m_sizeInternal;
}

/**
 * Execution of POD procedure.
 * Currently it computes the complete POD for static mesh: mean field of the snapshots,
 * correlation matrix, POD modes and reconstruction of database fields with
 * the retained modes.
 * Some reconstruction errors are evaluated (not yet).
 */
void POD::run()
{
    if (m_runMode == RunMode::COMPUTE) {
        // Evaluate mean field and mesh
        log::cout() << "pod : computing mean field and pod mesh... " << std::endl;
        evalMeanMesh();

        // Initialize ID list
        fillListID(m_filter);

        // Evaluate correlation matrices
        log::cout() << "pod : computing correlation matrix... " << std::endl;
        evalCorrelation();

        // Compute eigenvalues and eigenvectors
        log::cout() << "pod : computing eigenvectors... " << std::endl;
        evalEigen();

        // Compute POD modes
        log::cout() << "pod : computing pod modes... " << std::endl;
        evalModes();
    } else if (m_runMode == RunMode::RESTORE) {
        log::cout() << "pod : restore... " << std::endl;
        restore();
    }

    // Dump pod
    if (m_writeMode != WriteMode::NONE) {
        log::cout() << "pod : dumping... " << std::endl;
        dump();
    }

    // Evaluate reconstruction
    evalReconstruction();

}

/**
 * Evaluation of the mean field and the POD mesh, as the more refined envelope of the fluid
 * domain between all the meshes related to the database snapshots.
 * Currently implemented for static mesh only.
 */
void POD::evalMeanMesh()
{
    if (m_database.empty())
        throw std::runtime_error ("PODVolOctree database empty");

    // If the mesh is static do nothing on mesh and compute mean fields directly
    if (m_staticMesh) {
        evalMeanStaticMesh();
    } else{
        throw std::runtime_error("Dynamic mesh is not supported yet.");
    }

    // Set mesh POD to write only internal cells
#if BITPIT_ENABLE_MPI
    m_podkernel->getMesh()->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#else
    m_podkernel->getMesh()->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL);
#endif
}

/**
 * Evaluation of the mean field on static mesh (used as POD mesh).
 */
void POD::evalMeanStaticMesh()
{
    //Read first mesh and use as meshPOD
    VolumeKernel* meshr = m_podkernel->readMesh(m_database[0]);
    m_podkernel->setMesh(meshr);

    // Read first snapshot to set n scalar and vector fields TODO readNsfNvf only from one field
    // NOT initialize mean fields with first snapshot
    pod::PODField readf;
    readSnapshot(m_database[0], m_podkernel->getMesh(), readf);
    m_mean.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
    m_mean.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
    m_mean.scalar->fill(0.0);
    m_mean.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});

    m_filter.setStaticKernel(&(m_podkernel->getMesh()->getCells()));
    m_filter.fill(true);

    // Compute mean (starting from first snapshot) & filter+mask
    for (std::size_t i = 0; i < m_nSnapshots; i++) {
        log::cout() << "pod : evaluation mean - use snapshot " << i+1 << "/" << m_nSnapshots << std::endl;
        pod::PODField readf;
        readSnapshot(m_database[i], m_podkernel->getMesh(), readf);
        for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
            m_filter[cell.getId()] = (m_filter[cell.getId()]) && (readf.mask->at(cell.getId())); //N.B. mask == filter
            for (std::size_t j = 0; j < m_nScalarFields; j++)
                m_mean.scalar->at(cell.getId(),j) += readf.scalar->at(cell.getId(),j) / double(m_nSnapshots);
            for (std::size_t j = 0; j < m_nVectorFields; j++)
                m_mean.vector->at(cell.getId(),j) += readf.vector->at(cell.getId(),j) / double(m_nSnapshots);
        }
    }

    setSensorMask(m_filter);
}

/**
 * Fill list of active ID. Use boolean field passed as constant reference.
 * Fill with internal and ghosts cells.
 *
 * \param[in] mask is the Reference to the boolean field used as filter/mask to chose the active cells.
 */
void POD::fillListID(const PiercedStorage<bool> &mask)
{
    m_listActiveIDs.clear();
    std::vector<long> listGhost;
    for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
        long id = cell.getId();
        //TODO reset to commented when solid/fluid field on ghosts is dumped
        //if (mask[id]) {
        if (mask[id] && cell.isInterior()) {
            if (cell.isInterior()) {
                m_listActiveIDs.push_back(id);
            }else{
                listGhost.push_back(id);
            }
        }
    }

    m_sizeInternal = m_listActiveIDs.size();
    m_listActiveIDs.insert(m_listActiveIDs.end(), listGhost.begin(), listGhost.end());
}

/**
 * Evaluation of correlation matrix.
 */
void POD::evalCorrelation()
{
    initCorrelation();
    for (std::size_t i = 0; i < m_nSnapshots; i++){
        pod::PODField snapi;
        if (m_staticMesh)
            readSnapshot(m_database[i], m_podkernel->getMesh(), snapi);
        else
            readSnapshot(m_database[i], snapi);

        //TODO diff if dynamic mesh?...
        diff(snapi, m_mean);

        for (std::size_t j = i; j < m_nSnapshots; j++){
            log::cout() << "pod : evaluation " << i << "," << j << " term of correlation matrix " << std::endl;
            pod::PODField snapj;
            if (m_staticMesh)
                readSnapshot(m_database[j], m_podkernel->getMesh(), snapj);
            else
                readSnapshot(m_database[j], snapj);

            //TODO diff if dynamic mesh?...
            diff(snapj, m_mean);

            evalCorrelationTerm(i, snapi, j, snapj);

            if (!m_staticMesh)
                delete [] snapj.mesh;
        }
        if (!m_staticMesh)
            delete [] snapi.mesh;
    }

# if BITPIT_ENABLE_MPI
    for (std::size_t i = 0; i < m_nFields; i++ )
        MPI_Allreduce(MPI_IN_PLACE, m_correlationMatrices[i].data(), m_nSnapshots*m_nSnapshots, MPI_DOUBLE, MPI_SUM, m_communicator);
# endif

}

/**
 * Initialization of the correlation matrix.
 */
void POD::initCorrelation()
{
    m_correlationMatrices.clear();
    m_correlationMatrices.resize(m_nFields,std::vector<double>(m_nSnapshots*m_nSnapshots, 0.0));
}

/**
 * Evaluation of the correlation matrix.
 */
void POD::evalCorrelationTerm(int i, pod::PODField & snapi, int j, pod::PODField & snapj)
{

    if (m_staticMesh) {
        std::vector<long>::iterator it;
        for (it = m_listActiveIDs.begin(); it != m_listActiveIDs.begin()+m_sizeInternal; it++){
            long id = *it;
            if (m_nScalarFields){
                double* datasi = snapi.scalar->data(id);
                double* datasj = snapj.scalar->data(id);
                for (std::size_t ifield = 0; ifield < m_nScalarFields; ifield++){
                    m_correlationMatrices[ifield][i*m_nSnapshots+j] += (*datasi)*(*datasj)*snapi.mesh->evalCellVolume(id);
                    datasi++;
                    datasj++;
                }
            }
            if (m_nVectorFields){
                std::array<double,3>* datavi = snapi.vector->data(id);
                std::array<double,3>* datavj = snapj.vector->data(id);
                for (std::size_t ifield = m_nScalarFields; ifield < m_nFields; ifield++){
                    m_correlationMatrices[ifield][i*m_nSnapshots+j] += dotProduct((*datavi),(*datavj))*snapi.mesh->evalCellVolume(id);
                    datavi++;
                    datavj++;
                }
            }
        }

    } else {
        throw std::runtime_error("Dynamic mesh is not supported yet.");
    }
}

/**
 * Evaluation of snapshots reconstruction. Use the reconstruction database to read
 * the snapshots to be reconstructed.
 */
void POD::evalReconstruction()
{
    m_reconstructionCoeffs.clear();
    m_reconstructionCoeffs.resize(m_nFields, std::vector<double>(m_nModes,0.0));

    for (std::size_t i = 0; i < m_nReconstructionSnapshots; i++){
        pod::PODField snapi, reconi, erri;
        if (m_staticMesh)
            readSnapshot(m_databaseReconstruction[i], m_podkernel->getMesh(), snapi);
        else
            readSnapshot(m_databaseReconstruction[i], snapi);

        //TODO: evaluate l2norm
        reconstructFields(snapi, reconi);
        std::string rname = m_name + ".recon." + std::to_string(i);

        if (m_writeMode != WriteMode::NONE)
            dumpField(rname, reconi);

        //        {
        //            //TODO TEMPORARY IMPLEMENTATION
        //            //error field
        //            pod::PODField erri;
        //            erri.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
        //            erri.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
        //            erri.scalar->fill(0.0);
        //            erri.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
        //            for (long id : m_listActiveIDs){
        //                if (m_nScalarFields){
        //                    double* errs = erri.scalar->data(id);
        //                    double* recons = reconi.scalar->data(id);
        //                    double* vals = snapi.scalar->data(id);
        //                    for (std::size_t ifs = 0; ifs < m_nScalarFields; ifs++){
        //                        *errs  = std::abs(*recons - *vals);
        //                        errs++;
        //                        recons++;
        //                        vals++;
        //                    }
        //                }
        //
        //                if (m_nVectorFields){
        //                    std::array<double,3>* errv = erri.vector->data(id);
        //                    std::array<double,3>* reconv = reconi.vector->data(id);
        //                    std::array<double,3>* valv = snapi.vector->data(id);
        //                    for (std::size_t ifv = 0; ifv < m_nVectorFields; ifv++){
        //                        for (std::size_t j=0; j<3; j++){
        //                            (*errv)[j] = std::abs((*reconv)[j] - (*valv)[j]);
        //                        }
        //                        errv++;
        //                        reconv++;
        //                        valv++;
        //                    }
        //                }
        //            }
        //
        //            std::string errname = m_name + ".error." + std::to_string(i);
        //            dumpField(errname, erri, meshi);
        //
        //        }

        if (!m_staticMesh)
            delete [] snapi.mesh;

        // one coeff file for each reconstruction, each line is related to a field, each term is related to a mode
        std::ofstream outr(m_directory + m_name + ".recon.coeff."+ std::to_string(i) +".dat", std::ofstream::out);
        for (std::size_t j = 0; j < m_nFields; j++)
            outr << m_reconstructionCoeffs[j] << std::endl;

        outr.close();
    }
}

/**
 * Field reconstruction.
 * \param[in] field Original field.
 * \param[out] reconi Reconstructed field.
 */
void POD::reconstructFields(pod::PODField & field,  pod::PODField & reconi)
{
    //set the mesh of reconstructed field
    reconi.mesh = field.mesh;
    evalReconstructionCoeffs(field);
    buildFields(reconi);
}

/**
 * Evaluation of POD reconstruction coefficients.
 *
 * \param[in] snapi Original snapshot field.
 */
void POD::evalReconstructionCoeffs(pod::PODField & field)
{
    if (m_staticMesh)
        evalReconstructionCoeffsStaticMesh(field);
    else
        throw std::runtime_error("Dynamic mesh is not supported yet.");

# if BITPIT_ENABLE_MPI
    for (std::size_t i = 0; i < m_nFields; i++)
        MPI_Allreduce(MPI_IN_PLACE, m_reconstructionCoeffs[i].data(), m_nModes, MPI_DOUBLE, MPI_SUM, m_communicator);
# endif

}

/**
 * Evaluation of POD reconstruction coefficients on static mesh.
 *
 * \param[in] field Original field to be projected.
 */
void POD::evalReconstructionCoeffsStaticMesh(pod::PODField & field)
{

    diff(field, m_mean);
    m_reconstructionCoeffs.clear();
    m_reconstructionCoeffs.resize(m_nFields, std::vector<double>(m_nModes, 0.0));
    std::vector<std::vector<double> > rhs=m_reconstructionCoeffs;

    // Evaluate minimization matrices
    evalMinimizationMatrices();

    // Evaluate RHS
    for (std::size_t ir = 0; ir < m_nModes; ++ir){
        if (m_memoryMode == POD::MEMORY_LIGHT)
            readMode(ir);

        std::vector<long>::iterator it;
        for (it = m_listActiveIDs.begin(); it != m_listActiveIDs.begin()+m_sizeInternal; ++it){
            long id = *it;
            if (m_nScalarFields){
                double* modes = m_modes[ir].scalar->data(id);
                double* datas = field.scalar->data(id);
                for (std::size_t ifield = 0; ifield < m_nScalarFields; ifield++){
                    rhs[ifield][ir] += (*datas)*(*modes)*m_podkernel->getMesh()->evalCellVolume(id);
                    modes++;
                    datas++;
                }
            }
            if (m_nVectorFields){
                std::array<double,3>* modev = m_modes[ir].vector->data(id);
                std::array<double,3>* datav = field.vector->data(id);
                for (std::size_t ifield = m_nScalarFields; ifield < m_nFields; ifield++){
                    rhs[ifield][ir] += dotProduct((*datav),(*modev))*m_podkernel->getMesh()->evalCellVolume(id);
                    modev++;
                    datav++;
                }
            }
        }
    }

# if BITPIT_ENABLE_MPI
    for (std::size_t i = 0; i < m_nFields; i++)
        MPI_Allreduce(MPI_IN_PLACE, rhs[i].data(), m_nModes, MPI_DOUBLE, MPI_SUM, m_communicator);
# endif

    // Solve minimization
    solveMinimization(rhs);

    // Sum field and mean
    sum(field, m_mean);
}

/**
 * Reconstruct a field through POD. Use the m_reconstructionCoeffs already set.
 *
 * \param[out] reconi Resulting snapshot field.
 */
void POD::buildFields(pod::PODField & recon)
{
    recon.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
    recon.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
    recon.scalar->fill(0.0);
    recon.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});

    for (std::size_t ir = 0; ir < m_nModes; ir++){
        if (m_memoryMode == POD::MEMORY_LIGHT)
            readMode(ir);

        for (long id : m_listActiveIDs){
            if (m_nScalarFields){
                double* recons = recon.scalar->data(id);
                double* modes = m_modes[ir].scalar->data(id);
                for (std::size_t ifs = 0; ifs < m_nScalarFields; ifs++){
                    *recons  = *recons + *modes*m_reconstructionCoeffs[ifs][ir];
                    recons++;
                    modes++;
                }
            }

            if (m_nVectorFields){
                std::array<double,3>* reconv = recon.vector->data(id);
                std::array<double,3>* modev = m_modes[ir].vector->data(id);
                for (std::size_t ifv = 0; ifv < m_nVectorFields; ifv++){
                    *reconv = *reconv + *modev*m_reconstructionCoeffs[m_nScalarFields+ifv][ir];
                    reconv++;
                    modev++;
                }
            }
        }

        if (m_memoryMode == POD::MEMORY_LIGHT)
            m_modes[ir].clear();

    }

    sum(recon, m_mean);
}

/**
 * Initialize the minimization matrices.
 */
void POD::initMinimization()
{
    m_minimizationMatrices.clear();
    m_minimizationMatrices.resize(m_nFields,std::vector<double>(m_nModes*m_nModes, 0.0));
}

/**
 * Create or updated of the minimization matrices.
 */
void POD::evalMinimizationMatrices()
{
    if (m_toUpdate){
        initMinimization();


        if (m_reconstructionMode == ReconstructionMode::PROJECTION){
            for (std::size_t ifield = 0; ifield < m_nFields; ifield++){
                for (std::size_t ir = 0; ir < m_nModes; ++ir)
                    m_minimizationMatrices[ifield][ir*m_nModes+ir] = 1.0;
            }

        }else{

            for (std::size_t ir = 0; ir < m_nModes; ++ir){
                if (m_memoryMode == POD::MEMORY_LIGHT)
                    readMode(ir);

                for (std::size_t jr = ir; jr < m_nModes; ++jr){
                    if (m_memoryMode == POD::MEMORY_LIGHT)
                        readMode(jr);

                    std::vector<long>::iterator it;
                    for (it = m_listActiveIDs.begin(); it != m_listActiveIDs.begin()+m_sizeInternal; ++it){
                        long id = *it;
                        if (m_nScalarFields){
                            double* datasi = m_modes[ir].scalar->data(id);
                            double* datasj = m_modes[jr].scalar->data(id);
                            for (std::size_t ifield = 0; ifield < m_nScalarFields; ifield++){
                                m_minimizationMatrices[ifield][ir*m_nModes+jr] += (*datasi)*(*datasj)*m_podkernel->getMesh()->evalCellVolume(id);
                                datasi++;
                                datasj++;
                            }

                            if (ir != jr){
                                for (std::size_t ifield = 0; ifield < m_nScalarFields; ifield++){
                                    m_minimizationMatrices[ifield][jr*m_nModes+ir] = m_minimizationMatrices[ifield][ir*m_nModes+jr]; // sym: just in case...
                                }
                            }
                        }
                        if (m_nVectorFields){
                            std::array<double,3>* datavi = m_modes[ir].vector->data(id);
                            std::array<double,3>* datavj = m_modes[jr].vector->data(id);
                            for (std::size_t ifield = m_nScalarFields; ifield < m_nFields; ifield++){
                                m_minimizationMatrices[ifield][ir*m_nModes+jr] += dotProduct((*datavi),(*datavj))*m_podkernel->getMesh()->evalCellVolume(id);
                                datavi++;
                                datavj++;
                            }

                            if (ir != jr){
                                for (std::size_t ifield = m_nScalarFields; ifield < m_nFields; ifield++)
                                    m_minimizationMatrices[ifield][jr*m_nModes+ir] = m_minimizationMatrices[ifield][ir*m_nModes+jr];  // sym: just in case...
                            }
                        }
                    }
                }
            }

# if BITPIT_ENABLE_MPI
            for (std::size_t i = 0; i < m_nFields; i++ )
                MPI_Allreduce(MPI_IN_PLACE, m_minimizationMatrices[i].data(), m_nModes*m_nModes, MPI_DOUBLE, MPI_SUM, m_communicator);
# endif
        }
        m_toUpdate = false;
    }

}

/**
 * Solution of the minimization problem.
 *
 * \param[in] rhs right-hand side.
 */
void POD::solveMinimization(std::vector<std::vector<double> > & rhs)
{
    bool runSolver;
#if BITPIT_ENABLE_MPI
    runSolver = (m_rank == 0);
#else
    runSolver = true;
#endif

    if (runSolver) {
        for (std::size_t i = 0; i < m_nFields; ++i) {
            double A[m_nModes * m_nModes];
            int ipiv[m_nModes];

            for (std::size_t j = 0; j < m_nModes*m_nModes; ++j)
                A[j] = m_minimizationMatrices[i][j];

            lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, m_nModes, 1, A, m_nModes, ipiv, rhs[i].data(), m_nModes); //lapacke
            if (info != 0) {
                log::cout() << "WARNING!    algorithm convergence info " << info << std::endl;
                throw std::runtime_error("Unable to solve minimization problem.");
            }

            m_reconstructionCoeffs[i] = rhs[i];
        }
    }
}

/**
 * Evaluation of eigenvalues and eigenvectors of the correlation matrix.
 */
void POD::evalEigen()
{
    m_lambda.clear();
    m_lambda.resize(m_nFields,std::vector<double>(m_nSnapshots, 0.0));

    m_podCoeffs.clear();
    m_podCoeffs.resize(m_nFields, std::vector<std::vector<double>>(m_nSnapshots, std::vector<double>(m_nSnapshots,0.0)));

    bool runSolver;
#if BITPIT_ENABLE_MPI
    runSolver = (m_rank == 0);
#else
    runSolver = true;
#endif

    if (runSolver) {
        int N = m_nSnapshots*m_nSnapshots;
        for (std::size_t i = 0; i < m_nFields; ++i) {
            double Marr[N];
            double alambda[m_nSnapshots];
            for (std::size_t j = 0; j < m_nSnapshots*m_nSnapshots; ++j)
                Marr[j/m_nSnapshots+j%m_nSnapshots*m_nSnapshots] = m_correlationMatrices[i][j];

            lapack_int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', m_nSnapshots, Marr, m_nSnapshots, alambda); //lapacke
            if (info != 0) {
                log::cout() << "WARNING!    algorithm convergence info " << info << std::endl;
                throw std::runtime_error("Unable to solve minimization problem.");
            }

            // Check number of modes with energy level
            checkModeCount(alambda, i);

            for (std::size_t n = 0; n < m_nSnapshots; ++n)
                m_lambda[i][n] = alambda[m_nSnapshots - n - 1];

            for (std::size_t n = 0; n < m_nSnapshots; ++n) {
                for (std::size_t p = 0; p < m_nSnapshots; ++p)
                    m_podCoeffs[i][p][n] = Marr[N - m_nSnapshots * (p + 1) + n] * std::sqrt(abs(m_lambda[i][p]));
            }
        }

        // one coeff file for each mode, each line is related to a field, each term is related to a snapshot
        for (std::size_t ir=0; ir<m_nModes; ++ir) {
            std::ofstream out(m_directory + m_name + ".coeff."+std::to_string(ir)+".dat", std::ofstream::out);
            for (std::size_t i = 0; i < m_nFields; i++)
                out << m_podCoeffs[i][ir] << std::endl;

            out.close();
        }

        std::ofstream outl(m_directory + m_name + ".lambda.dat", std::ofstream::out);
        for (std::size_t i = 0; i < m_nFields; i++)
            outl << m_lambda[i] << std::endl;

        outl.close();
    }

    long bufferSize = m_nModes;
#if BITPIT_ENABLE_MPI
    MPI_Bcast(&bufferSize, 1, MPI_LONG, 0, m_communicator);
    m_nModes = bufferSize;
    for (std::size_t i = 0; i < m_nFields; i++) {
        for (std::size_t p=0; p<m_nModes; ++p)
            MPI_Bcast(m_podCoeffs[i][p].data(), m_nSnapshots, MPI_DOUBLE, 0, m_communicator);

        MPI_Bcast(m_lambda[i].data(), m_nSnapshots, MPI_DOUBLE, 0, m_communicator);
    }
#endif
}

/**
 * Check the number of modes with the number of snapshots and level of energy.
 *
 * \param[in] lambda Array of eigenvalues in ascending order.
 */
void POD::checkModeCount(double *alambda, std::size_t ifield)
{
    if (ifield == 0){
        _m_nr.clear();
        m_nModes = std::min(m_nModes, m_nSnapshots-1);
    }

    std::size_t nr = 0;
    double en = 0;
    double sumen = 0;
    std::vector<double> dlambda(m_nSnapshots);
    for (std::size_t n=0; n<m_nSnapshots; ++n){
        dlambda[n] = alambda[m_nSnapshots-n-1];
        sumen += std::abs(dlambda[n]);
    }
    while (nr < m_nModes && en < m_energyLevel){
        en += std::abs(dlambda[nr]) / sumen * 100;
        nr++;
    }
    _m_nr.push_back(nr);
    if (ifield == m_nFields-1)
        m_nModes = *(std::max_element(_m_nr.begin(), _m_nr.end()));
}

/**
 * Evaluation of POD modes.
 */
void POD::evalModes()
{
    m_modes.resize(m_nModes);

    // If m_staticMesh no need mapping on fields to compute POD modes
    if (m_staticMesh)
        evalModesStaticMesh();
    else
        throw std::runtime_error("Dynamic mesh is not supported yet.");
}

/**
 * Evaluation of POD modes for static mesh.
 */
void POD::evalModesStaticMesh()
{
    for (std::size_t ir = 0; ir < m_nModes; ir++){
        m_modes[ir].scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
        m_modes[ir].vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
        m_modes[ir].scalar->fill(0.0);
        m_modes[ir].vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
    }

    for (std::size_t is = 0; is < m_nSnapshots; is++){
        pod::PODField snapi;
        readSnapshot(m_database[is], m_podkernel->getMesh(), snapi);
        diff(snapi, m_mean);

        for (std::size_t ir = 0; ir < m_nModes; ir++){
            for (long id : m_listActiveIDs){
                if (m_nScalarFields){
                    double* modes = m_modes[ir].scalar->data(id);
                    double* datas = snapi.scalar->data(id);
                    for (std::size_t ifs = 0; ifs < m_nScalarFields; ifs++){
                        *modes = *modes + *datas*m_podCoeffs[ifs][ir][is] / m_lambda[ifs][ir];
                        modes++;
                        datas++;
                    }
                }
                if (m_nVectorFields){
                    std::array<double, 3>* modev = m_modes[ir].vector->data(id);
                    std::array<double, 3>* datav = snapi.vector->data(id);
                    for (std::size_t ifv = 0; ifv < m_nVectorFields; ifv++){
                        *modev = *modev + *datav*m_podCoeffs[m_nScalarFields+ifv][ir][is] / m_lambda[m_nScalarFields+ifv][ir];
                        modev++;
                        datav++;
                    }
                }
            }
            if (m_memoryMode == POD::MEMORY_LIGHT){
                dumpMode(ir);
                m_modes[ir].clear();
            }
        }
    }
}

/**
 * Read a snapshot from file defined on a given mesh.
 *
 * \param[in] snap Snapshot filename
 * \param[in] mesh Pointer to POD mesh
 * \param[out] fieldr Resulting POD field
 */
void POD::readSnapshot(pod::SnapshotFile snap, VolumeKernel* mesh, pod::PODField & fieldr)
{
    //TODO volumekernel
    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    std::string filename = std::string(snap.directory) + "/" + std::string(snap.name) + ".data";
    IBinaryArchive binaryReader(filename, dumpBlock);
    std::istream &dataStream = binaryReader.getStream();

    //Set mesh field to input mesh
    fieldr.mesh = mesh;

    // Restore solved cells
    fieldr.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &mesh->getCells()));
    fieldr.mask->restore(dataStream);

    // For each snapshot the global number of fields are read & updated again .

    // Restore scalar fields
    utils::binary::read(dataStream, m_nScalarFields);
    m_nameScalarFields.resize(m_nScalarFields);
    for (std::size_t i = 0; i < m_nScalarFields; ++i)
        utils::binary::read(dataStream, m_nameScalarFields[i]);

    if (m_nScalarFields) {
        fieldr.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &mesh->getCells()));
        fieldr.scalar->restore(dataStream);
    }

    // Restore vector fields
    utils::binary::read(dataStream, m_nVectorFields);
    m_nFields = m_nScalarFields + m_nVectorFields;
    m_nameVectorFields.resize(m_nVectorFields);
    for (std::size_t i = 0; i < m_nVectorFields; ++i)
        utils::binary::read(dataStream, m_nameVectorFields[i]);

    if (m_nVectorFields) {
        fieldr.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &mesh->getCells()));
        fieldr.vector->restore(dataStream);
    }
    binaryReader.close();

    // Write vtk with field in DEBUG mode
    if (m_writeMode == WriteMode::DEBUG) {
        std::vector<std::string > datastring;
        mesh->getVTK().setDirectory(m_directory);
        mesh->getVTK().setName(m_name + "." + snap.name);
        std::vector<uint8_t> mask(mesh->getCellCount());
        std::vector<std::vector<double>> fields(m_nScalarFields, std::vector<double>(mesh->getCellCount()));
        std::vector<std::vector<std::array<double,3>>> fieldv(m_nVectorFields, std::vector<std::array<double,3>>(mesh->getCellCount()));

        int i=0;
        for (const Cell &cell : mesh->getCells()) {
            long id = cell.getId();
            mask[i] = uint8_t(fieldr.mask->at(id));
            ++i;
        }
        mesh->getVTK().addData("mask", VTKFieldType::SCALAR, VTKLocation::CELL, mask);
        datastring.push_back("mask");

        for (std::size_t isf = 0; isf < m_nScalarFields; isf++) {
            int i = 0;
            for (const Cell &cell : mesh->getCells()) {
                long id = cell.getId();
                fields[isf][i] = fieldr.scalar->at(id, isf);
                ++i;
            }
            mesh->getVTK().addData(m_nameScalarFields[isf], VTKFieldType::SCALAR, VTKLocation::CELL, fields[isf]);
            datastring.push_back(m_nameScalarFields[isf]);
        }

        for (std::size_t ivf = 0; ivf < m_nVectorFields; ivf++) {
            int i = 0;
            for (const Cell &cell : mesh->getCells()) {
                long id = cell.getId();
                fieldv[ivf][i] = fieldr.vector->at(id, ivf);
                ++i;
            }
            std::string vname= m_nameVectorFields[ivf][0].substr(0,m_nameVectorFields[ivf][0].size()-2);
            mesh->getVTK().addData(vname, VTKFieldType::VECTOR, VTKLocation::CELL, fieldv[ivf]);
            datastring.push_back(vname);
        }

        // Set temporary write all cells if meshPOD
        if (mesh == m_podkernel->getMesh())
            mesh->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL);

        mesh->write();

        for (std::size_t is=0; is < datastring.size(); ++is)
            mesh->getVTK().removeData(datastring[is]);
    }

    // Reset temporary write internal cells only if meshPOD
    if (mesh == m_podkernel->getMesh()){
#if BITPIT_ENABLE_MPI
        mesh->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#else
        mesh->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL);
#endif
    }
}

/**
 * Read a snapshot from file. Read even the related mesh file.
 *
 * \param[in] snap Snapshot filename
 * \param[out] fieldr Resulting POD field
 */
void POD::readSnapshot(pod::SnapshotFile snap, pod::PODField & fieldr)
{
    //Read and set the mesh
    VolumeKernel* mesh = m_podkernel->readMesh(snap);
    fieldr.mesh = mesh;

    readSnapshot(snap, mesh, fieldr);
}

/**
 * Read a mode from file.
 *
 * \param[in] ir Index of POD mode to read. If ir == m_nModes read the mean field.
 */
void POD::readMode(std::size_t ir)
{
    if (!m_podkernel->getMesh())
        throw std::runtime_error ("POD: POD mesh not initialized before read a POD mode");

    if (ir > m_nModes)
        throw std::runtime_error ("POD: mode index > number of POD modes during read");

    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    if (ir == m_nModes){
        // Read mean
        std::string filename = m_directory + m_name + ".mean"+".data";
        IBinaryArchive binaryReader(filename, dumpBlock);
        std::istream &dataStream = binaryReader.getStream();

        // Restore solved cells (read for each POD mode)
        m_filter.restore(dataStream);

        // For each mode/mean the global number of fields are read & updated again .

        // Restore scalar fields
        utils::binary::read(dataStream, m_nScalarFields);
        m_nameScalarFields.resize(m_nScalarFields);
        for (std::size_t i = 0; i < m_nScalarFields; ++i)
            utils::binary::read(dataStream, m_nameScalarFields[i]);

        if (m_nScalarFields) {
            m_mean.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
            m_mean.scalar->restore(dataStream);
        }

        // Restore vector fields
        utils::binary::read(dataStream, m_nVectorFields);
        m_nameVectorFields.resize(m_nVectorFields);
        m_nFields = m_nScalarFields + m_nVectorFields;
        for (std::size_t i = 0; i < m_nVectorFields; ++i)
            utils::binary::read(dataStream, m_nameVectorFields[i]);

        if (m_nVectorFields) {
            m_mean.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
            m_mean.vector->restore(dataStream);
        }
        binaryReader.close();
    }
    else{
        // Read mode
        m_modes[ir].clear();

        std::string filename = m_directory + m_name + ".mode."+std::to_string(ir)+".data";
        IBinaryArchive binaryReader(filename, dumpBlock);
        std::istream &dataStream = binaryReader.getStream();

        // Restore solved cells (read for each POD mode)
        m_filter.restore(dataStream);

        // For each mode/mean the global number of fields are read & updated again .

        // Restore scalar fields
        utils::binary::read(dataStream, m_nScalarFields);
        m_nameScalarFields.resize(m_nScalarFields);
        for (std::size_t i = 0; i < m_nScalarFields; ++i)
            utils::binary::read(dataStream, m_nameScalarFields[i]);

        if (m_nScalarFields) {
            m_modes[ir].scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(m_nScalarFields, &m_podkernel->getMesh()->getCells()));
            m_modes[ir].scalar->restore(dataStream);
        }

        // Restore vector fields
        utils::binary::read(dataStream, m_nVectorFields);
        m_nFields = m_nScalarFields + m_nVectorFields;
        m_nameVectorFields.resize(m_nVectorFields);
        for (std::size_t i = 0; i < m_nVectorFields; ++i)
            utils::binary::read(dataStream, m_nameVectorFields[i]);

        if (m_nVectorFields) {
            m_modes[ir].vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(m_nVectorFields, &m_podkernel->getMesh()->getCells()));
            m_modes[ir].vector->restore(dataStream);
        }
        binaryReader.close();
    }
}

/**
 * Write a mode to file (dump).
 *
 * \param[in] ir Index of POD mode to write. If ir == m_nModes write the mean field.
 */
void POD::dumpMode(std::size_t ir)
{
    if (ir > m_nModes)
        throw std::runtime_error ("POD: mode index > number of pod modes during write");

    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    if (ir == m_nModes) {
        // Dumping mean as snapshot
        std::string header = "podmean." + std::to_string(ir);
        OBinaryArchive binaryWriter(m_directory + m_name + ".mean"+".data", ARCHIVE_VERSION, header, dumpBlock);
        std::ostream &dataStream = binaryWriter.getStream();
        m_filter.dump(dataStream);
        utils::binary::write(dataStream,std::size_t(m_nScalarFields));
        if (m_nScalarFields) {
            for (std::size_t i = 0; i < m_nScalarFields; i++)
                utils::binary::write(dataStream, m_nameScalarFields[i]);
            m_mean.scalar->dump(dataStream);
        }

        utils::binary::write(dataStream,std::size_t(m_nVectorFields));
        if (m_nVectorFields) {
            for (std::size_t i = 0; i < m_nVectorFields; i++)
                utils::binary::write(dataStream, m_nameVectorFields[i]);
            m_mean.vector->dump(dataStream);
        }
        binaryWriter.close();

    } else{
        // Dumping modes as snapshots
        std::string header = "podmode." + std::to_string(ir);
        OBinaryArchive binaryWriter(m_directory + m_name + ".mode."+std::to_string(ir)+".data", ARCHIVE_VERSION, header, dumpBlock);
        std::ostream &dataStream = binaryWriter.getStream();
        m_filter.dump(dataStream);
        utils::binary::write(dataStream,std::size_t(m_nScalarFields));
        if (m_nScalarFields) {
            for (std::size_t i = 0; i < m_nScalarFields; i++)
                utils::binary::write(dataStream, m_nameScalarFields[i]);

            m_modes[ir].scalar->dump(dataStream);
        }

        utils::binary::write(dataStream,std::size_t(m_nVectorFields));
        if (m_nVectorFields) {
            for (std::size_t i = 0; i < m_nVectorFields; i++)
                utils::binary::write(dataStream, m_nameVectorFields[i]);

            m_modes[ir].vector->dump(dataStream);
        }
        binaryWriter.close();
    }
}

/**
 * Write the POD results by dumping.
 */
void POD::dump()
{
    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    log::cout()  << "Dumping POD ... " << std::endl;
    // Dump the POD info (no reconstruction info)
    std::string header = "pod";
    OBinaryArchive binaryWriter(m_directory + m_name, ARCHIVE_VERSION, header, dumpBlock);
    std::ostream &dataStream = binaryWriter.getStream();
    utils::binary::write(dataStream, m_meshType);
    utils::binary::write(dataStream, m_nModes);
    utils::binary::write(dataStream, m_nSnapshots);
    for (std::size_t i = 0; i < m_nSnapshots; i++) {
        utils::binary::write(dataStream, m_database[i].directory);
        utils::binary::write(dataStream, m_database[i].name);
    }
    utils::binary::write(dataStream,m_staticMesh);
    utils::binary::write(dataStream,m_memoryMode);
    binaryWriter.close();

    // Dump the modes, mean and meshPOD
    if (m_memoryMode == MEMORY_NORMAL) {
        for (std::size_t ir = 0; ir < m_nModes; ir++)
            dumpMode(ir);

    }

    // Dumping mean as snapshot
    {
        dumpMode(m_nModes);
    }

    // Dumping mask
    {
        std::string header = "podmask";
        OBinaryArchive binaryWriter(m_directory + m_name + ".mask.data", ARCHIVE_VERSION, header, dumpBlock);
        std::ostream &dataStream = binaryWriter.getStream();
        m_sensorMask.dump(dataStream);
        binaryWriter.close();
    }

    // Dumping mesh
    {
        std::string header = "podmesh";
        OBinaryArchive binaryWriter(m_directory + m_name + ".mesh", ARCHIVE_VERSION, header, dumpBlock);
        m_podkernel->getMesh()->dump(binaryWriter.getStream());
        binaryWriter.close();
    }

    // Write debug info
    if (m_writeMode == WriteMode::DEBUG) {
        std::vector<std::string > datastring;
        m_podkernel->getMesh()->getVTK().setDirectory(m_directory);
        m_podkernel->getMesh()->getVTK().setName(m_name);
        std::vector<std::vector<std::vector<double>>> fields(m_nModes+1, std::vector<std::vector<double>> (m_nScalarFields, std::vector<double>(m_podkernel->getMesh()->getInternalCount(), 0.0)));
        std::vector<std::vector<std::vector<std::array<double,3>>>> fieldv(m_nModes+1, std::vector<std::vector<std::array<double,3>>>(m_nVectorFields, std::vector<std::array<double,3>>(m_podkernel->getMesh()->getInternalCount(), {{0.0, 0.0, 0.0}})));
        for (std::size_t isf = 0; isf < m_nScalarFields; isf++) {
            int i = 0;
            for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
                if (cell.isInterior()) {
                    long id = cell.getId();
                    fields[0][isf][i] = m_mean.scalar->at(id, isf);
                    ++i;
                }
            }
            m_podkernel->getMesh()->getVTK().addData("scalarField" + std::to_string(isf) + "_mean", VTKFieldType::SCALAR, VTKLocation::CELL, fields[0][isf]);
            datastring.push_back("scalarField" + std::to_string(isf) + "_mean");
        }

        for (std::size_t ivf = 0; ivf < m_nVectorFields; ivf++) {
            int i = 0;
            for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
                if (cell.isInterior()) {
                    long id = cell.getId();
                    fieldv[0][ivf][i] = m_mean.vector->at(id, ivf);
                    ++i;
                }
            }
            m_podkernel->getMesh()->getVTK().addData("vectorField"+std::to_string(ivf)+"_mean", VTKFieldType::VECTOR, VTKLocation::CELL, fieldv[0][ivf]);
            datastring.push_back("vectorField"+std::to_string(ivf)+"_mean");
        }

        for (std::size_t ir = 0; ir < m_nModes; ir++) {
            if (m_memoryMode == MEMORY_LIGHT)
                readMode(ir);

            for (std::size_t isf = 0; isf < m_nScalarFields; isf++) {
                int i=0;
                for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
                    if (cell.isInterior()) {
                        long id = cell.getId();
                        fields[ir+1][isf][i] = m_modes[ir].scalar->at(id, isf);
                        ++i;
                    }
                }
                m_podkernel->getMesh()->getVTK().addData("scalarField"+std::to_string(isf)+"_podMode"+std::to_string(ir), VTKFieldType::SCALAR, VTKLocation::CELL, fields[ir+1][isf]);
                datastring.push_back("scalarField"+std::to_string(isf)+"_podMode"+std::to_string(ir));
            }

            for (std::size_t ivf = 0; ivf < m_nVectorFields; ivf++) {
                int i=0;
                for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
                    if (cell.isInterior()) {
                        long id = cell.getId();
                        fieldv[ir+1][ivf][i] = m_modes[ir].vector->at(id, ivf);
                        ++i;
                    }
                }
                m_podkernel->getMesh()->getVTK().addData("vectorField"+std::to_string(ivf)+"_podMode"+std::to_string(ir), VTKFieldType::VECTOR, VTKLocation::CELL, fieldv[ir+1][ivf]);
                datastring.push_back("vectorField"+std::to_string(ivf)+"_podMode"+std::to_string(ir));
            }

            if (m_memoryMode == MEMORY_LIGHT)
                m_modes[ir].clear();
        }

        std::vector<uint8_t> mask(m_podkernel->getMesh()->getInternalCount(), 0);
        int i=0;
        for (const Cell &cell : m_podkernel->getMesh()->getCells()) {
            if (cell.isInterior()) {
                long id = cell.getId();
                mask[i] = m_filter[id];
                ++i;
            }
        }

        m_podkernel->getMesh()->getVTK().addData("filter", VTKFieldType::SCALAR, VTKLocation::CELL, mask);
        datastring.push_back("filter");
        m_podkernel->getMesh()->write();

        for (std::size_t is = 0; is < datastring.size(); ++is)
            m_podkernel->getMesh()->getVTK().removeData(datastring[is]);
    }
}

/**
 * Restore the POD results.
 */
void POD::restore()
{
    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    // Restore the POD info (no reconstruction info)
    log::cout()  << "Restore POD ... " << std::endl;

    std::string header = "pod";
    IBinaryArchive binaryReader(m_directory + m_name, dumpBlock);
    std::istream &dataStream = binaryReader.getStream();
    MeshType m_meshType_;
    utils::binary::read(dataStream, m_meshType_);
    setMeshType(m_meshType_);
    std::size_t nr_;
    utils::binary::read(dataStream, nr_);
    m_nModes = std::min(m_nModes, nr_);
    m_modes.resize(m_nModes);
    utils::binary::read(dataStream, m_nSnapshots);
    m_database.resize(m_nSnapshots);
    for (std::size_t i = 0; i < m_nSnapshots; i++) {
        utils::binary::read(dataStream, m_database[i].directory);
        utils::binary::read(dataStream, m_database[i].name);
    }
    utils::binary::read(dataStream,m_staticMesh);
    utils::binary::read(dataStream,m_memoryMode);
    binaryReader.close();

    // Restore mesh
    {
        pod::SnapshotFile podfile;
        podfile.directory = m_directory;
        podfile.name = m_name;
        m_podkernel->restoreMesh(podfile);
    }

    // Restore the modes, mean and mesh
    {
        m_filter.setStaticKernel(&m_podkernel->getMesh()->getCells());
        m_sensorMask.setStaticKernel(&m_podkernel->getMesh()->getCells());

        if (m_memoryMode == MEMORY_NORMAL) {
            for (std::size_t ir = 0; ir < m_nModes; ir++) {
                readMode(ir);
            }
        }
    }

    // Restore mean
    {
        readMode(m_nModes);
    }

    // Restore mask
    {
        IBinaryArchive binaryReader(m_directory + m_name + ".mask.data", dumpBlock);
        std::istream &dataStream = binaryReader.getStream();
        m_sensorMask.restore(dataStream);
        binaryReader.close();
    }

    // Lestore ID list
    fillListID(m_sensorMask);
    m_toUpdate = true;
}

/**
 * Write snapshots by dumping.
 *
 * \param[in] name Filename.
 * \param[in] snap Snapshot field.
 */
void POD::dumpField(const std::string &name, const pod::PODField &field) const
{
    log::cout()  << "Dumping snapshot ... " << std::endl;

    int dumpBlock = (m_nProcs > 1) ? m_rank : -1;

    // Dump field data
    {
        std::string header = name;
        OBinaryArchive binaryWriter(m_directory + name + ".data", ARCHIVE_VERSION, header, dumpBlock);
        std::ostream &dataStream = binaryWriter.getStream();
        m_filter.dump(dataStream);
        utils::binary::write(dataStream,std::size_t(m_nScalarFields));
        if (m_nScalarFields) {
            for (std::size_t i = 0; i < m_nScalarFields; i++)
                utils::binary::write(dataStream, m_nameScalarFields[i]);

            field.scalar->dump(dataStream);
        }

        utils::binary::write(dataStream,std::size_t(m_nVectorFields));
        if (m_nVectorFields) {
            for (std::size_t i = 0; i < m_nVectorFields; i++)
                utils::binary::write(dataStream, m_nameVectorFields[i]);

            field.vector->dump(dataStream);
        }
        binaryWriter.close();
    }

    // Dumping mesh
    {
        std::string header = name + "mesh";
        OBinaryArchive binaryWriter(m_directory + "/" + name+".mesh", ARCHIVE_VERSION, header, dumpBlock);
        field.mesh->dump(binaryWriter.getStream());
        binaryWriter.close();
    }

    if (m_writeMode == WriteMode::DEBUG) {
        std::vector<std::string > datastring;
        m_podkernel->getMesh()->getVTK().setDirectory(m_directory);
        field.mesh->getVTK().setName(name);
        std::vector<std::vector<double>> fields(m_nScalarFields, std::vector<double>(field.mesh->getInternalCount(), 0.0));
        std::vector<std::vector<std::array<double,3>>> fieldv(m_nVectorFields, std::vector<std::array<double,3>>(field.mesh->getInternalCount(), {{0.0, 0.0, 0.0}}));
        for (std::size_t isf = 0; isf < m_nScalarFields; isf++) {
            int i=0;
            for (const Cell &cell : field.mesh->getCells()) {
                if (cell.isInterior()) {
                    long id = cell.getId();
                    fields[isf][i] = field.scalar->at(id, isf);
                    ++i;
                }
            }
            field.mesh->getVTK().addData(m_nameScalarFields[isf], VTKFieldType::SCALAR, VTKLocation::CELL, fields[isf]);
            datastring.push_back(m_nameScalarFields[isf]);
        }

        for (std::size_t ivf = 0; ivf < m_nVectorFields; ivf++) {
            int i=0;
            for (const Cell &cell : field.mesh->getCells()) {
                if (cell.isInterior()) {
                    long id = cell.getId();
                    fieldv[ivf][i] = field.vector->at(id, ivf);
                    ++i;
                }
            }
            std::string vname= m_nameVectorFields[ivf][0].substr(0,m_nameVectorFields[ivf][0].size()-2);
            field.mesh->getVTK().addData(vname, VTKFieldType::VECTOR, VTKLocation::CELL, fieldv[ivf]);
            datastring.push_back(vname);
        }

        std::vector<uint8_t> mask(field.mesh->getInternalCount(), 0);
        int i=0;
        for (const Cell &cell : field.mesh->getCells()) {
            if (cell.isInterior()) {
                long id = cell.getId();
                mask[i] = m_filter[id];
                ++i;
            }
        }
        field.mesh->getVTK().addData("filter", VTKFieldType::SCALAR, VTKLocation::CELL, mask);
        datastring.push_back("filter");
        field.mesh->write();

        for (std::size_t is=0; is < datastring.size(); ++is)
            field.mesh->getVTK().removeData(datastring[is]);
    }
}

/**
 * Perform difference between a POD field and a POD mode.
 *
 * \param[in] a Field as POD object.
 * \param[in] b Mode as POD object.
 */
void POD::diff(pod::PODField& a, const pod::PODMode& b)
{
    std::size_t nsf = 0;
    if (m_nScalarFields){
        nsf = a.scalar->getFieldCount();
        if (nsf != b.scalar->getFieldCount())
            throw std::runtime_error ("POD: different field count");
    }

    std::size_t nvf = 0;
    if (m_nVectorFields){
        nvf = a.vector->getFieldCount();
        if (nvf != b.vector->getFieldCount())
            throw std::runtime_error ("POD: different field count");
    }

    //WARNING: it works only if nvf and nsf !=0 TODO
    for (long id : a.scalar->getKernel()->getIds()){
        double* datasa = a.scalar->data(id);
        double* datasb = b.scalar->data(id);
        for (std::size_t i = 0; i < nsf; i++){
            *datasa = *datasa - *datasb;
            datasa++;
            datasb++;
        }
        std::array<double,3>* datava = a.vector->data(id);
        std::array<double,3>* datavb = b.vector->data(id);
        for (std::size_t i = 0; i < nvf; i++)
            *datava = *datava - *datavb;
        datava++;
        datavb++;
    }
}

/**
 * Perform sum between a POD field and a POD mode.
 *
 * \param[in] a Field as POD object.
 * \param[in] b Mode as POD object.
 */
void POD::sum(pod::PODField& a, const pod::PODMode& b)
{
    std::size_t nsf = 0;
    if (m_nScalarFields){
        nsf = a.scalar->getFieldCount();
        if (nsf != b.scalar->getFieldCount())
            throw std::runtime_error ("POD: different field count");
    }

    std::size_t nvf = 0;
    if (m_nVectorFields){
        nvf = a.vector->getFieldCount();
        if (nvf != b.vector->getFieldCount())
            throw std::runtime_error ("POD: different field count");
    }

    for (long id : a.scalar->getKernel()->getIds()){
        double* datasa = a.scalar->data(id);
        double* datasb = b.scalar->data(id);
        for (std::size_t i = 0; i < nsf; i++){
            *datasa = *datasa + *datasb;
            datasa++;
            datasb++;
        }
        std::array<double,3>* datava = a.vector->data(id);
        std::array<double,3>* datavb = b.vector->data(id);
        for (std::size_t i = 0; i < nvf; i++)
            *datava = *datava + *datavb;
        datava++;
        datavb++;
    }
}

#if BITPIT_ENABLE_MPI
/**
 * Initializes the MPI communicator to be used for parallel communications.
 *
 * \param communicator is the communicator.
 */
void POD::initializeCommunicator(MPI_Comm communicator)
{
    // Communication can be set just once
    if (isCommunicatorSet())
        throw std::runtime_error ("POD communicator can be set just once");

    // The communicator has to be valid
    if (communicator == MPI_COMM_NULL)
        throw std::runtime_error ("POD communicator is not valid");

    // Create a copy of the user-specified communicator
    //
    // No library routine should use MPI_COMM_WORLD as the communicator;
    // instead, a duplicate of a user-specified communicator should always
    // be used.
    MPI_Comm_dup(communicator, &m_communicator);
}

/**
 * Returns the MPI communicator stored within LevelSetKernel.
 * @return MPI communicator.
 */
MPI_Comm POD::getCommunicator() const
{
    return m_communicator;
}

/**
 * Checks if the communicator to be used for parallel communications has
 * already been set.
 *
 * \result Returns true if the communicator has been set, false otherwise.
 */
bool POD::isCommunicatorSet() const
{

    return (getCommunicator() != MPI_COMM_NULL);
}

/**
 * Frees the MPI communicator associated to the patch
 */
void POD::freeCommunicator()
{
    if (!isCommunicatorSet())
        return;

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled)
        return;

    MPI_Comm_free(&m_communicator);
}
#endif

/**
 * Field reconstruction.
 *
 * \param[in,out] fields Original field used in the reconstruction. Only the values of \p targetCells are updated.
 * \param[in] mesh Pointer to snapshot field.
 * \param[in] scalarNames Map of names/ids of scalar fields in PiercedStorage.
 * \param[in] vectorNames Map of names/ids of vector fields in PiercedStorage.
 * \param[in] targetCells Pointer to list of cell to be updated with the reconstructed values (optional, default whole field).
 */
void POD::reconstructFields(PiercedStorage<double> &fields, const VolumeKernel *mesh,
        const std::map<std::string, std::size_t> & scalarNames,
        const std::map<std::array<std::string, 3>, std::array<std::size_t, 3>> & vectorNames,
        const std::unordered_set<long> * targetCells)
{
    std::vector<std::size_t> podscalarIds(scalarNames.size());
    std::vector<std::size_t> podvectorIds(vectorNames.size());

    std::vector<std::size_t> scalarIds(scalarNames.size());
    std::vector<std::array<std::size_t, 3>> vectorIds(vectorNames.size());

    std::size_t count = 0;
    for (auto val : scalarNames){
        std::string name = val.first;
        std::vector<std::string>::iterator found = find(m_nameScalarFields.begin(), m_nameScalarFields.end(), name);
        if (found == m_nameScalarFields.end())
            throw std::runtime_error ("POD: scalar field name not found in pod base");

        std::size_t ind = distance(m_nameScalarFields.begin(), found);
        podscalarIds[count] = ind;
        scalarIds[count] = val.second;
        count++;
    }
    count = 0;
    for (auto val : vectorNames){
        std::array<std::string, 3> name = val.first;
        std::vector<std::array<std::string,3>>::iterator found = find(m_nameVectorFields.begin(), m_nameVectorFields.end(), name);
        if (found == m_nameVectorFields.end())
            throw std::runtime_error ("POD: vector field name not found in pod base");

        std::size_t ind = distance(m_nameVectorFields.begin(), found);
        podvectorIds[count] = ind;
        vectorIds[count] = val.second;
        count++;
    }

    evalReconstructionCoeffs(fields, mesh, scalarIds, podscalarIds, vectorIds, podvectorIds);

    buildFields(fields, scalarIds, podscalarIds, vectorIds, podvectorIds, targetCells);
}

/**
 * Field reconstruction.
 *
 * \param[in,out] fields Original field used in the reconstruction. Only the values of \p targetCells are updated.
 * \param[in] mesh Pointer to snapshot field
 * \param[in] fields Map of names/ids of fields in PiercedStorage
 * \param[in] targetCells Pointer to list of cell to be updated with the reconstructed values (optional, default whole field).
 */
void POD::reconstructFields(PiercedStorage<double> &fields, const VolumeKernel *mesh,
        std::map<std::string, std::size_t> targetFields,
        const std::unordered_set<long> * targetCells)
{
    std::vector<std::size_t> podscalarIds;
    std::vector<std::size_t> podvectorIds;

    std::vector<std::size_t> scalarIds;
    std::vector<std::array<std::size_t, 3>> vectorIds;

    // Find scalar target fields
    std::size_t count = 0;
    for (std::string val : m_nameScalarFields) {
        std::map<std::string, std::size_t>::iterator found = targetFields.find(val);
        if (found != targetFields.end()) {
            std::size_t ind = found->second;
            podscalarIds.push_back(count);
            scalarIds.push_back(ind);
            targetFields.erase(found);
            count++;
        }
    }

    // Find vector target fields
    count = 0;
    for (std::array<std::string,3> valv : m_nameVectorFields) {
        std::size_t ic = 0;
        std::array<std::size_t,3> vectorarr;
        std::array<std::string,3> toerase;
        for (std::string val : valv) {
            std::map<std::string, std::size_t>::iterator found = targetFields.find(val);
            if (found != targetFields.end()) {
                std::size_t ind = found->second;
                vectorarr[ic] = ind;
                toerase[ic] = val;
                ic++;
            }

            if (ic == 3) {
                podvectorIds.push_back(count);
                vectorIds.push_back(vectorarr);
                for (std::size_t i = 0; i < 3; i++)
                    targetFields.erase(toerase[i]);
            }
        }
        count++;
    }

    if (targetFields.size() != 0) {
        std::string verror = " ";
        for (const auto &targetField : targetFields)
            verror += targetField.first + " ";

        throw std::runtime_error ("POD: fields not found in POD base: " + verror);
    }

    evalReconstructionCoeffs(fields, mesh, scalarIds, podscalarIds, vectorIds, podvectorIds);

    buildFields(fields, scalarIds, podscalarIds, vectorIds, podvectorIds, targetCells);
}

/**
 * Evaluation of POD reconstruction coefficients.
 *
 * \param[in] fields Original input field (optimad solver format).
 * \param[in] mesh Pointer to snapshot field.
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] scalarNames Names of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 * \param[in] vectorNames Names of vector fields in PiercedStorage.
 */
void POD::evalReconstructionCoeffs(PiercedStorage<double> &fields, const VolumeKernel *mesh,
        const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> & podscalarIds,
        const std::vector<std::array<std::size_t, 3>> & vectorIds, const std::vector<std::size_t> & podvectorIds)
{
    BITPIT_UNUSED(mesh);

    if (m_staticMesh)
        evalReconstructionCoeffsStaticMesh(fields, scalarIds, podscalarIds, vectorIds, podvectorIds);
    else
        throw std::runtime_error("Dynamic mesh is not supported yet.");

#if BITPIT_ENABLE_MPI
    for (std::size_t i = 0; i < m_nFields; i++)
        MPI_Allreduce(MPI_IN_PLACE, m_reconstructionCoeffs[i].data(), m_nModes, MPI_DOUBLE, MPI_SUM, m_communicator);
#endif
}

/**
 * Evaluation of POD reconstruction coefficients on static mesh.
 *
 * \param[in] fields Original input field (optimad solver format).
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] scalarNames Names of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 * \param[in] vectorNames Names of vector fields in PiercedStorage.
 */
void POD::evalReconstructionCoeffsStaticMesh(PiercedStorage<double> &fields,
        const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> &podscalarIds,
        const std::vector<std::array<std::size_t, 3>> &vectorIds, const std::vector<std::size_t> &podvectorIds)
{
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();

    // Difference field - m_mean
    for (long  id : m_listActiveIDs) {
        double *datag = fields.data(id);
        double *datams = m_mean.scalar->data(id);
        for (std::size_t i = 0; i < nsf; i++) {
            double *datagi = datag + scalarIds[i];
            double *datamsi = datams + podscalarIds[i];
            (*datagi) -= (*datamsi);
        }
        std::array<double,3>* datamv = m_mean.vector->data(id);
        for (std::size_t i = 0; i < nvf; i++) {
            std::array<double,3>* datamvi = datamv + podvectorIds[i];
            for (std::size_t j = 0; j < 3; j++) {
                double *datagi = datag + vectorIds[i][j];
                (*datagi) -= (*datamvi)[j];
            }
        }
    }

    m_reconstructionCoeffs.clear();
    m_reconstructionCoeffs.resize(m_nFields, std::vector<double>(m_nModes, 0.0));
    std::vector<std::vector<double> > rhs = m_reconstructionCoeffs;

    // Evaluate minimization matrices
    evalMinimizationMatrices();

    // Evaluate RHS
    for (std::size_t ir = 0; ir < m_nModes; ++ir) {
        if (m_memoryMode == POD::MEMORY_LIGHT)
            readMode(ir);

        std::vector<long>::iterator it;
        for (it = m_listActiveIDs.begin(); it != m_listActiveIDs.begin()+m_sizeInternal; it++) {
            long id = *it;
            double *datag = fields.data(id);
            if (nsf) {
                double *modes = m_modes[ir].scalar->data(id);
                for (std::size_t ifield = 0; ifield < nsf; ifield++) {
                    double *modesi = modes + podscalarIds[ifield];
                    double *datagi = datag + scalarIds[ifield];
                    rhs[ifield][ir] += (*datagi)*(*modesi)*m_podkernel->getMesh()->evalCellVolume(id);
                }
            }
            if (nvf) {
                std::array<double,3>* modev = m_modes[ir].vector->data(id);
                for (std::size_t ifield = 0; ifield < nvf; ifield++) {
                    std::array<double,3>* modevi = modev + podvectorIds[ifield];
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datagi = datag + vectorIds[ifield][j];
                        rhs[m_nScalarFields+ifield][ir] += ((*datagi)*(*modevi)[j])*m_podkernel->getMesh()->evalCellVolume(id);
                    }
                }
            }
        }
    }

#if BITPIT_ENABLE_MPI
    for (std::size_t i = 0; i < m_nFields; i++)
        MPI_Allreduce(MPI_IN_PLACE, rhs[i].data(), m_nModes, MPI_DOUBLE, MPI_SUM, m_communicator);
#endif

    // Solve
    solveMinimization(rhs);

    // Sum field and mean
    for (long id : m_listActiveIDs) {
        double *datag = fields.data(id);
        double *datams = m_mean.scalar->data(id);
        for (std::size_t i = 0; i < nsf; i++) {
            double *datagi = datag + scalarIds[i];
            double *datamsi = datams + podscalarIds[i];
            *datagi = *datagi + *datamsi;
        }
        std::array<double,3>* datamv = m_mean.vector->data(id);
        for (std::size_t i = 0; i < nvf; i++) {
            std::array<double,3>* datamvi = datamv + podvectorIds[i];
            for (std::size_t j = 0; j < 3; j++) {
                double *datagi = datag + vectorIds[i][j];
                *datagi = *datagi + (*datamvi)[j];
            }
        }
    }
}

/**
 * Reconstruct a field through POD. Use the m_reconstructionCoeffs already set.
 * Interface for hybrid code [optimad solver]. The reconstruction is performed
 * only on selected cells.
 *
 * \param[in,out] fields Input and resulting reconstructed field.
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 * \param[in] podscalarIds Ids of scalar fields in POD modes.
 * \param[in] podvectorIds Ids of vector fields in POD modes.
 * \param[in] targetCells Pointer to list of target cells for reconstruction (optional, default whole field).
 */
void POD::buildFields(PiercedStorage<double> &fields,
        const std::vector<std::size_t> &scalarIds, const std::vector<std::size_t> &podscalarIds,
        const std::vector<std::array<std::size_t, 3>> &vectorIds, const std::vector<std::size_t> &podvectorIds,
        const std::unordered_set<long> *targetCells)
{
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();
    if (nsf == 0 && nvf == 0)
        return;

    std::vector<long> targetCellsStorage;
    if (!targetCells) {
        std::unordered_set<long> targetCellsStorage;
        for (const Cell &cell : m_podkernel->getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    for (const long id : *targetCells) {
        double *recon = fields.data(id);
        for (std::size_t ifs = 0; ifs < m_nScalarFields; ifs++) {
            double *reconsi = recon + scalarIds[ifs];
            (*reconsi) = 0.;
        }
        for (std::size_t ifv = 0; ifv < m_nVectorFields; ifv++) {
            for (std::size_t j = 0; j < 3; j++) {
                double *reconvi = recon + vectorIds[ifv][j];
                (*reconvi) = 0.;
            }
        }
    }


    for (std::size_t ir = 0; ir < m_nModes; ir++) {
        if (m_memoryMode == POD::MEMORY_LIGHT)
            readMode(ir);

        for (const long id : *targetCells) {
            double *modes = m_modes[ir].scalar->data(id);
            double *recon = fields.data(id);
            for (std::size_t ifs = 0; ifs < m_nScalarFields; ifs++) {
                double *modesi = modes + podscalarIds[ifs];
                double *reconsi = recon + scalarIds[ifs];
                (*reconsi) += (*modesi) * m_reconstructionCoeffs[ifs][ir];
            }

            std::array<double,3> *modev = m_modes[ir].vector->data(id);
            for (std::size_t ifv = 0; ifv < m_nVectorFields; ifv++) {
                std::array<double,3>* modevi = modev + podvectorIds[ifv];
                for (std::size_t j = 0; j < 3; j++) {
                    double *reconvi = recon + vectorIds[ifv][j];
                    (*reconvi) += (*modevi)[j] * m_reconstructionCoeffs[m_nScalarFields + ifv][ir];
                }
            }
        }

        if (m_memoryMode == POD::MEMORY_LIGHT)
            m_modes[ir].clear();
    }

    // Sum field and mean
    for (const long id : *targetCells) {
        double *datag = fields.data(id);
        double *datams = m_mean.scalar->data(id);
        for (std::size_t i = 0; i < nsf; i++) {
            double *datagi = datag + scalarIds[i];
            double *datamsi = datams + podscalarIds[i];
            (*datagi) += (*datamsi);
        }

        std::array<double,3>* datamv = m_mean.vector->data(id);
        for (std::size_t i = 0; i < nvf; i++) {
            std::array<double,3>* datamvi = datamv + podvectorIds[i];
            for (std::size_t j = 0; j < 3; j++) {
                double *datagi = datag + vectorIds[i][j];
                (*datagi) += (*datamvi)[j];
            }
        }
    }
}

}
