/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2023 OPTIMAD engineering Srl
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

/**
 * \example POD_example_00005.cpp
 * 
 * \brief POD basis computation using voloctree.
 * This example 
 * - computes the POD basis starting from a database of simulations defined on the same mesh.
 * - Reads and writes the first snapshot and the reconstruction as VTK file.
 * - Projects the first snapshot on POD vector space generate by the modes.
 * - Reconstructs the first snapshots with the BuildFieldsWithCoeff function.
 * - Compares snapshot and reconstruction.
 * - Builds a new pod object restoring the modes dumped in the computation of the first POD object.
 * - Compares the reconstruction of the first POD object with the reconstruction obtained with the second.
 * 
 * <b>To run</b>: ./POD_example_00005 \n
 */ 

#include <array>
#include <vector>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"
#include "pod_voloctree.hpp"
#include "pod_kernel.hpp"
#include "piercedStorage.hpp"

using namespace bitpit;

/*
 * Build a PODField object with difference between two PODField objects.
 *
 * \param[in] field1, minuend PODField.
 * \param[in] field2, subtrahend PODField.
 * \param[in] listActIds, list of active ids of the mesh of field1 and field2.
 * \param[out] errorField, PODField with the difference between field1 and field2.
 */
pod::PODField buildErrorField (pod::PODField &field1, pod::PODField &field2, const std::unordered_set<long> listActIds)
{
    pod::PODField errorField;
    errorField.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(2, &field1.mesh->getCells()));
    errorField.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(1, &field1.mesh->getCells()));
    errorField.scalar->fill(0.0);
    errorField.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
    errorField.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &field1.mesh->getCells()));
    errorField.mask->fill(0.0);
    errorField.setMesh(field1.mesh);
    for (long id : listActIds) {
        for (std::size_t isf = 0; isf < 2; ++isf) {
            double field1s = field1.scalar->at(id, isf);
            double field2s = field2.scalar->at(id, isf);
            errorField.scalar->at(id, isf) += field1s-field2s;
        }
        std::array<double,3> field1sv = field1.vector->at(id,0);
        std::array<double,3> field2sv = field2.vector->at(id,0);
        errorField.vector->at(id, 0) += field1sv-field2sv;
        errorField.mask->at(id) = 1;
    }
    return errorField;
}

/*
 * Build a PODField object with the relative difference between two PODField objects.
 *
 * \param[in] field1, minuend PODField.
 * \param[in] field2, subtrahend PODField.
 * \param[in] pod, POD object over which the two fields have been defined.
 * \param[in] norm_type, L2 or Linf
 * \param[out] errorField, PODField with the difference between field1 and field2.
 */
pod::PODField buildRelErrorField (pod::PODField &field1, pod::PODField &field2, POD &pod, std::string norm_type )
{
    pod::PODField errorField;
    errorField.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(2, &field1.mesh->getCells()));
    errorField.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(1, &field1.mesh->getCells()));
    errorField.scalar->fill(0.0);
    errorField.vector->fill(std::array<double, 3>{{0.0, 0.0, 0.0}});
    errorField.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &field1.mesh->getCells()));
    errorField.mask->fill(0.0);
    errorField.setMesh(field1.mesh);
    std::vector<double> vec;
    if (norm_type == "L2") {
        vec = pod.fieldsl2norm(field1);
    }
    else if (norm_type == "Linf") {
        vec = pod.fieldsMax(field1);
    }
    const std::unordered_set<long> & listActIds = pod.getListActiveIDs();
    for (long id : listActIds) {
        for (std::size_t isf = 0; isf < 2; ++isf) {
            double field1s = field1.scalar->at(id, isf);
            double field2s = field2.scalar->at(id, isf);
            errorField.scalar->at(id, isf) += (field1s-field2s)/vec[isf];
        }
        std::array<double,3> field1sv = field1.vector->at(id,0);
        std::array<double,3> field2sv = field2.vector->at(id,0);
        errorField.vector->at(id, 0) += (field1sv-field2sv)/vec[2];
        errorField.mask->at(id) = 1;
    }
    return errorField;
}


/*
 * Print the L2 norm of each field of a PODField object.
 *
 * \param[in] field, PODField object.
 * \param[in] pod, POD object defined on the same mesh.
 * \param[in] field_name, string of the field name to display.
 */
void printL2norm (pod::PODField &field, POD &pod, std::string field_name)
{
    std::vector<double> vecL2 = pod.fieldsl2norm(field);
    std::vector<std::string> scalarNames = pod.getScalarNames();
    std::vector<std::array<std::string,3>> vectorNames= pod.getVectorNames();
    int N = scalarNames.size();
    for (int i=0; i<N; i++) {
        std::cout << "L2 norm of " << field_name << " " << scalarNames[i] << " is "<< vecL2[i] << std::endl;
    }
    int M = vectorNames.size();
    for (int i=N; i<N+M; i++) {
        std::cout << "L2 norm of " << field_name << " " << vectorNames[i-N][0].substr(0,vectorNames[i-N][0].size()-2) << " is "<< vecL2[i] << std::endl;
    }
}

/*
 * Print the L infinity norm of each field of a PODField object.
 *
 * \param[in] field, PODField object.
 * \param[in] pod, POD object defined on the same mesh.
 * \param[in] field_name, string of the field name to display.
 */
void printLinfnorm (pod::PODField &field, POD &pod, std::string field_name)
{
    std::vector<double> vecLinf = pod.fieldsMax(field);
    std::vector<std::string> scalarNames = pod.getScalarNames();
    std::vector<std::array<std::string,3>> vectorNames= pod.getVectorNames();
    int N = scalarNames.size();
    for (int i=0; i<N; i++) {
        std::cout << "L infinity norm of " << field_name << " " << scalarNames[i] << " is "<< vecLinf[i] << std::endl;
    }
    int M = vectorNames.size();
    for (int i=N; i<N+M; i++) {
        std::cout << "L infinity norm of " << field_name << " " << vectorNames[i-N][0].substr(0,vectorNames[i-N][0].size()-2) << " is "<< vecLinf[i] << std::endl;
    }
}

/*
 * Print a matrix.
 *
 * \param[in] mat, matrix of size M rows times N columns.
 */
void printMat (std::vector < std::vector<double>> mat)
{
    std::cout << "mat = " << std::endl;
    size_t M = mat.size();
    size_t N = mat[0].size();
    for (size_t i=0; i<M; i++) {
        for (size_t j=0; j<N; j++) {
            if (j == 0) {
                std::cout << "[ "<< mat[i][j] ;
            }
            else if (j==(N-1)) {
                std::cout << " , "  << mat[i][j] << " ]" << std::endl;
            }
            else {
                std::cout << " , " << mat[i][j] ;
            }
            if (N==1) {
                std::cout  << " ]" << std::endl;
            }
        }
    }
}

/**
 * Run the example.
 */ 
void run()
{
    /**<Create POD object.*/
    POD pod;

    /**<Add snapshots to database.*/
    for (int i=0; i<6; i++) {
        pod.addSnapshot("./data", "test_set2."+std::to_string(i));
    }

    /**<Set POD.*/
    pod.setMeshType(POD::MeshType::VOLOCTREE);
    pod.setStaticMesh(true);
    pod.setErrorMode(POD::ErrorMode::SINGLE);
    pod.setWriteMode(POD::WriteMode::DEBUG);
    pod.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    pod.setEnergyLevel(99.00);
    pod.setUseMean(false);
    pod.setDirectory("pod");
    pod.setName("pod.test.solver");

    /**<Compute the POD basis.*/
    pod.compute();

    /**<Read and write the first solution as VTK file */
    const pod::SnapshotFile snap0 ("./data", "test_set2.0");
    pod::PODField solution0;
    pod.readSnapshot(snap0, solution0);
    pod.write(solution0,"solution0");

    /**<Reconconstruct the first snapshot using the POD coefficients */
    std::vector < std::vector<double>> recon_coeff = pod.projectField(solution0);
    std::cout << "POD coefficients of the first snapshot " << std::endl;
    printMat(recon_coeff);
    pod::PODField recon0;
    pod.buildFieldsWithCoeff(recon_coeff, recon0);
    pod.write(recon0,"reconstruction0");

    /**<Compare the solution and reconstruction */
    const std::unordered_set<long> & listActIds = pod.getListActiveIDs();
    pod::PODField error0 = buildErrorField(solution0, recon0, listActIds);
    pod::PODField error0relL2 = buildRelErrorField(solution0, recon0, pod, "L2");
    pod::PODField error0relLinf = buildRelErrorField(solution0, recon0, pod, "Linf");
    printL2norm(solution0,pod,"solution");
    printLinfnorm(solution0,pod,"solution");
    printL2norm(error0,pod,"reconstruction error");
    printLinfnorm(error0,pod,"reconstruction error");
    printL2norm(error0relL2,pod,"relative error");
    printLinfnorm(error0relLinf,pod,"relative error");

    /**<Restore modes in a new pod object */
    POD pod2;
    // Set pod2.
    pod2.setWriteMode(POD::WriteMode::DEBUG);
    pod2.setUseMean(false);
    pod2.setDirectory("pod");
    pod2.setName("pod.test.solver");
    // restore modes.
    pod2.restore();
    std::cout << "the number of modes of pod 2 is " << pod2.getModeCount() << std::endl;
    std::cout << "the number of modes of pod is " << pod.getModeCount() << std::endl;
    // get modes from pod2.
    const std::vector<pod::PODMode> &vec_pod_modes = pod2.getModes();
    // get list of active ids from pod2
    const std::unordered_set<long> & listActIds_2 = pod2.getListActiveIDs();
    std::cout << "Number of active cells for pod2 object = ";
    std::cout << listActIds_2.size() << std::endl;
    // compare with list of active ids from pod1
    std::cout << "Number of active cells for pod object = ";
    std::cout << listActIds.size() << std::endl;

    /**<Compare the reconstruction of pod with the reconstruction of pod2 */
    // get reconstruction coefficients used for the computation of recon0
    double coeff11 = recon_coeff[0][0];
    double coeff12 = recon_coeff[0][1];
    double coeff21 = recon_coeff[1][0];
    double coeff22 = recon_coeff[1][1];
    double coeff31 = recon_coeff[2][0];
    double coeff32 = recon_coeff[2][1];
    double scalar1_diff = 0;
    double scalar2_diff = 0;
    double vector1_diff = 0;
    for (long id : listActIds_2) {
        double mode1 = vec_pod_modes[0].scalar->at(id, 0);
        double mode2 = vec_pod_modes[1].scalar->at(id, 0);
        scalar1_diff += coeff11*mode1+coeff12*mode2-recon0.scalar->at(id, 0);
        mode1 = vec_pod_modes[0].scalar->at(id, 1);
        mode2 = vec_pod_modes[1].scalar->at(id, 1);
        scalar2_diff += coeff21*mode1+coeff22*mode2-recon0.scalar->at(id, 1);
        std::array<double,3> mode1v = vec_pod_modes[0].vector->at(id, 0);
        std::array<double,3> mode2v = vec_pod_modes[1].vector->at(id, 0);
        std::array<double,3> diffv = recon0.vector->at(id, 0);
        diffv -= coeff31*mode1v+coeff32*mode2v;
        vector1_diff += std::sqrt( dotProduct((diffv),(diffv)));
    }
    std::cout << " " << std::endl;
    std::cout << "The absolute error between the reconstruction and the linear combination of the restored modes is " << std::endl;
    std::cout << scalar1_diff << " for the first scalar field, " << std::endl;
    std::cout << scalar2_diff << " for the second scalar field, "<< std::endl;
    std::cout << vector1_diff << " for the first vector field. "<< std::endl;
}

/**
 * Main program.
 */
int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#endif    

    /** Run the example */
    try {
        run();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif

}
