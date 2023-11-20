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
 * \example POD_example_00004.cpp
 * 
 * \brief POD basis computation using voloctree.
 *
 * This example computes the POD basis starting from a database of simulations
 * defined on the same mesh and reconstructs the first mode as a PODField object
 * using the buildFieldWithCoeff function.
 * Following the computation, both the POD Field with the mode and POD Mode object
 * corresponding to the first mode are written on file in the example folder.
 * 
 * <b>To run</b>: ./POD_example_00004 \n
 */ 

#include <array>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"

using namespace bitpit;

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

    /* Remark: the reconstruction of the first mode is equal to the first mode
     * only if the mean is not used in the POD algorithm
     * if the mean is used, the mean has to be subtracted to the reconstruction
     * to get the first mode */

    /**<Reconstruc the first mode. */
    std::cout << "the number of modes is = " << pod.getModeCount() << std::endl;
    std::size_t N_modes = pod.getModeCount();
    std::vector<std::string> names = pod.getScalarNames();
    std::vector<std::array<std::string,3>> namev = pod.getVectorNames();
    std::size_t N_sfields = names.size();
    std::size_t N_vfiedls = namev.size();
    std::size_t N_fields = N_sfields+N_vfiedls;
    /* set up of the coefficient matrix
     * each column contains the coefficient of a specific mode
     * each row contains the coefficient of a specific field, either scalar or vector.*/
    std::vector < std::vector<double>> coeff_mat;
    coeff_mat.clear();
    coeff_mat.resize(N_fields, std::vector<double>(N_modes, 0.0));
    for (std::size_t i = 0; i < N_fields; i++) {
        coeff_mat[i][0] = 1;
    }
    pod::PODField mode1_recon;
    pod.buildFieldsWithCoeff(coeff_mat, mode1_recon);
    printL2norm(mode1_recon, pod, "mode 1");
    std::vector < std::vector<double>> test_mat = pod.projectField(mode1_recon);
    std::cout << "the coefficient matrix of the projection on the first mode is " << std::endl;
    printMat(test_mat);

    /* <Write the first mode as VTK file */
    pod.write(mode1_recon, "mode1_recon");
    pod.write(0, "mode1");
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
