#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_surfunstructured.hpp"
#include "bitpit_voloctree.hpp"
#include "bitpit_levelset.hpp"
#include "bitpit_LA.hpp"
#include "bitpit_operators.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_RBF.hpp"
#include "bitpit_common.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iomanip>

using namespace std;
using namespace bitpit;

/*!
	\example RBF_example_00001.cpp

	\brief Generation of RBF object that represents a signed distance function

    bitpit for levelset RBF generation
    This file explains how to use the script linking the Levelset and RBF capabilities of bitpit to train an RBF
    collection that represents the signed distance field of STL-defined objects.

    Generator script usage:
    The `/examples/RBF_example_00001` executable is the one to use to generate RBFs from `.stl` files. The syntax to
    use it is `./RBF_example_00001 <path_to_stl_file> <name_of_stl> <parameters_file>.param`.
    1. <path_to_stl_file> is the path, absolute or relative, that contains the desired STL file.
    2. <name_of_stl> is the name of the STL file. The name must not include the `.stl` extension.
    3. <parameters_file>.param contains as many columns as parameters that must be specified. The first line contains the
       name of the parameters, and the second line contains the numerical values. All numbers and words on a same line
       must be separated by ` `. The parameters are `nb_subdivision`, `nb_adapations`, `radius_ratio`, `base_function`,
       `mesh_range`, `max_num_threads`, `tolerance`. The parameters have default values, so they don't all have to be in
       the parameter file.
        -`nb_subdivision` dictates the background grid initial refinement. Having a higher number of `nb_subdivisions`
        leads to a more precise but more costly RBF.
        -`nb_adaptations` dictates the number of adaptive refinement cycles.
         The RBF will have a higher concentration of nodes around high error zones, which can decrease the training time,
         the evaluation time and the file size.
        -`radius_ratio` dictates the radius that each node has an influence over. This values should be kept low enough so
        that sharp geometry features are not lost.
        -`base_function` dictates the base function type with an integer. The value `2` corresponds to a linear function
        and produces the lowest error in many cases.
        -`mesh_range` dictates how far outside of the STL bounding box the background mesh should extend. Often, a `0.1`
        value, which represents a 10% increase on each side, is a good choice.
        -`tolerance` dictates the expected tolerance when solving the system with all the nodes.

    The executable generates two files.
    1. RBF_<name_of_stl>.output: this file contains the location of every node and its associated weight, function and
    support radius.
    2. RBF_<name_of_stl>.vtu: this visualization file contains the levelset field of the considered object.

    <b>To run</b>: ./RBF_example_00001 \n
    <b>To run (with arguments)</b>: ./RBF_example_00001 ./ RBF_example_00001 RBF_example_00001.param \n
*/

void
parse_parameters(std::map<std::string, std::vector<double>> &map,
                 std::string                                 file,
                 const std::string                           delimiter)
{
    // fill a pair, first being a vector of vector name and the second being the
    // vector of vector associated with the vector name.
    std::ifstream myfile(file);
    // open the file.
    if (myfile.is_open())
    {
        std::string              line;
        std::vector<std::string> column_names;
        std::vector<double>      line_of_data;
        unsigned int             line_count = 0;

        while (std::getline(myfile, line))
        {
            // read the line and clean the resulting vector.
            std::vector<std::string> list_of_words_base;

            std::string s = line;
            size_t pos = 0;
            std::string token;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                list_of_words_base.push_back(token);
                s.erase(0, pos + delimiter.length());
            }
            std::vector<std::string> list_of_words_clean;
            for (std::size_t i = 0; i < list_of_words_base.size(); ++i)
            {
                if (list_of_words_base[i] != "")
                {
                    list_of_words_clean.push_back(list_of_words_base[i]);
                }
            }
            // check if the line is contained words or numbers.
            if (line_count != 0)
            {
                line_of_data.resize(list_of_words_clean.size());
                for (std::size_t i = 0; i < line_of_data.size(); i++)
                {
                    line_of_data[i] = std::stod(list_of_words_clean[i]);
                }
                for (std::size_t i = 0; i < line_of_data.size(); ++i)
                {
                    map[column_names[i]].push_back(line_of_data[i]);
                }
            }
            else
            {
                // the line contains words, we assume these are the columns names.
                column_names = list_of_words_clean;
                for (std::size_t i = 0; i < list_of_words_clean.size(); ++i)
                {
                    std::vector<double> base_vector;
                    map[column_names[i]] = base_vector;
                }
            }
            ++line_count;
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    // We add here the default values
    std::vector<std::string> names = {"nb_subdivision",
                                      "nb_adaptions",
                                      "radius_ratio",
                                      "base_function",
                                      "mesh_range",
                                      "tolerance",
                                      "scaling"};
    std::vector<double>      values = {16,
                                       0,
                                       1,
                                       2,
                                       0.2,
                                       1e-8,
                                       1.0};
    for (std::size_t i = 0; i < names.size(); i++)
    {
        if (map.find("f") == map.end())
            map[names[i]].push_back(values[i]);
    }
}

void run(std::string filename,
         std::string data_path,
         std::string parameter_file) {
    constexpr int dimensions(3);

    // Parsing the parameters
    std::map<std::string, std::vector<double>> parameters;
    parse_parameters(parameters, parameter_file," ");
    int nb_subdivision      = static_cast<int>(parameters["nb_subdivision"][0]);
    int nb_adaptions        = static_cast<int>(parameters["nb_adaptions"][0]);
    double radius_ratio     =                  parameters["radius_ratio"][0];
    int base_function       = static_cast<int>(parameters["base_function"][0]);
    double mesh_range       =                  parameters["mesh_range"][0];
    double TOL              =                  parameters["tolerance"][0];
    double scaling          =                  parameters["scaling"][0];

    std::vector<std::string> timers_name;
    std::vector<double> timers_values;

    timers_name.push_back("load_geometry");
    double time_start = MPI_Wtime();

    //LEVELSET PART
    //Input geometry
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> STL0(new bitpit::SurfUnstructured(dimensions - 1, MPI_COMM_NULL));
#else
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (dimensions - 1) );
#endif
    bitpit::log::cout() << " - Loading stl geometry" << std::endl;
    // Make sure that the STL format is in binary (not ASCII)
    try {
        STL0->importSTL(data_path + filename + ".stl", true);
    }
    catch (const std::bad_alloc) {
        STL0->importSTL(data_path + filename + ".stl", false);
    }
    STL0->deleteCoincidentVertices();
    STL0->initializeAdjacencies();
    STL0->getVTK().setName("levelset");
    std::array<double, dimensions> center{};
    STL0->scale(scaling,scaling,scaling,center);
    bitpit::log::cout() << "n. vertex: " << STL0->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL0->getCellCount() << std::endl;
    // Create initial octree mesh for levelset
    bitpit::log::cout() << " - Setting mesh" << std::endl;
    std::array<double, dimensions> stlMin, stlMax, meshMin, meshMax, delta;
    double h(0.), dh;
    STL0->getBoundingBox(stlMin, stlMax);
    bitpit::log::cout() << " Bounding box minima " << std::endl;
    bitpit::log::cout() << stlMin << std::endl;
    bitpit::log::cout() << " Bounding box maxima " << std::endl;
    bitpit::log::cout() << stlMax << std::endl;
    bitpit::log::cout() << " Bounding box center " << std::endl;
    bitpit::log::cout() << 0.5*(stlMax+stlMin) << std::endl;

    delta = stlMax - stlMin;
    meshMin = stlMin - mesh_range * delta;
    meshMax = stlMax + mesh_range * delta;
    for (int i = 0; i < dimensions; ++i) {
        h = std::max(h, meshMax[i] - meshMin[i]);
    }
    dh = h / nb_subdivision;
#if BITPIT_ENABLE_MPI
    bitpit::VolOctree mesh(dimensions, meshMin, h, dh, MPI_COMM_NULL);
#else
    bitpit::VolOctree mesh(dimensions, meshMin, h, dh);
#endif
    mesh.initializeAdjacencies();
    mesh.initializeInterfaces();
    mesh.update();
    mesh.getVTK().setName("RBF_" + filename);
#if BITPIT_ENABLE_MPI
    mesh.setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif

    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("compute_levelset");
    time_start = MPI_Wtime();


    // Set levelset configuration
    bitpit::LevelSet levelset;
    levelset.setPropagateSign(true);
    levelset.setSizeNarrowBand(sqrt(3.0) * h);
    levelset.setMesh(&mesh);

    int id0 = levelset.addObject(std::move(STL0), 0);
    const bitpit::LevelSetObject &object0 = levelset.getObject(id0);
    levelset.getObject(id0).enableVTKOutput(bitpit::LevelSetWriteField::VALUE);

    // Write levelset information
    mesh.write();
    bitpit::log::cout() << "Computed levelset within the narrow band... " << std::endl;


    // Adaptative Refinement
    std::vector<bitpit::adaption::Info> adaptionData_levelset;
    for (int r = 0; r < nb_adaptions; ++r) {
        for (auto &cell: mesh.getCells()) {
            long cellId = cell.getId();
            if (std::abs(object0.getValue(cellId)) < mesh.evalCellSize(cellId))
                mesh.markCellForRefinement(cellId);
        }
        adaptionData_levelset = mesh.update(true);
        levelset.update(adaptionData_levelset);
        mesh.write();
    }
    unsigned long nP_total = mesh.getCellCount();


    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("add_nodes_to_RBF");
    time_start = MPI_Wtime();


    // RBF PART
    bitpit::RBFBasisFunction basisFunction;
    switch (base_function) {
        case 0:
            basisFunction = bitpit::RBFBasisFunction::CUSTOM;
            break;
        case 1:
            basisFunction = bitpit::RBFBasisFunction::WENDLANDC2;
            break;
        case 2:
            basisFunction = bitpit::RBFBasisFunction::LINEAR;
            break;
        case 3:
            basisFunction = bitpit::RBFBasisFunction::GAUSS90;
            break;
        case 4:
            basisFunction = bitpit::RBFBasisFunction::GAUSS95;
            break;
        case 5:
            basisFunction = bitpit::RBFBasisFunction::GAUSS99;
            break;
        case 6:
            basisFunction = bitpit::RBFBasisFunction::C1C0;
            break;
        case 7:
            basisFunction = bitpit::RBFBasisFunction::C2C0;
            break;
        case 8:
            basisFunction = bitpit::RBFBasisFunction::C0C1;
            break;
        case 9:
            basisFunction = bitpit::RBFBasisFunction::C1C1;
            break;
        case 10:
            basisFunction = bitpit::RBFBasisFunction::C2C1;
            break;
        case 11:
            basisFunction = bitpit::RBFBasisFunction::C0C2;
            break;
        case 12:
            basisFunction = bitpit::RBFBasisFunction::C1C2;
            break;
        case 13:
            basisFunction = bitpit::RBFBasisFunction::C2C2;
            break;
        case 14:
            basisFunction = bitpit::RBFBasisFunction::COSINUS;
            break;
        default:
            basisFunction = bitpit::RBFBasisFunction::LINEAR;
            break;
    }
    std::vector<double> values;
    std::vector<double> weights;
    std::vector<double> radii;
    std::vector<std::array<double, dimensions>> nodes;
    values.resize(nP_total);
    weights.resize(nP_total);
    nodes.resize(nP_total);
    radii.resize(nP_total);

    bitpit::RBF RBFObject;
    RBFObject.setMode(bitpit::RBFMode::PARAM);
    RBFObject.setFunction(basisFunction);

    bitpit::log::cout() << "Adding nodes to the RBF" << std::endl;
    for (size_t it_RBF = 0; it_RBF < nP_total; it_RBF++) {
        nodes[it_RBF] = mesh.evalCellCentroid(it_RBF);
        values[it_RBF] = levelset.getObject(id0).getValue(it_RBF);
        RBFObject.addNode(nodes[it_RBF]);
        radii[it_RBF] = mesh.evalCellSize(it_RBF) * radius_ratio;
    }
    RBFObject.setSupportRadius(radii);

    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("fill_matrix");
    time_start = MPI_Wtime();

    // Training the RBF
    bitpit::log::cout() << "Training RBFObject" << std::endl;
    // Initializing the matrix and vector
    int nNZ = nP_total*(1+(2*(int)ceil(radius_ratio)^dimensions));
#if BITPIT_ENABLE_MPI==1
    SparseMatrix A(MPI_COMM_WORLD, true, nP_total, nP_total, nNZ);
#else
    SparseMatrix A(nP_total, nP_total, nNZ);
#endif

    std::vector<long> rowPattern;
    std::vector<double> rowValues;
    // Compute the values for A
    bitpit::log::cout() << "Filling the matrix" << std::endl;
    for (unsigned long i = 0; i < nP_total; i++) {
        for (unsigned long j = 0; j < nP_total; j++) {
            double v_ij = RBFObject.evalBasisPair(i,j);
            if (abs(v_ij) > TOL) {
                rowPattern.push_back(j);
                rowValues.push_back(v_ij);
            }
        }
        A.addRow(rowPattern, rowValues);
        rowPattern.clear();
        rowValues.clear();
    }
    A.assembly();

    SystemSolver system;
    system.assembly(A);

    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("fill_RHS");
    time_start = MPI_Wtime();

    double *b = system.getRHSRawPtr();
    // Compute the values for b
    bitpit::log::cout() << "Filling the RHS" << std::endl;
    for (unsigned long i = 0; i < nP_total; ++i) {
        b[i] = values[i];
    }
    system.restoreRHSRawPtr(b);

    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("solve_system");
    time_start = MPI_Wtime();

    bitpit::log::cout() << "Solving the system" << std::endl;

    double *initialSolution = system.getSolutionRawPtr();
    for (unsigned long i = 0; i < nP_total; ++i) {
        initialSolution[i] = values[i];
    }
    system.restoreSolutionRawPtr(initialSolution);

    system.solve();
    const double *x = system.getSolutionRawReadPtr();
    bitpit::log::cout() << "Solved system" << std::endl;

    timers_values.push_back(MPI_Wtime() - time_start);
    timers_name.push_back("add_weights_to_RBFObject");
    time_start = MPI_Wtime();

    bitpit::log::cout() << "Adding weights to RBFObject" << std::endl;
    for (unsigned long i = 0; i < nP_total; i++) {
        weights[i] = x[i];
    }
    system.restoreSolutionRawPtr(initialSolution);
    RBFObject.addData(weights);
    bitpit::log::cout() << "Added weights to RBFObject" << std::endl;
    bitpit::log::cout() << "Finished RBF training" << std::endl;

    std::size_t nMaxCellVertices = mesh.getCell(0).getVertexIds().size();

    timers_values.push_back(MPI_Wtime()-time_start);
    timers_name.push_back("output_RBF_and_analytical");
    time_start = MPI_Wtime();

    // RBF Output
    bitpit::log::cout() << "Outputting" << std::endl;
    std::vector<double> display_values;
    display_values.resize(nP_total);

    BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA dimensions>, nMaxCellVertices,
                            ReferenceElementInfo::MAX_ELEM_VERTICES)
    const PatchKernel::CellConstRange &cellWriteRange = mesh.getVTKCellWriteRange();
    for (const Cell &cell:cellWriteRange){
        const unsigned long global_id = cell.getId();
        mesh.getCellVertexCoordinates(global_id, vertexCoordinates);
        std::array<double, dimensions> point = cell.evalCentroid(vertexCoordinates);
        std::vector<double> temp_disp = RBFObject.evalRBF(point);
        display_values[global_id] = temp_disp[0]; //Only the first field is used, since there is only one
    }
    mesh.getVTK().addData<double>("RBF", VTKFieldType::SCALAR, VTKLocation::CELL, display_values);
    mesh.write();

    timers_values.push_back(MPI_Wtime()-time_start);
    timers_name.push_back("output_RBF_to_file");
    time_start = MPI_Wtime();

    //Outputting the combined RBF to a txt file
    ofstream fw("RBF_" + filename + ".output", std::ofstream::ate);
    if (fw.is_open()) {
        fw << "support_radius basis_function node_x node_y ";
        if (dimensions == 3)
        {
            fw << "node_z ";
        }
        fw << "weight \n";
        for (unsigned long line = 0; line < nP_total; line++) {
            // Set precision
            int PRECISION = 20; //number of decimals
            std::ostringstream streamObj_support_radius;
            std::ostringstream streamObj_x;
            std::ostringstream streamObj_y;
            std::ostringstream streamObj_z;
            std::ostringstream streamObj_weight;
            streamObj_support_radius << std::fixed;
            streamObj_x              << std::fixed;
            streamObj_y              << std::fixed;
            streamObj_z              << std::fixed;
            streamObj_weight         << std::fixed;
            streamObj_support_radius << std::setprecision(PRECISION);
            streamObj_x              << std::setprecision(PRECISION);
            streamObj_y              << std::setprecision(PRECISION);
            streamObj_z              << std::setprecision(PRECISION);
            streamObj_weight         << std::setprecision(PRECISION);
            streamObj_support_radius << radii[line];
            streamObj_x << nodes[line][0];
            streamObj_y << nodes[line][1];
            streamObj_z << nodes[line][2];
            streamObj_weight << weights[line];
            std::string strObj_support_radius = streamObj_support_radius.str();
            std::string strObj_x              = streamObj_x.str();
            std::string strObj_y              = streamObj_y.str();
            std::string strObj_z              = streamObj_z.str();
            std::string strObj_weight         = streamObj_weight.str();
            fw << strObj_support_radius << " " << base_function << " "<<
                  strObj_x << " " << strObj_y << " ";
            if (dimensions == 3)
            {
                fw << strObj_z << " ";
            }
            fw << strObj_weight << " ";
            if (line < nP_total -1)
                fw << "\n";
        }
        fw.close();
    } else bitpit::log::cout() << "Problem with opening file" << std::endl;
    bitpit::log::cout() << "Finished outputting" << std::endl;

    timers_values.push_back(MPI_Wtime()-time_start);

    int nb_timers = timers_name.size();
    bitpit::log::cout()<< std::endl << "Timers" << std::endl;
    for (int t = 0; t < nb_timers; t++)
    {
        bitpit::log::cout() << timers_name.at(t) << ":" << timers_values.at(t) << std::endl;
    }
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    int nProcs = 1;
    int rank   = 0;

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nProcs>1)
    {
        bitpit::log::cout() << "nProcs > 1 isn't supported at the moment" << std::endl;
        exit(1);
    }
#endif

    // Arguments
    std::vector<std::string> argList;
    for(int i=0;i<argc;i++)
        argList.emplace_back(argv[i]);

    // Default values
    std::string data_path      = "../test/integration_tests/levelset/data/";
    std::string filename       = "cube";
    std::string parameter_file = "RBF_example_00001.param";
    if (argc > 3)
    {
        data_path      = argList[1];
        filename       = argList[2];
        parameter_file = argList[3];
    }

    // Initialize the logger
	log::manager().initialize(log::MODE_COMBINE, true, nProcs, rank);
	log::cout() << log::fileVerbosity(log::INFO);
	log::cout() << log::disableConsole();

	// Run the example
    try {
        run(filename,
            data_path,
            parameter_file);
    }
    catch (const std::exception &exception) {
    log::cout() << exception.what();
    exit(1);
    }

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
    return 0;
}
