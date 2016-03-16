/*!
 *	\date			10/jul/2014
 *	\authors		Alessandro Alaia
 *	\authors		Haysam Telib
 *	\authors		Edoardo Lombardi
 *	\version		0.1
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Level Set Class Demos
 */

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

//Standard Template Library
# include <ctime>
# include <chrono>

// bitpit
# include <bitpit.hpp>
# include "BitP_Geom_LEVELSET.hpp"
//# include "User_LSData_LB.hpp"
//# include "User_VLSData_LB.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

/*!Demo for 3D level set of complex geometries on a Pablo octree mesh.
*/
int main( int argc, char *argv[]){

#if ENABLE_MPI==1
    MPI::Init(argc,argv);
#endif

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    uint8_t                 dimensions(3);


    // Input geometry
    bitpit::SurfUnstructured    STL(0);

    std::cout << " - Loading stl geometry" << std::endl;

    STL.importSTL("./data/cube1.stl", true);

    STL.deleteCoincidentVertices() ;
    STL.buildAdjacencies() ;

    STL.setName("geometry_002") ;
    STL.write() ;

    std::cout << "n. vertex: " << STL.getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL.getCellCount() << std::endl;



    // ========================================================================== //
    // CREATE MESH                                                                //
    // ========================================================================== //
    
    std::cout << " - Setting mesh" << std::endl;
    std::array<double,3>    meshMin, meshMax, delta ;
    double                  h(0), dh ;

    STL.getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    for( int i=0; i<3; ++i){
        h = max( h, meshMax[i]-meshMin[i] ) ;
    };

    dh = h / 16. ;
    bitpit::VolOctree    mesh(1, dimensions, meshMin, h, dh );
    mesh.update() ;


    // COMPUTE LEVEL SET in NARROW BAND
    std::chrono::time_point<std::chrono::system_clock>    start, end;
    int                                         elapsed_init, elapsed_refi(0);

    bitpit::LevelSetOctree          LSP(mesh);
    bitpit::LevelSetSegmentation    geometry(0,&STL);
    std::vector<bitpit::Adaption::Info> mapper ;
    std::vector<double>             LS ;
    std::vector<long>               SG ;
    std::vector<double>::iterator   itLS ;
    std::vector<long>::iterator     itSG ;


    mesh.addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS) ;
    mesh.addData("sg", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, SG) ;
    mesh.setName("levelset_004") ;
    mesh.setCounter() ;

    start = std::chrono::system_clock::now();
    LSP.compute( &geometry );
    end = std::chrono::system_clock::now();

    elapsed_init = chrono::duration_cast<chrono::milliseconds>(end-start).count();


    // Export level set ------------------------------------------------------- //
    if (mesh.getRank() == 0) cout << " - Exporting data" << endl;

    LS.resize(mesh.getCellCount() ) ;
    SG.resize(mesh.getCellCount() ) ;
    itLS = LS.begin() ;
    itSG = SG.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        *itLS = LSP.getLS(id) ;
        *itSG = geometry.getSupportSimplex(id) ;
        ++itLS ;
        ++itSG ;
    };

    mesh.write() ;

    //Refinement

    for( int i=0; i<1; ++i){

        std::cout << "rfinement loop " << i << std::endl ;
        std::cout << "rfinement loop " << i << std::endl ;
        std::cout << "rfinement loop " << i << std::endl ;

        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            if( std::abs(LSP.getLS(id)) < 100. ){
                mesh.markCellForRefinement(id) ;
            }
        }

        mapper = mesh.update(true) ;
        start = std::chrono::system_clock::now();
        LSP.update( &geometry, mapper) ;
        end = std::chrono::system_clock::now();

        elapsed_refi += chrono::duration_cast<chrono::milliseconds>(end-start).count();

        LS.resize(mesh.getCellCount() ) ;
        SG.resize(mesh.getCellCount() ) ;
        itLS = LS.begin() ;
        itSG = SG.begin() ;
        for( auto & cell : mesh.getCells() ){
            const long &id = cell.getId() ;
            *itLS = LSP.getLS(id) ;
            *itSG = geometry.getSupportSimplex(id) ;
            ++itLS ;
            ++itSG ;
        };
        mesh.write() ;
    }


    cout << "elapsed time initialization " << elapsed_init << " ms" << endl;
    cout << "elapsed time refinement     " << elapsed_refi << " ms" << endl;

#if ENABLE_MPI==1
    MPI::Finalize();
#endif

    return 0;

};


