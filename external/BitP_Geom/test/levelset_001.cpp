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

void Generate2DSurfMesh(
    bitpit::SurfUnstructured                   &mesh
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
const double                            _PI_ = 3.14159265358979;
const double                            R = 1.0;
const long                              N = 32;
double                                  theta;
double                                  dtheta = 2. * _PI_/((double) N);

// Counters
long                                    i;

// ========================================================================== //
// GENERATE MESH                                                              //
// ========================================================================== //
{
    // Local variables ------------------------------------------------------ //
    std::array<double, 3>                    point;
    std::vector<long>                        connect(2, bitpit::Element::NULL_ID);

    // Create vertex list --------------------------------------------------- //
    point[2] = 0.0;
    for (i = 0; i < N; ++i) {
        theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        mesh.addVertex(point);
    } //next i

    // Create simplex list -------------------------------------------------- //
    for (i = 0; i < N; ++i) {
        connect[0] = i;
        connect[1] = (i+1) % N;
        mesh.addCell(bitpit::ElementInfo::LINE, true, connect);
    } //next i
}

return; }

int main( int argc, char *argv[]){

#if NOMPI==0
    MPI::Init(argc,argv);
#endif


    int                    dimensions(2) ;

    // Input geometry
    bitpit::SurfUnstructured STL(0,1,dimensions);

    std::cout << " - Loading dgf geometry" << std::endl;

//    STL.importSTL("./data/naca0012.dgf") 

//    STL.deleteCoincidentVertices() ;
//    STL.buildAdjacencies() ;
    Generate2DSurfMesh( STL ) ;

    STL.setName("geometry_001") ;
    STL.write() ;

    std::cout << "n. vertex: " << STL.getVertexCount() << std::endl;
    std::cout << "n. simplex: " << STL.getCellCount() << std::endl;



    // create cartesian mesh around geometry 
    std::cout << " - Setting mesh" << std::endl;
    std::array<double,3>     meshMin, meshMax, delta ;
    std::array<int,3>        nc = {{64, 64, 0}} ;

    STL.getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);

    // Compute level set  in narrow band
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int                                      elapsed_seconds;

    bitpit::LevelSetCartesian       LSP(mesh);
    bitpit::LevelSetSegmentation    geometry(0,&STL);

    start = std::chrono::system_clock::now();
    LSP.compute( &geometry ) ;
    end = std::chrono::system_clock::now();

    elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    std::cout << " - Exporting data" << std::endl;
    mesh.update() ;
    std::vector<double> LS(mesh.getCellCount() ) ;
    std::vector<double>::iterator it = LS.begin() ;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        *it = LSP.getLS(id) ;
        ++it ;
    };

    mesh.addData("ls", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, LS) ;
    mesh.setName("levelset_001") ;
    mesh.write() ;

    std::cout << " - Exported data" << std::endl;



#if NOMPI==0
    MPI::Finalize();
#endif

    return 0;

};


