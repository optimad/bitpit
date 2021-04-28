#include <bitpit_surfunstructured.hpp>
#include <bitpit_patchkernel.hpp>
#include <bitpit_IO.hpp>
#include <exception>

//test function: will need dump of patch (patch.dump collection and 2 sub parts bxxxx) to perform restore
// and the partitioning map in cellRanks.dat file. launched with 2 processes (MPI)
void test() {

    bitpit::SurfUnstructured * patch = new bitpit::SurfUnstructured(2,3);
    MPI_Comm m_communicator;
    MPI_Comm_dup(MPI_COMM_WORLD, &m_communicator);
    int m_rank;
    MPI_Comm_rank(m_communicator, &m_rank);

    patch->setCommunicator(m_communicator); //without calling it first, the restore fails.

    //restore the patch from 2 proc dump collection file patch.dump.
    std::string filenameX = ("./data/patch");
    bitpit::IBinaryArchive binaryReader(filenameX, "dump", m_rank);
    patch->restore(binaryReader.getStream());
    binaryReader.close();

    //check the status of adjacencies for the first time.
    std::cout<<"bitpit adj build strategy :"<<m_rank<<"  "<<int(patch->getAdjacenciesBuildStrategy())<<std::endl;
    std::cout<<"bitpit adj dirty status   :"<<m_rank<<"  "<<int(patch->areAdjacenciesDirty())<<std::endl;

    //initializing Adjacencies
    patch->initializeAdjacencies();

    //recheck the status of adjacencies after initialization.
    std::cout<<"Recheck->bitpit adj build strategy :"<<m_rank<<"  "<<int(patch->getAdjacenciesBuildStrategy())<<std::endl;
    std::cout<<"Recheck->bitpit adj dirty status   :"<<m_rank<<"  "<<int(patch->areAdjacenciesDirty())<<std::endl;

    //write the restored patch.
    patch->write("00009_master");

    //read the partition map from file -> 0 rank will retain useful data
    std::unordered_map<long,int> partmap;
    long id;
    int rank;
    if(patch->getRank() == 0){

        std::ifstream input("./data/cellRanks.dat");
        while(!input.eof()){
            input>>id>>rank;
            partmap[id] = rank;

        }
        input.close();
    }

    //partition the patch -> the MPIWaitAny error is triggered here.
    patch->partition(partmap, false, true);

    // write the distributed patch.
    patch->write("00009_distributed");

    // delete the structure and exit.
    delete patch;
}

//main
int	main( int argc, char *argv[] ) {


#if BITPIT_ENABLE_MPI
    MPI_Init(&argc, &argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif
        try{
            /**<calling test*/
            test() ;
        }
        catch(std::exception & e){
            std::cout<<"test exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif

    return 0;
}
