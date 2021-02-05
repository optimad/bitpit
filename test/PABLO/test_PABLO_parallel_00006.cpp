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

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

#include <mpi.h>

#include <vector>

using namespace bitpit;


/*!
* Driver for handling user balance load balance.
*/
class UserDataLB : public bitpit::DataLBInterface<UserDataLB> {

public:
    static constexpr int DATA_SIZE = 300;

    struct Data {
        int a;
        double b[DATA_SIZE];
    };

    std::vector<Data> &data;
    std::vector<Data> &ghostdata;

    size_t fixedSize() const
    {
        return sizeof(int) + DATA_SIZE*sizeof(double);
    }

    void move(const uint32_t from, const uint32_t to)
    {
        data[to] = data[from];
    }

    template<class Buffer>
    void gather(Buffer & buff, const uint32_t e)
    {
        buff << data[e].a;
        for(int i=0; i<DATA_SIZE; i++)
        buff << data[e].b[i];
    }

    template<class Buffer>
    void scatter(Buffer & buff, const uint32_t e)
    {
        buff >> data[e].a;
        for(int i=0; i<DATA_SIZE; i++)
        buff >> data[e].b[i];
    }

    void assign(uint32_t stride, uint32_t length)
    {
        for(uint32_t i=0; i<length; i++)
        data[i] = data[i+stride];
    }
    void resize(uint32_t newSize)
    {
        data.resize(newSize);
    }
    void resizeGhost(uint32_t newSize)
    {
        ghostdata.resize(newSize);
    }
    void shrink() {}

    UserDataLB(std::vector<Data>& data_, std::vector<Data>& ghostdata_)
        : data(data_), ghostdata(ghostdata_)
    {}

    ~UserDataLB(){}
};

/*!
* Subtest 001
*
* Testing load balance with data on a two-dimensional tree.
*/
int subtest_001(int rank)
{
    /**<Instantation of a 3D pablo uniform object.*/
    bitpit::PabloUniform tree(2);

    /**<Set 2:1 balance for the octree.*/
    int idx = 0;
    tree.setBalance(idx,true);

    /** Set Periodic boundary conditions */
    tree.setPeriodic(0);
    tree.setPeriodic(2);

    /**<Compute the connectivity and write the octree.*/
    tree.computeConnectivity();

    /**<Refine globally.*/
    tree.adaptGlobalRefine();
    tree.adaptGlobalRefine();

    tree.setMarker(10,1);
    tree.setMarker(15,-1);

    tree.adapt();
    assert( !tree.checkToAdapt() );

    /**<First load balance without data.*/
    uint8_t levels = 4;
    tree.loadBalance(levels);

    /**<Adaption.*/
    tree.adapt();
    assert( !tree.checkToAdapt() );

    /**<Second load balance with data.*/
    {
      std::vector<UserDataLB::Data> oct_data(tree.getNumOctants()), ghost_data(tree.getNumGhosts());
      UserDataLB data_lb(oct_data,ghost_data);
      tree.loadBalance(data_lb, levels);
    }

    return 0;
}

/*!
* Subtest 001
*
* Testing load balance with data on a three-dimensional tree.
*/
int subtest_002(int rank)
{
    /**<Instantation of a 3D pablo uniform object.*/
    bitpit::PabloUniform tree(2);

    /**<Set 2:1 balance for the octree.*/
    int idx = 0;
    tree.setBalance(idx,true);

    /** Set Periodic boundary conditions */
    tree.setPeriodic(0);
    tree.setPeriodic(2);
    tree.setPeriodic(4);

    /**<Compute the connectivity and write the octree.*/
    tree.computeConnectivity();

    /**<Refine globally.*/
    tree.adaptGlobalRefine();
    tree.adaptGlobalRefine();
    tree.adaptGlobalRefine();

    tree.setMarker(10,1);
    tree.setMarker(15,-1);

    tree.adapt();
    assert( !tree.checkToAdapt() );

    /**<First load balance without data.*/
    uint8_t levels = 4;
    tree.loadBalance(levels);

    /**<Adaption.*/
    if( rank==0 ) {
      tree.setMarker(10, 1);
    } else if( rank==2 ) {
      tree.setMarker(15, -1);
    }
    tree.adapt();
    assert( !tree.checkToAdapt() );

    /**<Second load balance with data.*/
    {
      std::vector<UserDataLB::Data> oct_data(tree.getNumOctants()), ghost_data(tree.getNumGhosts());
      UserDataLB data_lb(oct_data,ghost_data);
      tree.loadBalance(data_lb, levels);
    }

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initialize the logger
    bitpit::log::manager().initialize(bitpit::log::SEPARATE, false, nProcs, rank);
    bitpit::log::cout() << fileVerbosity(bitpit::log::NORMAL);
    bitpit::log::cout() << consoleVerbosity(bitpit::log::QUIET);

    // Run the subtests
    log::cout() << "Testing load balance with data" << std::endl;

    int status;
    try {
        status = subtest_001(rank);
        if (status != 0) {
            return status;
        }

        status = subtest_002(rank);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
