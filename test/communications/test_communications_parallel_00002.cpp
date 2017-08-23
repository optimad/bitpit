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

#include "bitpit_IO.hpp"
#include "bitpit_communications.hpp"

using namespace bitpit;

/*!
 * Test for communications bigger than 2Gb.
 *
 * This test is meant to transfer a data packet bigger than 2Gb and check if
 * this type of communication is handled correctly. However, the side of the
 * exchanged data packets is limited by the constant MAX_MEMORY, which defines
 * the maximum memory available to the program, expressed in Mb. By defualt the
 * maximum available memory is limited to 10Mb. To really test communications
 * bigger than 2Gb the maximum avialble memory should be increased at least
 * to 2049 Mb.
 *
 */
int main(int argc, char *argv[]) {

    const double MAX_MEMORY = 10;

    MPI_Init(&argc, &argv);

    int nProcs;
    int    rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::GLOBAL);
    log::cout() << "Testing communications" << std::endl;

    DataCommunicator dataCommunicator(MPI_COMM_WORLD);

    // Data info
    typedef short DataType;

    const double SEND_SIZE_FROM_0_TO_1 = std::min(6.1 * 1024 * 1024 * 1024, MAX_MEMORY * 1024 * 1024);
    const double SEND_SIZE_FROM_N_TO_0 = std::min(256.0, MAX_MEMORY * 1024 * 1024);

    const DataType DEFAULT_VALUE = 21;
    const DataType LAST_VALUE    = 42;

    // Create the sends
    log::cout() << "Creating sends" << std::endl;

    int dstRank;
    if (rank == 0) {
        dstRank = 1;
    } else {
        dstRank = 0;
    }

    std::size_t nValues;
    if (rank == 0) {
        nValues = SEND_SIZE_FROM_0_TO_1 / sizeof(DataType);
    } else {
        nValues = SEND_SIZE_FROM_N_TO_0 / sizeof(DataType);
    }

    std::size_t dataSize = nValues * sizeof(DataType);

    log::cout() << "Initializing a send to " << dstRank << " of " << dataSize << " bytes" << std::endl;
    log::cout() << "  Number of values that will be send... " << nValues << std::endl;
    log::cout() << "  Data size that will be send... " << dataSize / (1024. * 1024) << " Mb " << std::endl;
    dataCommunicator.setSend(dstRank, dataSize);

    // Discorver the receives
    log::cout() << "Discovering receives" << std::endl;

    dataCommunicator.discoverRecvs();

    for (int rank : dataCommunicator.getRecvRanks()) {
        RecvBuffer &recvBuffer = dataCommunicator.getRecvBuffer(rank);
        std::size_t dataSize = recvBuffer.getSize();
        log::cout() << "  Discovered a receive from " << rank << " of " << dataSize << " bytes" << std::endl;
    }

    // Start receives
    log::cout() << "Starting the receives" << std::endl;

    dataCommunicator.startAllRecvs();

    // Start sending data
    log::cout() << "Starting data send" << std::endl;
    for (int rank : dataCommunicator.getSendRanks()) {
        log::cout() << "  Sending data to " << rank << std::endl;

        // Filling the buffer
        log::cout() << "    Filling the buffer" << std::endl;
        SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
        std::size_t nSendVaues = sendBuffer.getSize() / sizeof(DataType);
        for (std::size_t n = 0; n < nSendVaues - 1; ++n) {
            sendBuffer << DEFAULT_VALUE;
        }
        sendBuffer << LAST_VALUE;

        // Start sending data
        dataCommunicator.startSend(rank);
    }

    // Receive data
    log::cout() << "Waiting for data to arrive" << std::endl;

    int nCompletedRecvs = 0;
    while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
        int rank = dataCommunicator.waitAnyRecv();
        log::cout() << "Receiving data from " << rank << std::endl;

        RecvBuffer &recvBuffer = dataCommunicator.getRecvBuffer(rank);
        std::size_t nRecvValues = recvBuffer.getSize() / sizeof(DataType);
        log::cout() << "  Number of values that will be received... " << nRecvValues << std::endl;

        DataType value = 0;
        for (std::size_t n = 0; n < nRecvValues; ++n) {
            recvBuffer >> value;
        }

        log::cout() << " Last recevied value... " << (int) value <<std::endl;
        if (value != LAST_VALUE) {
            log::cout() << "Wrong value." << std::endl;
            log::cout() << "   Current data value : " << value << std::endl;
            log::cout() << "   Expected data value: " << LAST_VALUE << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 2);
        }

        ++nCompletedRecvs;
    }

    // Wait all sends
    dataCommunicator.startAllSends();

    // Finalize MPI
    MPI_Finalize();
}
