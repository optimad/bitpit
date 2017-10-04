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

#include <mpi.h>

#include "bitpit_IO.hpp"
#include "bitpit_communications.hpp"

using namespace bitpit;

/*!
 * Auxiliary function to evaluate the number of values to send
 */
int getSendCount(int dstRank)
{
	return (dstRank + 1);
}

/*!
 * Test for basic parallel communications.
 *
 * After setting the sends, the DataCommunicator is used to automatically
 * discorver and set the receives. Finally some data is communicated among
 * the processors.
 *
 * \param rank is the rank of the process
 * \param nProcs is the number of processes
 */
int subtest_001(int rank, int nProcs)
{
	DataCommunicator dataCommunicator(MPI_COMM_WORLD);

	// Create data to send
	std::vector<double> data(nProcs + 1);
	for (int i = 0; i < nProcs; ++i) {
		data[i] = i;
	}

	// Create the sends
	log::cout() << "Creating sends" << std::endl;

	for (int i = 0; i < nProcs; ++i) {
		int dstRank  = i;
		int dataSize = getSendCount(dstRank) * sizeof(data[0]);

		log::cout() << "Initializing a send to " << dstRank << " with size " << dataSize << std::endl;
		dataCommunicator.setSend(dstRank, dataSize);
	}

	// Discorver the receives
	log::cout() << "Discovering receives" << std::endl;

	dataCommunicator.discoverRecvs();

	for (int i = 0; i < nProcs; ++i) {
		RecvBuffer &recvBuffer = dataCommunicator.getRecvBuffer(i);
		int dataSize = recvBuffer.getSize();
		log::cout() << "Discovered a receive from " << i << " with size " << dataSize << std::endl;

		int expectedDataSize = getSendCount(rank) * sizeof(data[0]);
		if (dataSize != expectedDataSize) {
			log::cout() << "Wrong data size." << std::endl;
			log::cout() << "   Current data size : " << dataSize << std::endl;
			log::cout() << "   Expected data size: " << expectedDataSize << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	}

	// Start receives
	dataCommunicator.startAllRecvs();

	// Start sending data
	for (int i = 0; i < nProcs; ++i) {
		SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(i);
		int nSendVaues = sendBuffer.getSize() / sizeof(data[0]);
		for (int n = 0; n < nSendVaues; ++n) {
			sendBuffer << data[n];
		}

		dataCommunicator.startSend(i);
	}

	// Receive data
	int nCompletedRecvs = 0;
	while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
		int rank = dataCommunicator.waitAnyRecv();
		log::cout() << "Receiving data from " << rank << std::endl;

		RecvBuffer &recvBuffer = dataCommunicator.getRecvBuffer(rank);
		int nRecvValues = recvBuffer.getSize() / sizeof(double);
		for (int n = 0; n < nRecvValues; ++n) {
			double value;
			recvBuffer >> value;

			if (value != data[n]) {
				log::cout() << "Wrong data value." << std::endl;
				log::cout() << "   Current data value : " << value << std::endl;
				log::cout() << "   Expected data value: " << data[n] << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 2);
			}
		}

		++nCompletedRecvs;
	}

	// Wait all sends
	dataCommunicator.startAllSends();

	// Done
	return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	// Initialize the logger
	int nProcs;
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setVisibility(log::GLOBAL);

	// Run the subtests
	log::cout() << "Testing basic parallel communications" << std::endl;

	int status;
	try {
		status = subtest_001(rank, nProcs);
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		log::cout() << exception.what();
		exit(1);
	}

	MPI_Finalize();
}
