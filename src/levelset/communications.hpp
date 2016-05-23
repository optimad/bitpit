// ========================================================================== //
//                         - G.L.O.R.I.A. -                                   //
//                                                                            //
// Boundary routines for gloria solver                                        //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Andrea Iob
// Version  : v1.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

#if BITPIT_ENABLE_MPI==1

#ifndef __LEVELSET_COMMUNICATIONS_HPP__
#define __LEVELSET_COMMUNICATIONS_HPP__

#include <mpi.h>
#include <vector>
#include <unordered_map>

#include <bitpit_containers.hpp>
#include <bitpit_IO.hpp>
#include <bitpit_patchkernel.hpp>

class DataCommunicator
{

public:
	DataCommunicator(MPI_Comm communicator);

	void finalize();

	void setTag(int tag);

	void clearAllSends();
	void clearAllRecvs();

	void clearSend(int rank);
	void clearRecv(int rank);

	void setSend(int rank, long length = 0);
	void setRecv(int rank, long length = 0);

	void growSend(int rank, long growth);
	void growRecv(int rank, long growth);

	void setRecvsContinuous(bool enabled);
	bool areRecvsContinuous();

	int getSendCount();
	int getRecvCount();

	bitpit::OBinaryStream & getSendBuffer(int rank);
	bitpit::IBinaryStream & getRecvBuffer(int rank);

	void startSend(int dstRank);
	void startAllSends();

	void startRecv(int srcRank);
	void startAllRecvs();

	int waitAnySend();
	void waitSend(int rank);
	void waitAllSends();

	int waitAnyRecv();
	void waitRecv(int rank);
	void waitAllRecvs();

	void cancelSend(int rank);
	void cancelRecv(int rank);

	void cancelAllSends();
	void cancelAllRecvs();

private:
	MPI_Comm m_communicator;
	int m_rank;
	int m_tag;
	bool m_recvsContinuous;

    std::vector<int> m_recvRanks;
    std::unordered_map<int, int> m_recvIds;
    std::vector<MPI_Request> m_recvRequests;
    std::vector<bitpit::IBinaryStream> m_recvBuffers;

    std::vector<int> m_sendRanks;
    std::unordered_map<int, int> m_sendIds;
    std::vector<MPI_Request> m_sendRequests;
    std::vector<bitpit::OBinaryStream> m_sendBuffers;

};

# endif

# endif
