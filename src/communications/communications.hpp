/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#if BITPIT_ENABLE_MPI==1

#ifndef __BITPIT_COMMUNICATIONS_HPP__
#define __BITPIT_COMMUNICATIONS_HPP__

#include <mpi.h>
#include <vector>
#include <unordered_map>

#include "bitpit_containers.hpp"

#include "communications_buffers.hpp"

namespace bitpit {

class DataCommunicator
{

public:
    DataCommunicator(MPI_Comm communicator);

    void finalize();

    void setTag(int tag);
    int getTag() const;

    void clearAllSends();
    void clearAllRecvs();

    void clearSend(int rank);
    void clearRecv(int rank);

    void setSend(int rank, long length = 0);
    void setRecv(int rank, long length = 0);

    void resizeSend(int rank, long resize);
    void resizeRecv(int rank, long resize);

    void setRecvsContinuous(bool enabled);
    bool areRecvsContinuous();

    void discoverRecvs();
    void discoverRecvs(int discoverTag);

    int getSendCount();
    int getRecvCount();

    const std::vector<int> getSendRanks() const;
    const std::vector<int> getRecvRanks() const;

    SendBuffer & getSendBuffer(int rank);
    RecvBuffer & getRecvBuffer(int rank);

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
    static int DEFAULT_TAG;

    MPI_Comm m_communicator;
    int m_rank;
    int m_tag;
    bool m_recvsContinuous;

    std::vector<int> m_recvRanks;
    std::unordered_map<int, int> m_recvIds;
    std::vector<MPI_Request> m_recvRequests;
    std::vector<RecvBuffer> m_recvBuffers;

    std::vector<int> m_sendRanks;
    std::unordered_map<int, int> m_sendIds;
    std::vector<MPI_Request> m_sendRequests;
    std::vector<SendBuffer> m_sendBuffers;

};

}

# endif

# endif
