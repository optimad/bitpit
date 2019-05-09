/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
    static const int TAG_AUTO = -1;

    DataCommunicator(MPI_Comm communicator);
    ~DataCommunicator();

    const MPI_Comm & getCommunicator() const;

    void finalize(bool synchronous = false);

    void setTag(int exchangeTag);
    void setTags(int exchangeTag, int discoverTag, int notificationTag);
    void setExchangeTag(int tag);
    void setDiscoverTag(int tag);
    void setNotificationTag(int tag);
    int getTag() const;
    int getExchangeTag() const;
    int getDiscoverTag() const;
    int getNotificationTag() const;

    void clearAllSends(bool synchronous = false);
    void clearAllRecvs(bool synchronous = false);

    void clearSend(int rank);
    void clearRecv(int rank);

    void setSend(int rank, long length = 0);
    void setRecv(int rank, long length = 0);

    void resizeSend(int rank, long resize);
    void resizeRecv(int rank, long resize);

    void setRecvsContinuous(bool enabled);
    bool areRecvsContinuous();

    void discoverSends();
    void discoverRecvs();

    int getSendCount();
    int getRecvCount();

    const std::vector<int> & getSendRanks() const;
    const std::vector<int> & getRecvRanks() const;

    SendBuffer & getSendBuffer(int rank);
    RecvBuffer & getRecvBuffer(int rank);

    void startSend(int dstRank);
    void startAllSends();

    void startRecv(int srcRank);
    void startAllRecvs();

    int waitAnySend(const std::vector<int> &blackList = std::vector<int>());
    void waitSend(int rank);
    void waitAllSends();

    int waitAnyRecv(const std::vector<int> &blackList = std::vector<int>());
    void waitRecv(int rank);
    void waitAllRecvs();

    bool isSendActive(int rank);
    bool isRecvActive(int rank);

    void cancelSend(int rank);
    void cancelRecv(int rank);

    void cancelAllSends(bool synchronous = false);
    void cancelAllRecvs(bool synchronous = false);

private:
    MPI_Comm m_communicator;
    int m_rank;
    int m_exchangeTag;
    int m_discoverTag;
    int m_notificationTag;
    bool m_customExchangeTag;
    bool m_customDiscoverTag;
    bool m_customNotificationTag;
    bool m_recvsContinuous;

    std::vector<int> m_recvRanks;
    std::unordered_map<int, int> m_recvIds;
    std::vector<MPI_Request> m_recvRequests;
    std::vector<RecvBuffer> m_recvBuffers;

    std::vector<int> m_sendRanks;
    std::unordered_map<int, int> m_sendIds;
    std::vector<MPI_Request> m_sendRequests;
    std::vector<SendBuffer> m_sendBuffers;

    void _startSend(int dstRank);
    void _startRecv(int srcRank);

    MPI_Datatype getChunkDataType(int chunkSize) const;

};

}

# endif

# endif
