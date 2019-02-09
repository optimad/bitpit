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

#if BITPIT_ENABLE_MPI==1

#include <unistd.h>

#include "bitpit_IO.hpp"

#include "communications.hpp"
#include "communications_tags.hpp"

namespace bitpit {

/*!
    \class DataCommunicator
    \ingroup communications

    \brief The DataCommunicator class provides the infrastructure needed to
    exchange data among processors.
*/

/*!
    Creates a new communicator for data exchange.
*/
DataCommunicator::DataCommunicator(MPI_Comm communicator)
    : m_communicator(communicator), m_rank(-1),
    m_recvsContinuous(false)
{
    // Get MPI information
    MPI_Comm_rank(m_communicator, &m_rank);

    // Set a tag
    setTags(TAG_AUTO, TAG_AUTO, TAG_AUTO);
}

/*!
    Destructor.
*/
DataCommunicator::~DataCommunicator()
{
    if (!m_customExchangeTag) {
        if (m_rank == 0) {
            communications::tags().trash(m_exchangeTag);
        }
    }

    if (!m_customDiscoverTag) {
        if (m_rank == 0) {
            communications::tags().trash(m_discoverTag);
        }
    }

    if (!m_customNotificationTag) {
        if (m_rank == 0) {
            communications::tags().trash(m_notificationTag);
        }
    }

    if (!m_customExchangeTag || !m_customDiscoverTag || !m_customNotificationTag) {
        MPI_Barrier(m_communicator);
    }
}

/*!
	Gets the MPI communicator

	\return The MPI communicator.
*/
const MPI_Comm & DataCommunicator::getCommunicator() const
{
	return m_communicator;
}

/*!
    Finalizes the communicator

    \param synchronous if this parameter is set to true, the function will
    cancel all requests of the current rank and than will wait until all
    remote processes that exchange data with the current rank cancel their
    requests as well. This guarantees that, at the end of the function,
    both the current rank and the remote processes involved in a data
    exchange with the current rank have canceled their requests.
*/
void DataCommunicator::finalize(bool synchronous)
{
    // Cancels all sends
    cancelAllSends(synchronous);

    // Cancels all receives
    cancelAllRecvs(synchronous);
}

/*!
    Sets the tag to be used for the data exchange.

    By default, a unique tag for data exchange is generated in the constructor.
    However, using this function, it is possible to assign a custom tag.

    \param exchangeTag is the custom tag to be used for data exchange
*/
void DataCommunicator::setTag(int exchangeTag)
{
    setExchangeTag(exchangeTag);
}

/*!
    Sets the tags to be used for the communications.

    By default, unique tags are generated in the constructor. However, using
    this function, it is possible to assign custom tags.

    \param exchangeTag is the custom tag to be used for data exchange
    \param discoverTag is the custom tag to be used for data size discover
    \param notificationTag is the custom tag to be used for notifications
*/
void DataCommunicator::setTags(int exchangeTag, int discoverTag, int notificationTag)
{
    setExchangeTag(exchangeTag);
    setDiscoverTag(discoverTag);
    setNotificationTag(notificationTag);
}

/*!
    Sets the tag to be used for data exchange

    By default, a unique tag is generated in the constructor. However, using
    this function, it is possible to assign a custom tag.

    \param tag is the custom tag to be used for data exchange
*/
void DataCommunicator::setExchangeTag(int tag)
{
    m_customExchangeTag = (tag != TAG_AUTO);
    if (m_customExchangeTag) {
        m_exchangeTag = tag;
    } else {
        if (m_rank == 0) {
            m_exchangeTag = communications::tags().generate();
        }

        MPI_Bcast(&m_exchangeTag, 1, MPI_INT, 0, m_communicator);
    }
}

/*!
    Sets the tag to be used for data size discover

    By default, a unique tag is generated in the constructor. However, using
    this function, it is possible to assign a custom tag.

    \param tag is the custom tag to be used for data size discover
*/
void DataCommunicator::setDiscoverTag(int tag)
{
    m_customDiscoverTag = (tag != TAG_AUTO);
    if (m_customDiscoverTag) {
        m_discoverTag = tag;
    } else {
        if (m_rank == 0) {
            m_discoverTag = communications::tags().generate();
        }

        MPI_Bcast(&m_discoverTag, 1, MPI_INT, 0, m_communicator);
    }
}

/*!
    Sets the tag to be used for notifications

    By default, a unique tag is generated in the constructor. However, using
    this function, it is possible to assign a custom tag.

    \param tag is the custom tag to be used for notifications
*/
void DataCommunicator::setNotificationTag(int tag)
{
    m_customNotificationTag = (tag != TAG_AUTO);
    if (m_customNotificationTag) {
        m_notificationTag = tag;
    } else {
        if (m_rank == 0) {
            m_notificationTag = communications::tags().generate();
        }

        MPI_Bcast(&m_notificationTag, 1, MPI_INT, 0, m_communicator);
    }
}

/*!
    Gets the tag to be used for data exchange communications

    \result The tag to be used for data exchange communications.
*/
int DataCommunicator::getTag() const
{
    return getExchangeTag();
}

/*!
    Gets the tag to be used for data exchange communications

    \result The tag to be used for data exchange communications.
*/
int DataCommunicator::getExchangeTag() const
{
    return m_exchangeTag;
}

/*!
    Gets the tag to be used for data size discover communications

    \result The tag to be used for data size discover communications.
*/
int DataCommunicator::getDiscoverTag() const
{
    return m_discoverTag;
}

/*!
    Gets the tag to be used for notifications

    \result The tag to be used for notifications.
*/
int DataCommunicator::getNotificationTag() const
{
    return m_notificationTag;
}

/*!
    Set the recevies in "continuous" mode.

    When the recives are in "continuous" mode they will be restarted as soon
    as they end. In this way there is always a buffer ready for receiveing
    data.

    Calling this function will cancels all current receives.

    \param enabled if set to true enables the "continuous" mode.
*/
void DataCommunicator::setRecvsContinuous(bool enabled)
{
    if (m_recvsContinuous == enabled) {
        return;
    }

    cancelAllRecvs(true);

    int nRecvBuffers = m_recvBuffers.size();
    for (int k = 0; k < nRecvBuffers; ++k) {
        size_t size = m_recvBuffers[k].getSize();
        m_recvBuffers[k] = RecvBuffer(size, enabled);
    }

    m_recvsContinuous = enabled;

    startAllRecvs();
}

/*!
    Checks if the recevies are in "continuous" mode.

    \result Returns true if the receives are in "continuous" mode,
    false otherwise.
*/
bool DataCommunicator::areRecvsContinuous()
{
    return m_recvsContinuous;
}

/*!
    Discover the sends inspecting the receives that the user has already set.
*/
void DataCommunicator::discoverSends()
{
    // Cancel current sends
    clearAllSends();

    // Send the data sizes with a synchronous send
    int nRecvs = getRecvCount();

    std::vector<long> discoverSizes(nRecvs, 0);
    std::vector<MPI_Request> discoverRequests(nRecvs, MPI_REQUEST_NULL);
    for (int i = 0; i < nRecvs; ++i) {
        int rank = m_recvRanks[i];
        RecvBuffer &buffer = m_recvBuffers[i];
        discoverSizes[i] = buffer.getSize();
        MPI_Issend(discoverSizes.data() + i, 1, MPI_LONG, rank, m_discoverTag, m_communicator, discoverRequests.data() + i);
    }

    // Receive the data sizes and set the sends
    MPI_Request exchangeCompletedRequest = MPI_REQUEST_NULL;
    while (true) {
        // If there are messagea available receive them and set the sends
        int messageAvailable = 1;
        while (messageAvailable) {
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, m_discoverTag, m_communicator, &messageAvailable, &status);
            if (messageAvailable) {
                long dataSize;
                MPI_Recv(&dataSize, 1, MPI_LONG, status.MPI_SOURCE, m_discoverTag, m_communicator, MPI_STATUS_IGNORE);
                setSend(status.MPI_SOURCE, dataSize);
            }
        }

        // If all the sends are complete notify it
        if (exchangeCompletedRequest == MPI_REQUEST_NULL) {
            int discoverSendsCompleted;
            MPI_Testall(discoverRequests.size(), discoverRequests.data(), &discoverSendsCompleted, MPI_STATUSES_IGNORE);
            if (discoverSendsCompleted) {
                MPI_Ibarrier(m_communicator, &exchangeCompletedRequest);
            }
        }

        // If all sends are completed, check if also the other processes have
        // completed the sends. Sice these are synchronous sends, they will
        // be makred as completed only when the corresponding receive has
        // completed. When all processes have completed the send/recevies
        // all sizes have been exchanged.
        if (exchangeCompletedRequest != MPI_REQUEST_NULL) {
            int exchangeCompleted = 0;
            MPI_Test(&exchangeCompletedRequest, &exchangeCompleted, MPI_STATUS_IGNORE);
            if (exchangeCompleted) {
                break;
            }
        }
    }

    // At this point we are sure that all the processes have sent the data
    // sizes and the data sizes have reached the destination, however it is
    // not guaranteed that the data sizes have already been received (the
    // message is on the destination process, but we don't know if it has
    // been processed and received). Consequently, we can't be sure that all
    // sends have been set.
    //
    // For the above reason it is not possible to call this function multiple
    // times without receiving the messages or changes the exchanges between
    // one call an another. Nor it is possible to call a discover recives
    // followed by a discover sends (or viceversa) without receviing the
    // messages or canceling them in between the two calls.
}

/*!
    Discover the receives inspecting the sends that the user has already set.
*/
void DataCommunicator::discoverRecvs()
{
    // Cancel current receives
    clearAllRecvs();

    // Send the data sizes with a synchronous send
    int nSends = getSendCount();

    std::vector<long> discoverSizes(nSends, 0);
    std::vector<MPI_Request> discoverRequests(nSends, MPI_REQUEST_NULL);
    for (int i = 0; i < nSends; ++i) {
        int rank = m_sendRanks[i];
        SendBuffer &buffer = m_sendBuffers[i];
        discoverSizes[i] = buffer.getSize();
        MPI_Issend(discoverSizes.data() + i, 1, MPI_LONG, rank, m_discoverTag, m_communicator, discoverRequests.data() + i);
    }

    // Receive the data sizes and set the receives
    MPI_Request exchangeCompletedRequest = MPI_REQUEST_NULL;
    while (true) {
        // If there are messagea available receive them and set the receives
        int messageAvailable = 1;
        while (messageAvailable) {
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, m_discoverTag, m_communicator, &messageAvailable, &status);
            if (messageAvailable) {
                long dataSize;
                MPI_Recv(&dataSize, 1, MPI_LONG, status.MPI_SOURCE, m_discoverTag, m_communicator, MPI_STATUS_IGNORE);
                setRecv(status.MPI_SOURCE, dataSize);
            }
        }

        // If all the sends are complete notify it
        if (exchangeCompletedRequest == MPI_REQUEST_NULL) {
            int discoverSendsCompleted;
            MPI_Testall(discoverRequests.size(), discoverRequests.data(), &discoverSendsCompleted, MPI_STATUSES_IGNORE);
            if (discoverSendsCompleted) {
                MPI_Ibarrier(m_communicator, &exchangeCompletedRequest);
            }
        }

        // If all sends are completed, check if also the other processes have
        // completed the sends. Sice these are synchronous sends, they will
        // be makred as completed only when the corresponding receive has
        // completed. When all processes have completed the send/recevies
        // all sizes have been exchanged.
        if (exchangeCompletedRequest != MPI_REQUEST_NULL) {
            int exchangeCompleted = 0;
            MPI_Test(&exchangeCompletedRequest, &exchangeCompleted, MPI_STATUS_IGNORE);
            if (exchangeCompleted) {
                break;
            }
        }
    }

    // At this point we are sure that all the processes have sent the data
    // sizes and the data sizes have reached the destination, however it is
    // not guaranteed that the data sizes have already been received (the
    // message is on the destination process, but we don't know if it has
    // been processed and received). Consequently, we can't be sure that all
    // receives have been set.
    //
    // For the above reason it is not possible to call this function multiple
    // times without receiving the messages or changes the exchanges between
    // one call an another. Nor it is possible to call a discover recives
    // followed by a discover sends (or viceversa) without receviing the
    // messages or canceling them in between the two calls.
}

/*!
    Clear the send associated to the specified processor.

    \param rank is the ranks associated to the send that will be cleared
*/
void DataCommunicator::clearSend(int rank)
{
    if (m_sendIds.count(rank) == 0) {
        return;
    }

    // Cancel waiting send
    cancelSend(rank);

    // Remove the send associated to the specified rank
    int id = m_sendIds[rank];
    m_sendIds.erase(rank);
    for (auto &entry : m_sendIds) {
        if (entry.second > id) {
            entry.second--;
        }
    }

    m_sendRanks.erase(m_sendRanks.begin() + id);
    m_sendRequests.erase(m_sendRequests.begin() + id);
    m_sendBuffers.erase(m_sendBuffers.begin() + id);
}

/*!
    Clear the receive associated to the specified processor.

    \param rank is the ranks associated to the receive that will be cleared
*/
void DataCommunicator::clearRecv(int rank)
{
    if (m_recvIds.count(rank) == 0) {
        return;
    }

    // Cancel waiting recv
    cancelRecv(rank);

    // Remove the recv associated to the specified rank
    int id = m_recvIds[rank];
    m_recvIds.erase(rank);
    for (auto &entry : m_recvIds) {
        if (entry.second > id) {
            entry.second--;
        }
    }

    m_recvRanks.erase(m_recvRanks.begin() + id);
    m_recvRequests.erase(m_recvRequests.begin() + id);
    m_recvBuffers.erase(m_recvBuffers.begin() + id);
}

/*!
    Clear all the sends.

    \param synchronous if this parameter is set to true, the function will
    cancel the send requests of the current rank and than will wait until
    all remote processes that sends data to the current rank cancel their
    send requests as well. This guarantees that, at the end of the function,
    both the current rank and the remote processes involved in a data exchange
    with the current rank have canceled their send requests.
*/
void DataCommunicator::clearAllSends(bool synchronous)
{
    // Cancel waiting sends
    cancelAllSends(synchronous);

    // Clear the sends
    m_sendRanks.clear();
    m_sendIds.clear();
    m_sendRequests.clear();
    m_sendBuffers.clear();
}

/*!
    Clear all the receives.

    \param synchronous if this parameter is set to true, the function will
    cancel the receive requests of the current rank and than will wait until
    all remote processes that receives data from the current rank cancel their
    receive requests as well. This guarantees that, at the end of the function,
    both the current rank and all remote processes involved in a data exchange
    with it have canceled their receive requests.
*/
void DataCommunicator::clearAllRecvs(bool synchronous)
{
    // Cancel waiting receives
    cancelAllRecvs(synchronous);

    // Clear data associated to the recevies
    m_recvRanks.clear();
    m_recvIds.clear();
    m_recvRequests.clear();
    m_recvBuffers.clear();
}

/*!
    Set send information for the specified rank

    \param rank is the rank of the processor associated to the send
    \param length is the length, expressed in bytes, of the data to be sent
*/
void DataCommunicator::setSend(int rank, long length)
{
    // Clear the send associated to the rank
    clearSend(rank);

    // Set send info
    int id = m_sendIds.size();
    m_sendIds[rank] = id;

    m_sendRanks.push_back(rank);
    m_sendRequests.push_back(MPI_REQUEST_NULL);
    m_sendBuffers.emplace_back(length);
}

/*!
    Set recevie information for the specified rank

    \param rank is the rank of the processor associated to the receive
    \param length is the length, expressed in bytes, of the data to be received
*/
void DataCommunicator::setRecv(int rank, long length)
{
    // Clear the recv associated to the rank
    clearRecv(rank);

    // Set recv info
    int id = m_recvIds.size();
    m_recvIds[rank] = id;

    m_recvRanks.push_back(rank);
    m_recvRequests.push_back(MPI_REQUEST_NULL);
    m_recvBuffers.emplace_back(length, m_recvsContinuous);

    // If the receives are continous start the receive
    if (areRecvsContinuous()) {
        _startRecv(rank);
    }
}

/*!
    Resize the send associated to the specified rank

    \param rank is the rank of the processor associated to the send
    \param size is the size, expressed in bytes, of the send
*/
void DataCommunicator::resizeSend(int rank, long size)
{
    // If there is no send associate to the specified rank we have to set
    // a new send from scratch
    if (m_sendIds.count(rank) == 0) {
        setSend(rank, size);
        return;
    }

    // Cancel the send associated to the processor
    cancelSend(rank);

    // Resize the buffer
    int id = m_sendIds[rank];
    m_sendBuffers[id].setSize(size);
}

/*!
    Resize the receive associated to the specified rank

    \param rank is the rank of the processor associated to the receive
    \param size is the size, expressed in bytes, of the receive
*/
void DataCommunicator::resizeRecv(int rank, long size)
{
    // If there is no receive associate to the specified rank we have to set
    // a new receive from scratch
    if (m_recvIds.count(rank) == 0) {
        setRecv(rank, size);
        return;
    }

    // Cancel the receive associated to the processor
    cancelRecv(rank);

    // Resize the buffer
    int id = m_recvIds[rank];
    m_recvBuffers[id].setSize(size);
}

/*!
    Counts the number of sends that will be performed.

    \result The number of sends that will be performed.
*/
int DataCommunicator::getSendCount()
{
    return m_sendBuffers.size();
}

/*!
    Counts the number of receives that will be performed.

    \result The number of receives that will be performed.
*/
int DataCommunicator::getRecvCount()
{
    return m_recvBuffers.size();
}

/*!
    Get a constant reference to the list of ranks for which a send has been set.

    \result A constant reference to the list of ranks for with a send has been
    set.
*/
const std::vector<int> & DataCommunicator::getSendRanks() const
{
    return m_sendRanks;
}

/*!
    Get a constant reference to the list of ranks for which a receive has been
    set.

    \result A constant reference to the list of ranks for with a receive has
    been set.
*/
const std::vector<int> & DataCommunicator::getRecvRanks() const
{
    return m_recvRanks;
}

/*!
    Gets the send buffer associated with the requested rank

    \param rank is the rank for which the buffer is requested
    \result The send buffer associated with the requested rank.
*/
SendBuffer & DataCommunicator::getSendBuffer(int rank)
{
    int id = m_sendIds.at(rank);

    return m_sendBuffers[id];
}

/*!
    Gets the receive buffer associated with the specified rank

    \param rank is the rank for which the buffer is requested
    \result The receive buffer associated with the specified rank.
*/
RecvBuffer & DataCommunicator::getRecvBuffer(int rank)
{
    int id = m_recvIds.at(rank);

    return m_recvBuffers[id];
}

/*!
    Starts sending the data to the specified rank

    \param dstRank is the destination rank
*/
void DataCommunicator::startSend(int dstRank)
{
    // Wait for the previous send to finish
    waitSend(dstRank);

    // If the buffer is a double buffer, swap it
    int id = m_sendIds.at(dstRank);
    SendBuffer &sendBuffer = m_sendBuffers[id];
    if (sendBuffer.isDouble()) {
        sendBuffer.swap();
    }

    // Start the send
    _startSend(dstRank);
}

/*!
    Starts sending the data to all the ranks
*/
void DataCommunicator::startAllSends()
{
    for (int rank : m_sendRanks) {
        startSend(rank);
    }
}

/*!
    Internal function that starts sending the data to the specified rank

    \param dstRank is the destination rank
*/
void DataCommunicator::_startSend(int dstRank)
{
    // Get the buffer
    int id = m_sendIds.at(dstRank);
    SendBuffer &sendBuffer = m_sendBuffers[id];
    OBinaryStream &buffer = sendBuffer.getBack();

    // Start the send
    int chunkSize = buffer.getChunkSize();
    MPI_Datatype chunkDataType = getChunkDataType(chunkSize);

    MPI_Isend(buffer.data(), buffer.getChunkCount(), chunkDataType, dstRank, m_exchangeTag,
              m_communicator, &m_sendRequests[id]);
}

/*!
    Starts receiving the data from the specified rank

    \param srcRank is the source rank
*/
void DataCommunicator::startRecv(int srcRank)
{
    // Wait for the previous receive to finish
    waitRecv(srcRank);

    // Start the recevier
    _startRecv(srcRank);
}

/*!
    Starts receiving the data from all the ranks
*/
void DataCommunicator::startAllRecvs()
{
    for (int rank : m_recvRanks) {
        startRecv(rank);
    }
}

/*!
    Internal function taht starts receiving the data from the specified rank

    \param srcRank is the source rank
*/
void DataCommunicator::_startRecv(int srcRank)
{
    // Reset the position of the buffer
    int id = m_recvIds.at(srcRank);
    IBinaryStream &buffer = m_recvBuffers[id].getBack();
    buffer.seekg(0);

    // Start the receive
    int chunkSize = buffer.getChunkSize();
    MPI_Datatype chunkDataType = getChunkDataType(chunkSize);

    MPI_Irecv(buffer.data(), buffer.getChunkCount(), chunkDataType, srcRank, m_exchangeTag,
              m_communicator, &m_recvRequests[id]);
}

/*!
    Waits for any send to completes and returns the associated rank.

    If there are no active sends, the call returns MPI_UNDEFINED.

    \param blackList is a list of ranks whose sends don't have to be waited
    for
    \result The rank of the completed send or MPI_UNDEFINED if there was
    no active sends.
*/
int DataCommunicator::waitAnySend(const std::vector<int> &blackList)
{
    // Exclude blackListed ranks
    std::vector<MPI_Request> requestList(m_sendRequests);
    for (const int rank : blackList) {
        int id = m_sendIds.at(rank);
        requestList[id] = MPI_REQUEST_NULL;
    }

    // Wait for a send to complete
    int id;
    MPI_Waitany(requestList.size(), requestList.data(), &id, MPI_STATUS_IGNORE);
    if (id == MPI_UNDEFINED) {
        return MPI_UNDEFINED;
    }

    m_sendRequests[id] = requestList[id];

    // Reset the position of the buffer
    m_sendBuffers[id].seekg(0);

    // Return the rank associated to the completed send
    return m_sendRanks[id];
}

/*!
    Waits for the send associate to the sepcified rank to complete.

    \param rank is the rank associated to the send to wait for
*/
void DataCommunicator::waitSend(int rank)
{
    // Wait for the send to complete
    int id = m_sendIds.at(rank);
    auto request = m_sendRequests[id];
    if (request == MPI_REQUEST_NULL) {
        return;
    }

    MPI_Wait(&m_sendRequests[id], MPI_STATUS_IGNORE);

    // Reset the position of the buffer
    m_sendBuffers[id].seekg(0);
}

/*!
    Waits for all the sends to complete.
*/
void DataCommunicator::waitAllSends()
{
    if (m_sendRequests.size() == 0) {
        return;
    }

    // Wait for all sends to complete
    MPI_Waitall(m_sendRequests.size(), m_sendRequests.data(), MPI_STATUS_IGNORE);

    // Reset the position of the buffers
    for (auto &buffer : m_sendBuffers) {
        buffer.seekg(0);
    }
}

/*!
    Waits for any receive to completes and returns the associated rank.

    If there are no active recevies, the call returns MPI_UNDEFINED.

    \param blackList is a list of ranks whose recevies don't have to be waited
    for
    \result The rank of the completed receive or MPI_UNDEFINED if there was
    no active receives.
*/
int DataCommunicator::waitAnyRecv(const std::vector<int> &blackList)
{
    // Exclude blackListed ranks
    std::vector<MPI_Request> requestList(m_recvRequests);
    for (const int rank : blackList) {
        int id = m_recvIds.at(rank);
        requestList[id] = MPI_REQUEST_NULL;
    }

    // Wait for a receive to complete
    int id;
    MPI_Waitany(requestList.size(), requestList.data(), &id, MPI_STATUS_IGNORE);
    if (id == MPI_UNDEFINED) {
        return MPI_UNDEFINED;
    }

    m_recvRequests[id] = requestList[id];

    // If the buffer is a double buffer, swap it
    RecvBuffer &recvBuffer = m_recvBuffers[id];
    if (recvBuffer.isDouble()) {
        recvBuffer.swap();
    }

    // Rank of the request
    int rank = m_recvRanks[id];

    // Restart the recevie
    if (areRecvsContinuous()) {
        _startRecv(rank);
    }

    // Return the rank associated to the completed receive
    return rank;
}

/*!
    Waits for the receive associate to the sepcified rank to complete.

    \param rank is the rank associated to the receive to wait for
*/
void DataCommunicator::waitRecv(int rank)
{
    // Wait for the receive to complete
    int id = m_recvIds.at(rank);
    auto request = m_recvRequests[id];
    if (request == MPI_REQUEST_NULL) {
        return;
    }

    MPI_Wait(&m_recvRequests[id], MPI_STATUS_IGNORE);

    // If the buffer is a double buffer, swap it
    RecvBuffer &recvBuffer = m_recvBuffers[id];
    if (recvBuffer.isDouble()) {
        recvBuffer.swap();
    }

    // Restart the recevie
    if (areRecvsContinuous()) {
        _startRecv(rank);
    }
}

/*!
    Waits for all the receives to complete.
*/
void DataCommunicator::waitAllRecvs()
{
    if (m_recvRequests.size() == 0) {
        return;
    }

    // Wait for all the receives to complete
    MPI_Waitall(m_recvRequests.size(), m_recvRequests.data(), MPI_STATUS_IGNORE);

    // Swap double buffers
    for (RecvBuffer &buffer : m_recvBuffers) {
        if (buffer.isDouble()) {
            buffer.swap();
        }
    }

    // Restart all the receives
    if (areRecvsContinuous()) {
        for (int rank : m_recvRanks) {
            _startRecv(rank);
        }
    }
}

/*!
    Checkis if the send associated to the sepcified rank is active.

    A send is considered active if the associated request is not set to
    'MPI_REQUEST_NULL'.

    \param rank is the rank associated to the send
*/
bool DataCommunicator::isSendActive(int rank)
{
    int id = m_sendIds[rank];

    return (m_sendRequests[id] != MPI_REQUEST_NULL);
}

/*!
    Checkis if the receive associated to the sepcified rank is active.

    A receive is considered active if the associated request is not set to
    'MPI_REQUEST_NULL'.

    \param rank is the rank associated to the receive
*/
bool DataCommunicator::isRecvActive(int rank)
{
    int id = m_recvIds[rank];

    return (m_recvRequests[id] != MPI_REQUEST_NULL);
}

/*!
    Cancels the send associated to the specified rank.

    \param rank is the rank associated to the send to cancel
*/
void DataCommunicator::cancelSend(int rank)
{
    if (m_sendIds.count(rank) == 0) {
        return;
    }

    int id = m_sendIds[rank];
    if (m_sendRequests[id] == MPI_REQUEST_NULL) {
        return;
    }

    MPI_Cancel(&m_sendRequests[id]);
    MPI_Request_free(&m_sendRequests[id]);
}

/*!
    Cancels the receive associated to the specified rank.

    \param rank is the rank associated to the receive to cancel
*/
void DataCommunicator::cancelRecv(int rank)
{
    if (m_recvIds.count(rank) == 0) {
        return;
    }

    int id = m_recvIds[rank];
    if (m_recvRequests[id] == MPI_REQUEST_NULL) {
        return;
    }

    MPI_Cancel(&m_recvRequests[id]);
    MPI_Request_free(&m_recvRequests[id]);
}

/*!
    Cancels all the sends.

    \param synchronous if this parameter is set to true, the function will
    cancel the send requests of the current rank and than will wait until
    all remote processes that sends data to the current rank cancel their
    send requests. This guarantees that, at the end of the function, both
    the current rank and all remote processes involved in a data exchange
    with it have canceled their send requests. When synchronous mode is
    enabled send and receive requests have to match among all the processes.
*/
void DataCommunicator::cancelAllSends(bool synchronous)
{
    // If no synchronization is needed just call the standard overload
    if (!synchronous) {
        cancelAllSends();
        return;
    }

    // Notify that the sends have been canceled
    //
    // The tag of the notifications has to be different from the tag used for
    // data size discovery. We may have a scenario like the following:
    //
    //    discoverRecvs(); # <-- This call is using the discover tag
    //
    //    startAllRecvs();
    //
    //    <User fill send buffers>
    //
    //    startAllSends();
    //
    //    waitAllRecvs();
    //
    //    waitAllSends();
    //
    //    cancelAllSends(true); # <-- This call is using the notification tag
    //
    //    cancelAllrecvs(true); # <-- This call is using the notification tag
    //
    // If a process completes all its exchanges, it may reach the function
    // that cancels the sends while other processes are still discovering
    // the receives. Thats's because data size discovery is a non-blocking
    // operation. At the end of data size discovery we know that the data
    // sizes for the recevies has reached the destination processes, but it
    // not guaranteed that the messages have been processed. What can happend
    // is that a process can exit from the function while other processes are
    // still receving/processing the recevied messages (so they are still
    // setting the receives). The notification issued when canceling the
    // sends may than be catched by the processes that are still discovering
    // the receives. To avoid this discover tags and notification tags have
    // to be different.

    // Start auxiliary receives for handling synchronization
    //
    // We will receive a notification from all the processes for which a
    // receive has been set. The notification will contain the status of
    // the send request on that processor that matches the receive of this
    // rank. A status equal to zero means that the send has been successfully
    // canceled.
    int nRecvs = getRecvCount();
    std::vector<MPI_Request> recvNotificationRequests(nRecvs);

    std::vector<int> remoteSendStatus(nRecvs, 1);
    for (int i = 0; i < nRecvs; ++i) {
        int rank = m_recvRanks[i];
        MPI_Irecv(remoteSendStatus.data() + i, 1, MPI_INT, rank, m_notificationTag, m_communicator, &recvNotificationRequests[i]);
    }

    // Cancel all the sends
    //
    // We will send a notification to all processes for wich a send has been
    // set. The notification will contain the status of the send request. A
    // status equal to zero means that the send request has been successfully
    // canceled.
    int nSends = getSendCount();
    std::vector<MPI_Request> sendNotificationRequests(nSends);
    for (int i = 0; i < nSends; ++i) {
        int rank = m_sendRanks[i];

        // Cancel the send
        cancelSend(rank);

        // Notify the other process that the send has been canceled
        int sendStatus = 0;
        MPI_Isend(&sendStatus, 1, MPI_INT, rank, m_notificationTag, m_communicator, &sendNotificationRequests[i]);
    }

    // Synchronize with other processes
    //
    // Wait until all processes involved in data exchanging cancel their
    // send requests.
    std::size_t nMissingNotifications = recvNotificationRequests.size();
    while (nMissingNotifications > 0) {
        int i;
        MPI_Status status;
        MPI_Waitany(recvNotificationRequests.size(), recvNotificationRequests.data(), &i, &status);
        if (remoteSendStatus[i] == 1) {
            log::cout() << "Unable to properly cancel the sends.";
            MPI_Abort(m_communicator, 1);
        }

        --nMissingNotifications;
    }

    // Wait until the notifications have been sent
    MPI_Waitall(sendNotificationRequests.size(), sendNotificationRequests.data(), MPI_STATUS_IGNORE);

#if BITPIT_ENABLE_DEBUG
    // Wait until all processes cancel the sends
    MPI_Barrier(m_communicator);
#endif
}

/*!
    Cancels all the receives.

    \param synchronous if this parameter is set to true, the function will
    cancel the receive requests of the current rank and than will wait until
    all remote processes that receives data from the current rank cancel their
    receive requests. This guarantees that, at the end of the function, both
    the current rank and all remote processes involved in a data exchange with
    it have canceled their receive requests. When synchronous mode is enabled
    send and receive requests have to match among all the processes.
*/
void DataCommunicator::cancelAllRecvs(bool synchronous)
{
    // If no synchronization is needed just call the standard overload
    if (!synchronous) {
        cancelAllRecvs();
        return;
    }

    // Notify that the sends have been canceled
    //
    // The tag of the notifications has to be different from the tag used for
    // data size discovery. See comment in 'cancelAllSends' function for an
    // explanation.

    // Start auxiliary receives for handling synchronization
    //
    // We will receive a notification from all the processes for which a
    // send has been set. The notification will contain the status of the
    // receive request on that processor that matches the send of this rank.
    // A status equal to zero means that the receive has been successfully
    // canceled.
    int nSends = getSendCount();
    std::vector<MPI_Request> recvNotificationRequests(nSends);

    std::vector<int> remoteReceiveStatus(nSends, 1);
    for (int i = 0; i < nSends; ++i) {
        int rank = m_sendRanks[i];
        MPI_Irecv(remoteReceiveStatus.data() + i, 1, MPI_INT, rank, m_notificationTag, m_communicator, &recvNotificationRequests[i]);
    }

    // Cancel all the receives
    //
    // We will send a notification to all processes for wich a receive has been
    // set. The notification will contain the status of the receive request. A
    // status equal to zero means that the receive request has been successfully
    // canceled.
    int nRecvs = getRecvCount();
    std::vector<MPI_Request> sendNotificationRequests(nRecvs);
    for (int i = 0; i < nRecvs; ++i) {
        int rank = m_recvRanks[i];

        // Cancel the send
        cancelRecv(rank);

        // Notify other process that the received has been canceled
        int receiveStatus = 0;
        MPI_Isend(&receiveStatus, 1, MPI_INT, rank, m_notificationTag, m_communicator, &sendNotificationRequests[i]);
    }

    // Synchronize with other processes
    //
    // Wait until all processes involved in data exchanging cancel their
    // receive requests.
    std::size_t nMissingNotifications = recvNotificationRequests.size();
    while (nMissingNotifications > 0) {
        int i;
        MPI_Waitany(recvNotificationRequests.size(), recvNotificationRequests.data(), &i, MPI_STATUS_IGNORE);
        if (remoteReceiveStatus[i] == 1) {
            log::cout() << "Unable to properly cancel the receives.";
            MPI_Abort(m_communicator, 1);
        }

        --nMissingNotifications;
    }

    // Wait until the notifications have been sent
    MPI_Waitall(sendNotificationRequests.size(), sendNotificationRequests.data(), MPI_STATUS_IGNORE);

#if BITPIT_ENABLE_DEBUG
    // Wait until all processes cancel the receives
    MPI_Barrier(m_communicator);
#endif
}

/*!
    Cancels all the sends.
*/
void DataCommunicator::cancelAllSends()
{
    for (int rank : m_sendRanks) {
        cancelSend(rank);
    }
}

/*!
    Cancels all the receives.
*/
void DataCommunicator::cancelAllRecvs()
{
    for (int rank : m_recvRanks) {
        cancelRecv(rank);
    }
}

/*!
    Get the MPI data type associate to a data chunk.

    \param chunkSize is the size of the data chunk
    \result The MPI data type associate to a data chunk.
*/
MPI_Datatype DataCommunicator::getChunkDataType(int chunkSize) const
{
    MPI_Datatype chunkDataType;
    if (chunkSize == 1) {
        chunkDataType = MPI_CHAR;
    } else {
        MPI_Type_contiguous(chunkSize, MPI_CHAR, &chunkDataType);
        MPI_Type_commit(&chunkDataType);
    }

    return chunkDataType;
}

}

#endif
