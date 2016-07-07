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

#include "bitpit_IO.hpp"

#include "communications.hpp"

namespace bitpit {

/*!
    \ingroup communications
    @{
*/

/*!
    \class DataCommunicator

    \brief The DataCommunicator class provides the infrastructure needed to
    exchange data among processors.
*/

int DataCommunicator::DEFAULT_TAG = 0;

/*!
    Creates a new communicator for data exchange.
*/
DataCommunicator::DataCommunicator(MPI_Comm communicator)
    : m_communicator(communicator), m_rank(-1),
    m_tag(DEFAULT_TAG), m_recvsContinuous(false)
{
    // Get MPI information
    MPI_Comm_rank(m_communicator, &m_rank);
}

/*!
    Finalizes the communicator
*/
void DataCommunicator::finalize()
{
    // Cancels all sends
    cancelAllSends();

    // Cancels all receives
    cancelAllRecvs();
}

/*!
    Sets the tag to be used for the communications

    \param tag is the tag to use
*/
void DataCommunicator::setTag(int tag)
{
    m_tag = tag;
}

/*!
    Gets the tag to be used for the communications

    \result The tag to be used for the communications.
*/
int DataCommunicator::getTag() const
{
    return m_tag;
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

    cancelAllRecvs();

    int nRecvBuffers = m_recvBuffers.size();
    for (int k = 0; k < nRecvBuffers; ++k) {
        size_t capacity = m_recvBuffers[k].capacity();
        m_recvBuffers[k] = RecvBuffer(capacity, enabled);
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

    For the communications needed to discover the sends it will be used
    se same tag set in the communicator.
*/
void DataCommunicator::discoverSends()
{
    discoverSends(m_tag);
}

/*!
    Discover the sends inspecting the receives that the user has already set.

    \param discoverTag is the tag to be used for the communications needed
    to discover the sends
*/
void DataCommunicator::discoverSends(int discoverTag)
{
    // Cancel current sends
    clearAllSends();

    // Send the size of the messages that the processors want to receive
    for (auto &entry : m_recvIds) {
        int rank = entry.first;
        RecvBuffer &buffer = getRecvBuffer(rank);
        long dataSize = buffer.capacity();

        MPI_Request dataSizeRequest;
        MPI_Isend(&dataSize, 1, MPI_LONG, rank, discoverTag, m_communicator, &dataSizeRequest);

        // MPI_Isend initiates an asynchronous (background) data transfer.
        // The actual data transfer might not happen unless one of the
        // MPI_Wait* or MPI_Test* calls has been made on the request.
        int completeFlag;
        MPI_Test(&dataSizeRequest, &completeFlag, MPI_STATUS_IGNORE);
    }

    // Raise a barrier to make sure that all the sends starts
    MPI_Barrier(m_communicator);

    // Receive the data size of the messages that the processors want to receive
    std::unordered_map<int, long> dataSizes;
    while (true) {
        // Probe for messages
        int messageAvailable;
        MPI_Status status;

        MPI_Iprobe(MPI_ANY_SOURCE, discoverTag, m_communicator, &messageAvailable, &status);
        if (!messageAvailable) {
            break;
        }

        // Receive the data size that will be received from the source
        long dataSize;
        MPI_Recv(&dataSize, 1, MPI_LONG, status.MPI_SOURCE, discoverTag, m_communicator, MPI_STATUS_IGNORE);
        dataSizes[status.MPI_SOURCE] = dataSize;
    }

    // Wait that all processors correctly receive the data to communicate.
    //
    // Without the barrier some processors may start the receives while other
    // are still waiting to receive the data sizes. Since the probe for the
    // data size will accept message for all sources, if a processor start
    // receiveing data it will intefere with the other processors still
    // receiving the data sizes.
    if (discoverTag == m_tag) {
        MPI_Barrier(m_communicator);
    }

    // Set the sends
    for (auto &entry : dataSizes) {
        setSend(entry.first, entry.second);
    }
}

/*!
    Discover the receives inspecting the sends that the user has already set.

    For the communications needed to discover the receives it will be used
    se same tag set in the communicator.
*/
void DataCommunicator::discoverRecvs()
{
	discoverRecvs(m_tag);
}

/*!
    Discover the receives inspecting the sends that the user has already set.

    \param discoverTag is the tag to be used for the communications needed
    to discover the receives
*/
void DataCommunicator::discoverRecvs(int discoverTag)
{
	// Cancel current receives
	clearAllRecvs();

	// Send the size of the messages that will be send
	for (auto &entry : m_sendIds) {
		int rank = entry.first;
		SendBuffer &buffer = getSendBuffer(rank);
		long dataSize = buffer.capacity();

		MPI_Request dataSizeRequest;
		MPI_Isend(&dataSize, 1, MPI_LONG, rank, discoverTag, m_communicator, &dataSizeRequest);

		// MPI_Isend initiates an asynchronous (background) data transfer.
		// The actual data transfer might not happen unless one of the
		// MPI_Wait* or MPI_Test* calls has been made on the request.
		int completeFlag;
		MPI_Test(&dataSizeRequest, &completeFlag, MPI_STATUS_IGNORE);
	}

	// Raise a barrier to make sure that all the sends starts
	MPI_Barrier(m_communicator);

	// Receive the data size of the sends
	std::unordered_map<int, long> dataSizes;
	while (true) {
		// Probe for messages
		int messageAvailable;
		MPI_Status status;

		MPI_Iprobe(MPI_ANY_SOURCE, discoverTag, m_communicator, &messageAvailable, &status);
		if (!messageAvailable) {
			break;
		}

		// Receive the data size that will be received from the source
		long dataSize;
		MPI_Recv(&dataSize, 1, MPI_LONG, status.MPI_SOURCE, discoverTag, m_communicator, MPI_STATUS_IGNORE);
		dataSizes[status.MPI_SOURCE] = dataSize;
	}

	// Wait that all processors correctly receive the data to communicate.
	//
	// Without the barrier some processors may start the receives while other
	// are still waiting to receive the data sizes. Since the probe for the
	// data size will accept message for all sources, if a processor start
	// receiveing data it will intefere with the other processors still
	// receiving the data sizes.
	if (discoverTag == m_tag) {
		MPI_Barrier(m_communicator);
	}

	// Set the receives
	for (auto &entry : dataSizes) {
		setRecv(entry.first, entry.second);
	}
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
*/
void DataCommunicator::clearAllSends()
{
    // Cancel waiting sends
    cancelAllSends();

    // Clear the sends
    m_sendRanks.clear();
    m_sendIds.clear();
    m_sendRequests.clear();
    m_sendBuffers.clear();
}

/*!
    Clear all the receives.
*/
void DataCommunicator::clearAllRecvs()
{
    // Cancel waiting receives
    cancelAllRecvs();

    // Clear data associated to the recevies
    m_recvRanks.clear();
    m_recvIds.clear();
    m_recvRequests.clear();
    m_recvBuffers.clear();
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
    m_recvBuffers.emplace_back(length);

    // If the receives are continous start the receive
    if (areRecvsContinuous()) {
        startRecv(rank);
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
    m_sendBuffers[id].setCapacity(size);
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
    m_recvBuffers[id].setCapacity(size);
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
    Get the list of ranks for with a send has been set.

    \result The list of ranks for with a send has been set.
*/
const std::vector<int> DataCommunicator::getSendRanks() const
{
    return m_sendRanks;
}

/*!
    Get the list of ranks for with a receive has been set.

    \result The list of ranks for with a receive has been set.
*/
const std::vector<int> DataCommunicator::getRecvRanks() const
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

    // Id of the buffer
    int id = m_sendIds.at(dstRank);

    // If the buffer is a double buffer, swap it
    SendBuffer &sendBuffer = m_sendBuffers[id];
    if (sendBuffer.isDouble()) {
        sendBuffer.swap();
    }

    // Start the send
    OBinaryStream &buffer = sendBuffer.getBack();

    MPI_Isend(buffer.rawData(), buffer.capacity(), MPI_CHAR, dstRank, m_tag,
            m_communicator, &m_sendRequests[id]);
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
    Starts receiving the data from the specified rank

    \param srcRank is the source rank
*/
void DataCommunicator::startRecv(int srcRank)
{
    // Wait for the previous receive to finish
    waitRecv(srcRank);

    // Reset the position of the buffer
    int id = m_recvIds.at(srcRank);
    IBinaryStream &buffer = m_recvBuffers[id].getBack();
    buffer.seekg(0);

    // Start the receive
    MPI_Irecv(buffer.rawData(), buffer.capacity(), MPI_CHAR, srcRank, m_tag,
            m_communicator, &m_recvRequests[id]);
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
    Waits for any send to completes and returns the associated rank.

    If there are no active sends, the call returns MPI_UNDEFINED.

    \result The rank of the completed send or MPI_UNDEFINED if there was
    no active sends.
*/
int DataCommunicator::waitAnySend()
{
    // Wait for a send to complete
    int id;
    MPI_Waitany(m_sendRequests.size(), m_sendRequests.data(), &id, MPI_STATUS_IGNORE);
    if (id == MPI_UNDEFINED) {
        return MPI_UNDEFINED;
    }

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

    \param if set to true
    \result The rank of the completed receive or MPI_UNDEFINED if there was
    no active receives.
*/
int DataCommunicator::waitAnyRecv()
{
    // Wait for a receive to complete
    int id;
    MPI_Waitany(m_recvRequests.size(), m_recvRequests.data(), &id, MPI_STATUS_IGNORE);
    if (id == MPI_UNDEFINED) {
        return MPI_UNDEFINED;
    }

    // If the buffer is a double buffer, swap it
    RecvBuffer &recvBuffer = m_recvBuffers[id];
    if (recvBuffer.isDouble()) {
        recvBuffer.swap();
    }

    // Rank of the request
    int rank = m_recvRanks[id];

    // Restart the recevie
    if (areRecvsContinuous()) {
        startRecv(rank);
    }

    // Return the rank associated to the completed receive
    return rank;
}

/*!
    Waits for the receive associate to the sepcified rank to complete.

    \param rank is the rank associated to the receive to wait for
    \param
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
        startRecv(rank);
    }
}

/*!
    Waits for all the receives to complete.
*/
void DataCommunicator::waitAllRecvs()
{
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
        startAllRecvs();
    }
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

}

/*!
* @}
*/

#endif
