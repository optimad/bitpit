// ========================================================================== //
//                         - G.L.O.R.I.A. -                                   //
//                                                                            //
// Routins for exchanging data among the processors                           //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Andrea Iob
// Version  : v1.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

#if BITPIT_ENABLE_MPI==1

#include "communications.hpp"

using namespace bitpit;

/*!
	\class DataCommunicator

	\brief The DataCommunicator class provides the infrastructure needed to
	exchange data among processors.
*/

/*!
	Creates a new communicator for data exchange.
 */
DataCommunicator::DataCommunicator(MPI_Comm communicator)
	: m_communicator(communicator), m_rank(-1),
	  m_tag(MPI_ANY_TAG), m_recvsContinuous(false)
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
	Set the recevies in "continuous" mode.

	When the recives are in "continuous" mode they will be restarted as soon
	as they end. In this way there is always a buffer ready for receiveing
	data.

	\param enabled if set to true enables the "continuous" mode.
 */
void DataCommunicator::setRecvsContinuous(bool enabled)
{
	m_recvsContinuous = enabled;
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
	// Cancel waiting sends
	cancelAllRecvs();

	// Clear the recevies
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
    m_recvBuffers.emplace_back(length);
}

/*!
	Grow the send associated to the specified rank

	\param rank is the rank of the processor associated to the send
	\param growth is the growth, expressed in bytes, of the send
 */
void DataCommunicator::growSend(int rank, long growth)
{
	// If there is no send associate to the specified rank we have to set
	// a new send from scratch
	if (m_sendIds.count(rank) == 0) {
		setSend(rank, growth);
		return;
	}

	// Cancel the send associated to the processor
	cancelSend(rank);

	// Grow the buffer
	int id = m_sendIds[rank];
	m_sendBuffers[id].setCapacity(m_sendBuffers[id].capacity() + growth);
}

/*!
	Grow the receive associated to the specified rank

	\param rank is the rank of the processor associated to the receive
	\param growth is the growth, expressed in bytes, of the receive
 */
void DataCommunicator::growRecv(int rank, long growth)
{
	// If there is no receive associate to the specified rank we have to set
	// a new receive from scratch
	if (m_recvIds.count(rank) == 0) {
		setRecv(rank, growth);
		return;
	}

	// Cancel the receive associated to the processor
	cancelRecv(rank);

	// Grow the buffer
	int id = m_recvIds[rank];
	m_recvBuffers[id].setCapacity(m_recvBuffers[id].capacity() + growth);
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
	Gets the send buffer associated with the requested rank

	\param rank is the rank for which the buffer is requested
	\result The send buffer associated with the requested rank.
 */
bitpit::OBinaryStream & DataCommunicator::getSendBuffer(int rank)
{
	int id = m_sendIds.at(rank);

	return m_sendBuffers[id];
}

/*!
	Gets the receive buffer associated with the specified rank

	\param rank is the rank for which the buffer is requested
	\result The receive buffer associated with the specified rank.
 */
bitpit::IBinaryStream & DataCommunicator::getRecvBuffer(int rank)
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

	// Start the send
	int id = m_sendIds.at(dstRank);
	OBinaryStream &buffer = m_sendBuffers[id];

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
	IBinaryStream &buffer = m_recvBuffers[id];
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

	// Restart the recevie
	if (areRecvsContinuous()) {
		int rank = m_recvRanks[id];
		startRecv(rank);
	}

	// Return the rank associated to the completed receive
	return m_recvRanks[id];
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

#endif
