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

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include <mpi.h>
#include <chrono>
#include<unordered_set>
#include "patch_kernel.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace chrono;

namespace bitpit {

/*!
	\ingroup patchkernel
	@{

	\class Patch
*/

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void PatchKernel::setCommunicator(MPI_Comm communicator)
{
	// Communication can be set just once
	if (isCommunicatorSet()) {
		throw std::runtime_error ("Patch communicator can be set just once");
	}

	// The communicator has to be valid
	if (communicator == MPI_COMM_NULL) {
		throw std::runtime_error ("Patch communicator is not valid");
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	MPI_Comm_dup(communicator, &m_communicator);

	// Get MPI information
	MPI_Comm_size(m_communicator, &m_nProcessors);
	MPI_Comm_rank(m_communicator, &m_rank);

	// Set parallel data for the VTK output
	setParallel(m_nProcessors, m_rank);
}

/*!
	Checks if the communicator to be used for parallel communications has
	already been set.

	\result Returns true if the communicator has been set, false otherwise.
*/
bool PatchKernel::isCommunicatorSet() const
{
	return (getCommunicator() != MPI_COMM_NULL);
}

/*!
	Gets the MPI communicator associated to the patch

	\return The MPI communicator associated to the patch.
*/
const MPI_Comm & PatchKernel::getCommunicator() const
{
	return m_communicator;
}

/*!
	Gets the MPI rank associated to the patch

	\return The MPI rank associated to the patch.
*/
int PatchKernel::getRank() const
{
	return m_rank;
}

/*!
	Gets the MPI processors in the communicator associated to the patch

	\return The MPI processors in the communicator associated to the patch
*/
int PatchKernel::getProcessorCount() const
{
	return m_nProcessors;
}

/*!
	Partitions the patch.

	\param cellRanks are the ranks of the cells after the partitioning.
*/
void PatchKernel::partition(const std::vector<int> &cellRanks)
{
	// Build the send map
	std::unordered_map<int, std::vector<long>> sendMap;

	auto cellItr = cellBegin();
	for (int k = 0; k < getInternalCount(); ++k) {
		const int &rank = cellRanks[k];
		if (rank == getRank()) {
			cellItr++;
			continue;
		}

		sendMap[rank].push_back(cellItr->getId());

		cellItr++;
	}

	// Local sender-receiver pairs
	std::vector<int> localPairs;
	localPairs.reserve(2 * sendMap.size());
	for (const auto &entry : sendMap) {
		localPairs.push_back(getRank());
		localPairs.push_back(entry.first);
	}

	// Exchange the size of the communication
	int globalPairsSize = localPairs.size();
	std::vector<int> globalPairsSizes(getProcessorCount());
	MPI_Allgather(&globalPairsSize, 1, MPI_INT, globalPairsSizes.data(), 1, MPI_INT, getCommunicator());

	std::vector<int> globalPairsOffsets(getProcessorCount());
	globalPairsOffsets[0] = 0;
	for (int i = 1; i < getProcessorCount(); ++i) {
		globalPairsOffsets[i] = globalPairsOffsets[i-1] + globalPairsSizes[i-1];
	}

	// Global sender-receiver pairs
	std::vector<int> globalPairs;
	globalPairs.resize(globalPairsOffsets.back() + globalPairsSizes.back());

	MPI_Allgatherv(localPairs.data(), localPairs.size(), MPI_INT, globalPairs.data(),
				   globalPairsSizes.data(), globalPairsOffsets.data(), MPI_INT,
                   getCommunicator());

	// Exchange the cells
	std::vector<long> emptyCellList;
	for (size_t i = 0; i < globalPairs.size(); i += 2) {
		int srcRank = globalPairs[i];
		int dstRank = globalPairs[i+1];

		std::vector<long> *cellList;
		if (srcRank == getRank()) {
			cellList = &(sendMap[dstRank]);
		} else {
			cellList = &emptyCellList;
		}

		sendCells(srcRank, dstRank, *cellList);
	}
}

/*!
    Move cells with specified IDs from process with rank snd_rank (sender) to
    process with rank rcv_rank (receiver).
    If the rank the process currently hosting the mesh is neither the sender or
    the receiver, a notification is received in case ghost cells has changed hosting.

    \param[in] snd_rank sender rank
    \param[in] rcv_rank receiver rank
    \param[in] cell_list list of cells to be moved
 */
void PatchKernel::sendCells(const unsigned short &snd_rank, const unsigned short &rcv_rank, const vector< long > &cell_list)
{

// ========================================================================== //
// SCOPE VARIABLES                                                            //
// ========================================================================== //

// Debug variables
/*DEBUG*/stringstream                           out_name;
/*DEBUG*/ofstream                               out;
/*DEBUG*/high_resolution_clock::time_point      t0, t1;
/*DEBUG*/duration<double>                       time_span;

/*DEBUG*/{
/*DEBUG*/    out_name << "DEBUG_rank_" << m_rank << ".log";
/*DEBUG*/    out.open(out_name.str(), ifstream::app);
/*DEBUG*/    out_name.str("");
/*DEBUG*/    out << endl;
/*DEBUG*/    out << "** RANK #" << m_rank << "** NEW COMMUNICATION **" << endl;
/*DEBUG*/    out << endl;
/*DEBUG*/}

// ========================================================================== //
// SENDER                                                                     //
// ========================================================================== //
if (m_rank == snd_rank)
{

    // Scope variables ====================================================== //

    // Variable required for communicators
    long                                        buff_size;

    // Variables required for cells communication
    long                                        n_cells = cell_list.size();
    unordered_map<long, long>                   cell_map;

    // Variables required for ghosts communication
    long                                        n_ghosts;
    vector< pair<long, pair<long, short> > >    ghost_list;
    unordered_map<long, long>                   ghost_map;
    unordered_map<long, long>                   sender_ghost_new_ids;

    // Variables required for vertex communication
    long                                        n_vertex;
    vector<long>                                vertex_list;
    unordered_map<long, long>                   vertex_map;
    
    // Initialize data structures =========================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), initializing data structure {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {

        // Scope variables -------------------------------------------------- //
        long                                    cell_count = 0;
        vector<long>::const_iterator            i, e = cell_list.cend();

        // Initialize lists of cells for communication ---------------------- //
        cell_count = 0;
        for ( i = cell_list.cbegin(); i != e; ++i ) {
            cell_map[ *i ] = cell_count;
            ++cell_count;
        } //next i

/*DEBUG*/{
/*DEBUG*/    out << "    sending " << n_cells;
// /*DEBUG*/    out << " cell(s): " << cell_list
/*DEBUG*/    out << endl;
// /*DEBUG*/    out << "    cell_map: ";
// /*DEBUG*/    unordered_map<long, long>::const_iterator      ii;
// /*DEBUG*/    for (ii = cell_map.begin(); ii != cell_map.end(); ++ii) {
// /*DEBUG*/        out << "(" << ii->first << "->" << ii->second << ") ";
// /*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

        // Initialize lists for ghost communication ------------------------- //
        ghost_list.reserve( cell_list.size() );

        // Initialize lists for vertex communication ------------------------ //
        vertex_list.reserve( 4*cell_list.size() );

    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Create list of ghosts ================================================ //
    // The list of ghost stores the ghost cell already existing on the current
    // rank (but not the ghost cells towards the receiver) and the candidate
    // ghosts. Candidate ghosts are cell within the cell list (i.e. those
    // cells which will be sent to the receiver), which have at least one neighbour
    // which is a internal cell.
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), creating list of ghosts {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                     n_neighs;
        int                                     k, j;
        vector<long>                            neighs;
        long                                    ghost_counter;
        vector<long>::const_iterator            i, e;
        unordered_map<short, unordered_map<long, long> >::const_iterator   m;
        unordered_map<long, long>::const_iterator                          n;
        unordered_map<long, bool>::const_iterator                          l;
        unordered_map<long, long>::const_iterator                          ii, ee;

        // Extract receiver's ghosts ---------------------------------------- //
        ghost_counter = 0;
        e = cell_list.cend();
        for (i = cell_list.cbegin(); i < e; ++i) {
            neighs = findCellNeighs(*i);
            n_neighs = neighs.size();
            for (k = 0; k < n_neighs; ++k) {
                long neigh_idx = neighs[k];

                // Discard cells that are internal for the receiver
                if (cell_map.count(neigh_idx) > 0) {
                    continue;
                }

                // Discard cells that are a ghost for the receiver
                if ( find(m_ghost2id[rcv_rank].begin(),
                          m_ghost2id[rcv_rank].end(),
                          UnaryPredicate<const long, long>(neigh_idx) ) != m_ghost2id[rcv_rank].end() ) {
                    continue;
                }

                // Update ghost list for interior cells
                if ( m_cells[neigh_idx].isInterior() ) {
                    ghost_list.push_back( pair<long, pair<long, short> >(neigh_idx, pair<long, short>(neigh_idx, snd_rank) ) );
                    ghost_map[neigh_idx] = ghost_counter + n_cells;
                    ++ghost_counter;
                }

                // Update ghost list for ghost cells (w.r.t to another process)
                else {
                    for (m = m_ghost2id.begin(); m != m_ghost2id.end(); ++m) {
                        n = find( m->second.begin(), m->second.end(), UnaryPredicate<const long, long>(neigh_idx) );
                        if ( n != m->second.end() ) {
                            ghost_list.push_back( pair<long, pair<long, short> >(neigh_idx, pair<long, short>(n->first, m->first) ) );
                            ghost_map[neigh_idx] = ghost_counter + n_cells;
                            ++ghost_counter;
                        }
                   } //next m
                }
	    }
	}
        n_ghosts = ghost_counter;
/*DEBUG*/{
/*DEBUG*/    out << "    sending " << n_ghosts;
// /*DEBUG*/    out << " ghost(s): ";
// /*DEBUG*/    vector< pair< long, pair<long, short> > >::iterator    kk;
// /*DEBUG*/    for (kk = ghost_list.begin(); kk != ghost_list.end(); ++kk) {
// /*DEBUG*/        out << "[" << kk->first << ", (" << kk->second.first << ", " << kk->second.second << ")], ";
// /*DEBUG*/    } //next kk
/*DEBUG*/    out << endl;
// /*DEBUG*/    out << "    ghost_map: ";
// /*DEBUG*/    unordered_map<long, long>::const_iterator ii;
// /*DEBUG*/    for (ii = ghost_map.begin(); ii != ghost_map.end(); ++ii) {
// /*DEBUG*/           out << "(" << ii->first << "->" << ii->second << ") ";
// /*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

        // Build list of sender's ghost with repsect to the recevider ------- //
        //
        // If the sender is sending all the cells, there's no point in finding
        // sender's ghost: after the communication the sender will contain no
        // cells, neither internal nor ghosts.
        if (cell_list.size() < (unsigned long) m_nInternals) {
            ee = ghost_map.cend();
            for (ii = ghost_map.cbegin(); ii != ee; ++ii) {
                const long recv_ghost_idx = ii->first;

                // Sender ghosts can be identified as the sent cells that are
                // neighbours of the receiver's ghosts and are not receiver's
                // ghost themselves.
                neighs = findCellNeighs(recv_ghost_idx);
                n_neighs = neighs.size();
                for (j = 0; j < n_neighs; ++j) {
                    long send_guess_ghost = neighs[j];
                    if ( sender_ghost_new_ids.count(send_guess_ghost) > 0 ) {
                        continue;
                    } else if ( ghost_map.count(send_guess_ghost) > 0 ) {
                        continue;
                    } else if ( cell_map.count(send_guess_ghost) == 0 ) {
                        continue;
                    }

                    sender_ghost_new_ids.insert({{send_guess_ghost, Element::NULL_ID}});
                }
            }
        }

    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Create list of vertices ============================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), creating list of vertices {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                     j, n_vertices;
        long                                    vertex_counter;
        long                                    vertex_idx;
        Cell                                   *cell_;
        vector< long >::const_iterator          i, e = cell_list.cend();
        vector< pair< long, pair< long, short > > >::const_iterator     ii, ee = ghost_list.cend();

        // Extract vertex list from cells ----------------------------------- //
        vertex_counter = 0;
        for ( i = cell_list.cbegin(); i != e; ++i ) {
            cell_ = &m_cells[*i];
            n_vertices = cell_->getVertexCount();
            for ( j = 0; j < n_vertices; ++j ) {
                vertex_idx = cell_->getVertex( j );
                if ( vertex_map.count(vertex_idx) == 0 ) {
                    vertex_list.push_back(vertex_idx);
                    vertex_map[vertex_idx] = vertex_counter;
                    ++vertex_counter;
                }
            } //next j
        } //next i

        // Extract vertex list from ghosts ---------------------------------- //
        for ( ii = ghost_list.cbegin(); ii != ee; ++ii ) {
            cell_ = &m_cells[ii->first];
            n_vertices = cell_->getVertexCount();
            for ( j = 0; j < n_vertices; ++j ) {
                vertex_idx = cell_->getVertex( j );
                if ( vertex_map.count(vertex_idx) == 0 ) {
                    vertex_list.push_back(vertex_idx);
                    vertex_map[vertex_idx] = vertex_counter;
                    ++vertex_counter;
                }
            } //next j
        } //next i
        n_vertex = vertex_counter;

/*DEBUG*/{
/*DEBUG*/    out << "    sending " << vertex_counter;
// /*DEBUG*/    out << " vertices " << vertex_list;
/*DEBUG*/    out << endl;
// /*DEBUG*/    out << "    vertex_map: ";
// /*DEBUG*/    unordered_map<long, long>::const_iterator jj;
// /*DEBUG*/    for (jj = vertex_map.begin(); jj != vertex_map.end(); ++jj) {
// /*DEBUG*/           out << "(" << jj->first << "->" << jj->second << ") ";
// /*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Communicate vertices ================================================= //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating vertices {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        vector<long>::const_iterator            i, e = vertex_list.cend();
        
        // Initialize communication buffer ---------------------------------- //
        buff_size = 3 * n_vertex * sizeof(double) + sizeof(long);
        OBinaryStream                           com_buff( buff_size );

        // Stream coords to communication buffer ---------------------------- //
        com_buff << n_vertex;
        for ( i = vertex_list.cbegin(); i != e; ++i ) {
            com_buff << m_vertices[*i];
        } //next i

        // Send communication buffer ---------------------------------------- //
        MPI_Send( &buff_size, 1, MPI_LONG, rcv_rank, 0, m_communicator );
        MPI_Send( com_buff.get_buffer(), buff_size, MPI_CHAR, rcv_rank, 1, m_communicator );

/*DEBUG*/{
/*DEBUG*/    out << "    " << buff_size << " bytes sent" << endl;
/*DEBUG*/    out << endl;
/*DEBUG*/}
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Send cells =========================================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating cells {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                     j, k, l;
        int                                     n_vertices, n_faces, n_adj, n_itf;
        long                                    vertex_idx, neigh_idx;
        Cell                                   *cell_;
        vector<long>                            connect_;
        vector< vector<long> >                  interfs_;
        vector<long>::const_iterator            i, e;
        unordered_map<long, long>::const_iterator     m;

        // Compute buffer size required to communicate cells ---------------- //
        buff_size = sizeof(long);
        e = cell_list.cend();
        for ( i = cell_list.cbegin(); i != e; ++i ) {
            buff_size += m_cells[ *i ].getBinarySize();
        } //next i

        // Initialize communication buffer ---------------------------------- //
        OBinaryStream           com_buff( buff_size );

        // Fill communication buffer ---------------------------------------- //
        com_buff << n_cells;
        e = cell_list.cend();
        for (i = cell_list.cbegin(); i != e; ++i) {

            // Pointer to cells
            cell_ = &m_cells[ *i ];
// /*DEBUG*/   {
// /*DEBUG*/       out << "    sending cell: " << endl;
// /*DEBUG*/       cell_->display(out, 6);
// /*DEBUG*/   }

            // Modify connectivity to match sending order
            n_vertices = cell_->getVertexCount();
            connect_.resize( n_vertices );
            for ( j = 0; j < n_vertices; ++j ) {
                vertex_idx = cell_->getVertex( j );
                cell_->setVertex( j, vertex_map[ vertex_idx ] );
                connect_[ j ] = vertex_idx;
            } //next j

            // Modify adjacencies to match sending order
            n_faces = cell_->getFaceCount();
            interfs_.resize(n_faces);
            for ( j = 0; j < n_faces; ++j ) {
                n_adj = cell_->getAdjacencyCount( j );
                interfs_[ j ].resize( n_adj, Element::NULL_ID );
                l = 0;
                for ( k = 0; k < n_adj; ++k ) {
                    neigh_idx = cell_->getAdjacency( j, k );
                    interfs_[ j ][ l ] = neigh_idx;
                    ++l;
                    if ( neigh_idx != Element::NULL_ID ) {
                        if ( m_cells[ neigh_idx ].isInterior() ) {
                            if ( cell_map.find( neigh_idx ) != cell_map.end() )        cell_->setAdjacency( j, k, cell_map[ neigh_idx ] );
                            else if ( ghost_map.find( neigh_idx ) != ghost_map.end() ) cell_->setAdjacency( j, k, ghost_map[ neigh_idx ] );
                            else {
                                if (n_adj == 1) cell_->setAdjacency( j, k, Element::NULL_ID );
                                else {
                                    cell_->deleteAdjacency( j, k );
                                    --n_adj;
                                    --k;
                                }
                            }
                        }
                        else {
                            m = find( m_ghost2id[rcv_rank].begin(), m_ghost2id[rcv_rank].end(), UnaryPredicate<const long, long>(neigh_idx) );
                            if ( m != m_ghost2id[rcv_rank].end() ) cell_->setAdjacency( j, k, m->first + n_cells + n_ghosts );
                            else {
                                if (n_adj == 1) cell_->setAdjacency( j, k, Element::NULL_ID );
                                else {
                                    cell_->deleteAdjacency( j, k );
                                    --n_adj;
                                    --k;
                                }
                            }
                        }
                    }
                } //next k
            } //next j

            // Store cell into communication buffer
            com_buff << *cell_;
// /*DEUBG*/   {
// /*DEUBG*/       out << "    re-arranged as:" << endl;
// /*DEUBG*/       cell_->display(out, 6);
// /*DEUBG*/   }

            // Restore original connectivity
            for (j = 0; j < n_vertices; ++j) {
                cell_->setVertex(j, connect_[j]);
            } //next j

            // Restore original adjacencies
            for ( j = 0; j < n_faces; ++j ) {
                n_adj = cell_->getAdjacencyCount( j );
                for (k = 0; k < n_adj; ++k) {
                    cell_->setAdjacency( j, k, interfs_[j][k] );
                } //next k
                n_itf = interfs_[j].size();
                for ( k = n_adj; k < n_itf; ++k ) {
                    cell_->pushAdjacency( j, interfs_[j][k] );
                } //next k
            } //next j

        } //next i

        // Communicate cells ------------------------------------------------ //
        MPI_Send(&buff_size, 1, MPI_LONG, rcv_rank, 2, m_communicator);
        MPI_Send(com_buff.get_buffer(), buff_size, MPI_CHAR, rcv_rank, 3, m_communicator);
/*DEBUG*/out << "    " << buff_size << " bytes sent" << endl;

        // Receive new IDs -------------------------------------------------- //

        // Receive buffer size
        MPI_Recv(&buff_size, 1, MPI_LONG, rcv_rank, 6, m_communicator, MPI_STATUS_IGNORE);

        // Initialize input memory stream
        IBinaryStream                           feedback( buff_size );

        // Initialize buffer for notification
        long                                    notification_size = (2 * n_cells + 1) * sizeof(long);
        OBinaryStream                           notification( notification_size );

        // Receive new IDs
        MPI_Recv( feedback.get_buffer(), buff_size, MPI_CHAR, rcv_rank, 7, m_communicator, MPI_STATUS_IGNORE );

        // Update ghost cells
        notification << n_cells;
        e = cell_list.cend();
        for ( i = cell_list.cbegin(); i != e; ++i ) {
            feedback >> neigh_idx;
            if ( sender_ghost_new_ids.count(*i) > 0 ) {
                sender_ghost_new_ids[*i] = neigh_idx;
            }
            //m_ghost2id[rcv_rank][neigh_idx] = *i;
            notification << *i << neigh_idx;
        } //next i
/*DEBUG*/{
// /*DEBUG*/    out << "  receiving new IDs, ";
// /*DEBUG*/    out << "  ghost2idx[" << rcv_rank << "] = {";
// /*DEBUG*/    for (m = m_ghost2id[rcv_rank].begin(); m != m_ghost2id[rcv_rank].end(); ++m) {
// /*DEBUG*/        out << " (" << m->first << ", " << m->second << ")";
// /*DEBUG*/    } //next m
// /*DEBUG*/    out << " }" << endl << endl;
/*DEBUG*/}

        // Notify other processes ------------------------------------------- //
        int                     nproc;
        MPI_Comm_size(m_communicator, &nproc);
        for (j = 0; j < nproc; ++j) {
            if ( (j != snd_rank) && (j != rcv_rank) && (m_ghost2id.find(j) != m_ghost2id.end()) ) {
                MPI_Send(&notification_size, 1, MPI_LONG, j, 8+j, m_communicator);
                MPI_Send(notification.get_buffer(), notification_size, MPI_CHAR, j, 8+nproc+j, m_communicator);
            }
        } //next j
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Send ghost cells ===================================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating ghosts {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                     j, k;
        int                                     n_vertices, n_faces, n_adj, n_itf;
        long                                    vertex_idx, neigh_idx, cell_idx;
        vector<long>                            connect_;
        vector< vector<long> >                  interfs_;
        Cell                                   *cell_;
        vector< pair< long, pair< long, short > > >::const_iterator             i, e;
        unordered_map<long, long>::const_iterator                               m;

        // Compute buffer size ---------------------------------------------- //
        buff_size = sizeof(long);
        e = ghost_list.cend();
        for ( i = ghost_list.cbegin(); i != e; ++i ) {
            buff_size += m_cells[ i->first ].getBinarySize() + sizeof(short);
        } //next i

        // Initialize communication buffer ---------------------------------- //
        OBinaryStream                   com_buff(buff_size);

        // Fill communication buffer ---------------------------------------- //
        com_buff << n_ghosts;
        e = ghost_list.cend();
        for ( i = ghost_list.cbegin(); i != e; ++i ) {

            // Pointer to cell
            cell_ = &m_cells[ i->first ];

            // Modify cell ID in case of ghost on another process
            if ( i->second.second != snd_rank ) {
                cell_idx = cell_->getId();
                cell_->setId( i->second.first);
            }
/*DEBUG*/   {
// /*DEBUG*/       out << "    sending ghost:" << endl;
// /*DEBUG*/       cell_->display(out, 6);
/*DEBUG*/   }

            // Modify connectivity to match sending order
            n_vertices = cell_->getVertexCount();
            connect_.resize(n_vertices);
            for ( j = 0; j < n_vertices; ++j ) {
                vertex_idx = cell_->getVertex( j );
                connect_[ j ] = vertex_idx;
                cell_->setVertex( j, vertex_map[ vertex_idx ] );
            } // next j

            // Modify adjacencies to match sending order
            n_faces = cell_->getVertexCount();
            interfs_.resize( n_faces );
            for ( j = 0; j < n_faces; ++j ) {
                n_adj = cell_->getAdjacencyCount( j );
                interfs_[j].resize( n_adj, Element::NULL_ID );
                for ( k = 0; k < n_adj; ++k ) {
                    neigh_idx = cell_->getAdjacency( j, k );
                    interfs_[ j ][ k ] = neigh_idx;
                    if (neigh_idx != Element::NULL_ID) {
                        if ( m_cells[ neigh_idx ].isInterior() ) {
                            if ( cell_map.find( neigh_idx ) != cell_map.end() )        cell_->setAdjacency( j, k, cell_map[ neigh_idx ] );
                            else if ( ghost_map.find( neigh_idx ) != ghost_map.end() ) cell_->setAdjacency( j, k, ghost_map[ neigh_idx ] );
                            else {
                                if (n_adj == 1) cell_->setAdjacency( j, k, Element::NULL_ID );
                                else {
                                    cell_->deleteAdjacency( j, k );
                                    --n_adj;
                                    --k;
                                }
                            }
                        }
                        else {
                            m = find( m_ghost2id[rcv_rank].begin(), m_ghost2id[rcv_rank].end(), UnaryPredicate<const long, long>(neigh_idx) );
                            if ( m != m_ghost2id[rcv_rank].end() ) cell_->setAdjacency( j, k, m->first + n_cells + n_ghosts );
                            else {
                                if (n_adj == 1) cell_->setAdjacency( j, k, Element::NULL_ID );
                                else {
                                    cell_->deleteAdjacency( j, k );
                                    --n_adj;
                                    --k;
                                }
                            }
                        }
                    }
                } //next k
            } //next j

            // Fill communication buffer
            com_buff << i->second.second;
            com_buff << *cell_;
/*DEBUG*/   {
// /*DEBUG*/       out << "    remapped as:" << endl;
// /*DEBUG*/       cell_->display(out, 6);
/*DEBUG*/   }

            // Restore original id
            if ( i->second.second != snd_rank ) {
                cell_->setId( cell_idx );
            }

            // Restore original connectivity
            for ( j = 0; j < n_vertices; ++j ) {
                cell_->setVertex( j, connect_[j] );
            } //next j

            // Restore original adjacencies
            for ( j = 0; j < n_faces; ++j ) {
                n_adj = cell_->getAdjacencyCount( j );
                for ( k = 0; k < n_adj; ++k ) {
                    cell_->setAdjacency( j, k, interfs_[j][k] );
                } //next k
                n_itf = interfs_[j].size();
                for ( k = n_adj; k < n_itf; ++k ) {
                    cell_->pushAdjacency( j, interfs_[j][k] );
                }
            } //next j
        } //next i

        // Communicate ghosts ----------------------------------------------- //
        MPI_Send(&buff_size, 1, MPI_LONG, rcv_rank, 4, m_communicator);
        MPI_Send(com_buff.get_buffer(), buff_size, MPI_CHAR, rcv_rank, 5, m_communicator);
/*DEBUG*/out << "    " << buff_size << " bytes sent" << endl << endl;
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Update ghost lists =================================================== //
    // m_ghost2id stores the list of ghost divided by rank, that is:
    // if the cell with ID "i" is a ghost cell on the current rank, the cell 
    // exists on process with rank "neigh_rank" and hase ID "j" in that partition,
    // thus m_ghost2id[neigh_rank][j] stores "i".
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), updating ghosts lists {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        bool                                                    flag_keep;
        int                                                     j;
        int                                                     n_neighs;
        long                                                    ghost_idx;
        long                                                    ghost_send_idx, ghost_recv_idx;
        vector<long>                                            neighs;
        unordered_map<long, long>::iterator                     ii, ee;
        unordered_map<long, long>::const_iterator               i, e;
        vector<long>::const_iterator                            m, n;

        // Initialize sender's ghosts --------------------------------------- //
        e = sender_ghost_new_ids.cend();
        for (i = sender_ghost_new_ids.cbegin(); i != e; ++i) {
            ghost_send_idx = i->first;
            ghost_recv_idx = i->second;

            moveInternal2Ghost(ghost_send_idx);
            m_ghost2id[rcv_rank][ghost_recv_idx] = ghost_send_idx;
        }

        // Delete sent cells that are not ghosts -----------------------------//
        for (const auto &cell : cell_list) {
            if (m_cells[cell].isInterior()) {
                deleteCell(cell, true, true);
            }
        }

        // Delete previous ghost cells of other processors -------------------//
        for (auto &rank_ghost : m_ghost2id) {
            ee = rank_ghost.second.end();
            ii = rank_ghost.second.begin();
            while (ii != ee) {
                ghost_idx = ii->second;
                neighs = findCellNeighs(ghost_idx);
                n_neighs = neighs.size();

                // We want to keep only ghosts sells that have at least one
                // internal neighbour
                flag_keep = false;
                for (j = 0; j < n_neighs; ++j) {
                    if ( m_cells[neighs[j]].isInterior() ) {
                        flag_keep = true;
                        break;
                    }
                } //next j

                if ( !flag_keep )  {
// /*DEBUG*/           out << "      deleting ghost: " << ghost_idx << endl;
                    deleteCell(ghost_idx, true, true);
                    ii = rank_ghost.second.erase(ii);
                }
                else {
                    ++ii;
                }
            }
        } //next cell

        // Flush cells -------------------------------------------------------//
        m_cells.flush();

/*DEBUG*/{        
// /*DEBUG*/    out << "  ghost2idx = {";
// /*DEBUG*/    for (i = m_ghost2id[rcv_rank].begin(); i != m_ghost2id[snd_rank].end(); ++i) {
// /*DEBUG*/        out << " (" << i->first << ", " << i->second << ")";
// /*DEBUG*/    }
// /*DEBUG*/    out << " }" << endl << endl;
// /*DEBUG*/    out << "  nCells: " << m_nInternals << endl;
// /*DEBUG*/    out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/}
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Remove isolated vertices ============================================= //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), removing isolated vertices {" << endl << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        // none

        // Remove Isolated vertices ----------------------------------------- //
        deleteOrphanVertices();
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

/*DEBUG*/{
/*DEBUG*/    ostream        *msg = reinterpret_cast<ostream*>(&out);
/*DEBUG*/    out << "* stats (rank#" << snd_rank << "): " << endl;
/*DEBUG*/    displayTopologyStats(*msg);
/*DEBUG*/    out << endl;
/*DEBUG*/    out << "* sender (rank #" << snd_rank << ") completed its tasks" << endl;
/*DEBUG*/}
}

// ========================================================================== //
// RECEIVER                                                                   //
// ========================================================================== //
if (m_rank == rcv_rank)
{

    // Scope variables ====================================================== //

    // Variables required for communication
    long                                        buff_size;

    // Variables required for vertex communication
    long                                        n_vertex;
    vector<long>                                v_local_mapping;

    // Variables required for cells communication
    long                                        n_cells;
    unordered_map<long, long>                   c_local_mapping;

    // Variables required for ghosts communication
    long                                        n_ghosts;

    // Receive vertices ===================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving vertices {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        long                                    i;
        Vertex                                  vertex;
        VertexIterator                          it;

        // Receive communication buffer ------------------------------------- //

        // Receive buffer size
        MPI_Recv(&buff_size, 1, MPI_LONG, snd_rank, 0, m_communicator, MPI_STATUS_IGNORE); 
/*DEBUG*/out << "    receiving " << buff_size << " bytes" << endl;

        // Initialize memory stream
        IBinaryStream                   com_buff( buff_size );

        // Update vertex list ----------------------------------------------- //

        // Receive vertex coordinate list
        MPI_Recv( com_buff.get_buffer(), buff_size, MPI_CHAR, snd_rank, 1, m_communicator, MPI_STATUS_IGNORE);

        // Read number of vertices stored in buffer
        com_buff >> n_vertex;

        // Insert vertices
        v_local_mapping.resize(n_vertex, 0);
        for (i = 0; i < n_vertex; ++i) {
            com_buff >> vertex;
            it = addVertex( move(vertex) );
            v_local_mapping[i] = it->getId();
        } //next i
/*DEBUG*/out << "    received " << n_vertex;
// /*DEBUG*/out << " vertices " << v_local_mapping;
/*DEBUG*/out << endl;
        
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Receive cells ======================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving cells {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                             j;
        int                                             n_vertices;
        long                                            i, cell_idx, ghost_idx;
        long                                            cell_count = 0;
        long                                            feedback_size;
        Cell                                            cell;
        CellIterator                                    it;

        // Receive communication buffer ------------------------------------- //

        // Receive size of communication buffer
        MPI_Recv(&buff_size, 1, MPI_LONG, snd_rank, 2, m_communicator, MPI_STATUS_IGNORE);
/*DEBUG*/out << "    receiving " << buff_size << " bytes" << endl;

        // Initialize memory stream
        IBinaryStream                           com_buff(buff_size);

        // Update cell list ------------------------------------------------- //

        // Receive cells
        MPI_Recv(com_buff.get_buffer(), buff_size, MPI_CHAR, snd_rank, 3, m_communicator, MPI_STATUS_IGNORE);

        // Extract cells from communication buffer
        com_buff >> n_cells;

        // Initialize communication buffer for feedback
        feedback_size = n_cells * sizeof(long);
        OBinaryStream                           feedback( feedback_size );

        // Add cells to cell list
        reserveCells( m_nInternals + m_nGhosts + n_cells );
        for ( i = 0; i < n_cells; ++i ) {

            // Stream cell data
            com_buff >> cell;
            cell_idx = cell.getId();
            cell.setInterior( true );

/*DEBUG*/   {
// /*DEBUG*/       out << "    received cell:" << endl;
// /*DEBUG*/       cell.display(out, 6);
/*DEBUG*/   }

            // Update connectivity
            n_vertices = cell.getVertexCount();
            for (j = 0; j < n_vertices; ++j) {
                cell.setVertex( j, v_local_mapping[ cell.getVertex( j ) ] );
            } //next j

            // Insert cell in cells list & update ghost lists
            if ( m_ghost2id[snd_rank].find(cell_idx) != m_ghost2id[snd_rank].end() ) {
                ghost_idx = m_ghost2id[snd_rank][cell_idx];
                cell.setId( ghost_idx );
                feedback << ghost_idx;
/*DEBUG*/       {
// /*DEBUG*/           out << "    (already exists as ghost), remapped as:" << endl;
// /*DEBUG*/           cell.display(out, 6);
/*DEBUG*/       }
                m_cells[ghost_idx] = move( cell );
                moveGhost2Internal(ghost_idx);
                m_ghost2id[snd_rank].erase( cell_idx );
                c_local_mapping[ cell_count ] = ghost_idx;
                ++cell_count;
            }
            else {

                // Add cell before last internal cell
// /*DEBUG*/       out << "    (adding new cell), remapped as:" << endl;
                it = addCell( move(cell), generateCellId() );
                feedback << long( it->getId() );
// /*DEBUG*/       it->display(out, 6);
                c_local_mapping[ cell_count ] = it->getId();
                ++cell_count;
            }
        } //next i

/*DEBUG*/{
/*DEBUG*/    out << "    received " << cell_count << " cells: ";
// /*DEBUG*/    unordered_map<long, long>::const_iterator    ii;
// /*DEBUG*/    for (ii = c_local_mapping.begin(); ii != c_local_mapping.end(); ++ii) {
// /*DEBUG*/           out << "(" << ii->first << ", " << ii->second << "), ";
// /*DEBUG*/    } //next i
/*DEBUG*/    out << endl;
/*DEBUG*/}

        // Communicate new IDs to sender ------------------------------------ //
/*DEBUG*/{
// /*DEBUG*/    out << "    sending new IDs: ";
// /*DEBUG*/    unordered_map<long, long>::const_iterator      m;
// /*DEBUG*/    for (m = c_local_mapping.begin(); m != c_local_mapping.end(); ++m) {
// /*DEBUG*/        out << " (" << m->first << ", " << m->second << ") ";
// /*DEBUG*/    } //next m
// /*DEBUG*/    out << endl << endl;
/*DEBUG*/}

        // Communicate buffer size
        MPI_Send( &feedback_size, 1, MPI_LONG, snd_rank, 6, m_communicator );

        // Communicate buffer
        MPI_Send( feedback.get_buffer(), feedback_size, MPI_CHAR, snd_rank, 7, m_communicator );

    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Receive ghosts ======================================================= //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving ghosts {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {

        // Scope variables -------------------------------------------------- //
        short                                   ghost_rank;
        int                                     j;
        int                                     n_vertices;
        long                                    i;
        long                                    cell_idx;
        long                                    ghost_count = 0;
        Cell                                    cell;
        CellIterator                            it;

        // Receive communication buffer ------------------------------------- //

        // Receive size of communication buffer
        MPI_Recv(&buff_size, 1, MPI_LONG, snd_rank, 4, m_communicator, MPI_STATUS_IGNORE);
/*DEBUG*/out << "    receiving " << buff_size << " bytes" << endl;

        // Initialize memory stream
        IBinaryStream                  com_buff(buff_size);

        // Update cells list ------------------------------------------------ //

        // Receive cell list
        MPI_Recv(com_buff.get_buffer(), buff_size, MPI_CHAR, snd_rank, 5, m_communicator, MPI_STATUS_IGNORE);

        // Extract cells from communication buffer
        com_buff >> n_ghosts;

        // Add cells to cell list
        reserveCells(m_nInternals + m_nGhosts + n_ghosts);

        for (i = 0; i < n_ghosts; ++i) {

            // Stream cell data 
            com_buff >> ghost_rank;
            com_buff >> cell;
            cell_idx = cell.getId();
/*DEBUG*/   {
// /*DEBUG*/       out << "    receiving ghost (existing on rank: " << ghost_rank << ")" << endl;
// /*DEBUG*/       cell.display(out, 6);
/*DEBUG*/   }

            // Update connectivity
            n_vertices = cell.getVertexCount();
            for ( j = 0; j < n_vertices; ++j ) {
                cell.setVertex( j, v_local_mapping[ cell.getVertex( j ) ] );
            } //next j

            // Update cell list
            if ( m_ghost2id[ghost_rank].find(cell_idx) == m_ghost2id[ghost_rank].end() ) {

// /*DEBUG*/       out << "    (non-existent ghost cell)" << endl;

                // Cell has to be added after the last internal cell
                cell.setInterior( false );
                it = addCell( move(cell), generateCellId() );

                c_local_mapping[ n_cells + ghost_count ] = it->getId();
                m_ghost2id[ghost_rank][ cell_idx ] = it->getId();
                ++ghost_count;

/*DEBUG*/       {
// /*DEBUG*/           out << "    remapped as:" << endl;
// /*DEBUG*/           it->display(out, 6);
/*DEBUG*/       }
            }
            else {
// /*DEBUG*/       out << "    (cell already exists as ghost)" << endl;
                cell.setId(m_ghost2id[ghost_rank][cell_idx]);
                cell.setInterior( false );
/*DEBUG*/       {
// /*DEBUG*/           out << "    remapped as:" << endl;
// /*DEBUG*/           cell.display(out, 6);
/*DEBUG*/       }

                // Cell has to be overwritten
                m_cells[m_ghost2id[ghost_rank][cell_idx]] = move(cell);

                c_local_mapping[ n_cells + ghost_count ] = m_ghost2id[ghost_rank][cell_idx];
                ++ghost_count;
            }

        } //next i

/*DEBUG*/{
// /*DEBUG*/    unordered_map<long, long>::const_iterator      m;
// /*DEBUG*/    out << "    c_local_mapping is now: {";
// /*DEBUG*/    for ( m = c_local_mapping.begin(); m != c_local_mapping.end(); ++m) {
// /*DEBUG*/        out << " (" << m->first << ", " << m->second << ")";
// /*DEBUG*/    } //next m
// /*DEBUG*/    out << " }" << endl << endl;
/*DEBUG*/}
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Update adjacencies =================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), updating adjacencies {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        int                                             n_faces, n_adj;
        int                                             j, k;
        long                                            cell_idx, neigh_idx;
        unordered_map<long, long>::const_iterator       m, e;
        Cell                                           *cell_;

        // Update adjacencies ----------------------------------------------- //
        e = c_local_mapping.cend();
        for ( m = c_local_mapping.cbegin(); m != e; ++m ) {
            cell_idx = m->second;
            cell_ = &m_cells[cell_idx];
// /*DEBUG*/   out << "    updating adjacencies for cell ID " << cell_idx << endl;
            n_faces = cell_->getFaceCount();
            for ( j = 0; j < n_faces; ++j ) {
                n_adj = cell_->getAdjacencyCount( j );
                for ( k = 0; k < n_adj; ++k ) {
                    neigh_idx = cell_->getAdjacency( j, k );
                    if ( neigh_idx >= 0) {
                        if ( neigh_idx >= n_cells + n_ghosts ) cell_->setAdjacency( j, k, neigh_idx - n_cells - n_ghosts );
                        else                                   cell_->setAdjacency( j, k, c_local_mapping[neigh_idx] );
                    }
                } //next k
            } //next j
// /*DEBUG*/   out << "    cell ID " << cell_idx << " is now" << endl;
// /*DEBUG*/   cell_->display(out, 6);
        } //next m

// /*DEBUG*/out << "  ghost2idx = {";
// /*DEBUG*/for (m = m_ghost2id[snd_rank].begin(); m != m_ghost2id[snd_rank].end(); ++m) {
// /*DEBUG*/    out << " (" << m->first << ", " << m->second << ")";
// /*DEBUG*/}
// /*DEBUG*/out << " }" << endl << endl;

    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Remove duplicated vertices =========================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), removing duplicated vertices {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        // none

        // Remove duplicated vertices --------------------------------------- //
        out << "    n. vertices is: " << getVertexCount() << endl;
        deleteCoincidentVertices();
        deleteOrphanVertices();
        out << "    (after cleaning), n.vertices is: " << getVertexCount() << endl << endl;
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

    // Update adjacencies =================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), updating adjacencies {" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    {
        // Scope variables -------------------------------------------------- //
        vector<long>                                    cell_list(n_cells + n_ghosts, -1);
        unordered_map<long, long>::const_iterator       i, e;
        long                                            j;

        // Update Adjacencies ----------------------------------------------- //
        j = 0;
        e = c_local_mapping.cend();
        for (i = c_local_mapping.cbegin(); i != e; ++i) {
            cell_list[j] = i->second;
            ++j;
        } //next i

        // Update Adjacencies ----------------------------------------------- //
// /*DEBUG*/out << "    updating adjacencies of simplicies: " << cell_list << endl << endl;
        updateAdjacencies(cell_list);
/*DEBUG*/{
// /*DEBUG*/    for (int ii = 0; ii < cell_list.size(); ++ii) {
// /*DEBUG*/           out << "    cell ID: " << cell_list[ii] << " is now: " << endl;
// /*DEBUG*/           m_cells[cell_list[ii]].display(out, 4);
// /*DEBUG*/    } //next ii
/*DEBUG*/}
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;

/*DEBUG*/{
/*DEBUG*/    ostream        *msg = reinterpret_cast<ostream*>(&out);
/*DEBUG*/    out << "* stats (rank#" << rcv_rank << "): " << endl;
/*DEBUG*/    displayTopologyStats(*msg);
/*DEBUG*/    out << endl;
/*DEBUG*/    out << "* receiver (rank #" << rcv_rank << ") completed its tasks" << endl << endl;
/*DEBUG*/}

}

// ========================================================================== //
// NOTIFY OTHER PROCESS (MULTIPLE ADJACENCIES)                                //
// ========================================================================== //
if ( (m_rank != snd_rank) && (m_rank != rcv_rank) )
{

    // Scope variables ====================================================== //
    bool                                                                waiting;
    int                                                                 nproc;
    long                                                                buff_size;
    long                                                                n_cells;
    long                                                                i;
    long                                                                cell_idx, new_cell_idx, my_cell_idx;
    IBinaryStream                                                       buffer;
    unordered_map< short, unordered_map<long, long> >::const_iterator   m, e;

    // Check if something has to be received ================================ //
    MPI_Comm_size(m_communicator, &nproc);
    waiting = false;
    e = m_ghost2id.cend();
    for (m = m_ghost2id.cbegin(); m != e; ++m) {
        waiting |= (m->first == snd_rank);
    };

    // Receive cell list ==================================================== //
/*DEBUG*/out << "* process (rank#" << m_rank << ") receving notification" << endl;
/*DEBUG*/t0 = high_resolution_clock::now();
    if (waiting) {

        // Receive buffer size
        MPI_Recv(&buff_size, 1, MPI_LONG, snd_rank, 8 + m_rank, m_communicator, MPI_STATUS_IGNORE);

        // Initialize stream from memory
        IBinaryStream           buffer( buff_size );

        // Receive buffer
        MPI_Recv(buffer.get_buffer(), buff_size, MPI_CHAR, snd_rank, 8 + nproc + m_rank, m_communicator, MPI_STATUS_IGNORE);

        // Update ghost list
/*DEBUG*/{
/*DEBUG*/    out << "* process (rank#" << m_rank << ") updating ghosts lists" << endl;
/*DEBUG*/    out << "* process (rank#" << m_rank << ") is waiting for notification" << endl;
/*DEBUG*/}
        if ( m_ghost2id.find( snd_rank ) != e ) {
            buffer >> n_cells;
/*DEGUB*/   out << "  number of cells is " << n_cells << endl;
            for (i = 0; i < n_cells; ++i) {
                buffer >> cell_idx;
                buffer >> new_cell_idx;
                out << "    being notifyied that cell " << cell_idx << " has now ID: " << new_cell_idx << " on rank " << rcv_rank << endl;
                if ( m_ghost2id[snd_rank].find( cell_idx ) != m_ghost2id[snd_rank].end() ) {
                    my_cell_idx = m_ghost2id[snd_rank][cell_idx];
                    m_ghost2id[snd_rank].erase( cell_idx );
                    m_ghost2id[rcv_rank][ new_cell_idx ] = my_cell_idx;
                }
            } //next i
        }
    }
/*DEBUG*/t1 = high_resolution_clock::now();
/*DEBUG*/time_span = duration_cast<duration<double>>(t1 - t0);
/*DEBUG*/out<< "}" << endl << "  (" << time_span.count() << " sec.)" << endl;
}

/*DEBUG*/{
/*DEBUG*/    out << "* display mesh infos {" << endl;
/*DEBUG*/    displayCells(out);
/*DEBUG*/    out << "}" << endl;
/*DEBUG*/}
/*DEBUG*/{
/*DEBUG*/    out << "* display ghost infos {" << endl;
/*DEBUG*/    unordered_map<short, unordered_map<long, long>>::const_iterator    ii;
/*DEBUG*/    unordered_map<long, long>::const_iterator                          jj;
/*DEBUG*/    for (ii = m_ghost2id.cbegin(); ii != m_ghost2id.cend(); ++ii) {
/*DEBUG*/        out << "  rank: " << ii->first << endl;
/*DEBUG*/        for (jj = ii->second.cbegin(); jj != ii->second.cend(); ++jj) {
/*DEBUG*/            out << "    {" << jj->second << ", " << jj->first << "}" << endl;
/*DEBUG*/        } //next j
/*DEBUG*/    } //next ii
/*DEBUG*/    out << "}" << endl;
/*DEBUG*/    out.close();
/*DEBUG*/}

return; }

/*!
	@}
*/

}

#endif
