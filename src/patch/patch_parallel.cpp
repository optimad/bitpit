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
#if ENABLE_MPI==1

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include <mpi.h>
#include "patch.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

namespace bitpit {

/*!
	\ingroup patch
	@{

	\class Patch
*/

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void Patch::setCommunicator(MPI_Comm communicator)
{
	// Free previous communicator
	if (m_communicator != MPI_COMM_NULL) {
		MPI_Comm_free(&m_communicator);
		m_communicator = MPI_COMM_NULL;
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	if (communicator != MPI_COMM_NULL) {
		MPI_Comm_dup(communicator, &m_communicator);
	}

	// Get MPI information
	if (m_communicator) {
		MPI_Comm_size(m_communicator, &m_nProcessors);
		MPI_Comm_rank(m_communicator, &m_rank);
	} else {
		m_rank        = 0;
		m_nProcessors = 1;
	}

	// Set parallel data for the VTK output
	setParallel(m_nProcessors, m_rank);
}

/*!
	Unsets the MPI communicator to be used for parallel communications.
*/
void Patch::unsetCommunicator()
{
	setCommunicator(MPI_COMM_NULL);
}

/*!
	Gets the MPI communicator associated to the patch

	\return The MPI communicator associated to the patch.
*/
const MPI_Comm & Patch::getCommunicator() const
{
	return m_communicator;
}

/*!
	Gets the MPI rank associated to the patch

	\return The MPI rank associated to the patch.
*/
int Patch::getRank() const
{
	return m_rank;
}

/*!
	Gets the MPI processors in the communicator associated to the patch

	\return The MPI processors in the communicator associated to the patch
*/
int Patch::getProcessorCount() const
{
	return m_nProcessors;
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
void Patch::sendCells(const unsigned short &snd_rank, const unsigned short &rcv_rank, const vector< long > &cell_list)
{

// ========================================================================== //
// SCOPE VARIABLES                                                            //
// ========================================================================== //

// Debug variables
/*DEBUG*/stringstream                           out_name;
/*DEBUG*/ofstream                               out;

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

    // Variables required for vertex communication
    long                                        n_vertex;
    vector<long>                                vertex_list;
    unordered_map<long, long>                   vertex_map;
    
    // Initialize data structures =========================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), initializing data structure" << endl;
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
/*DEBUG*/    out << "    sending " << n_cells << " cell(s): " << cell_list << endl;
/*DEBUG*/    out << "    cell_map: ";
/*DEBUG*/    unordered_map<long, long>::const_iterator      ii;
/*DEBUG*/    for (ii = cell_map.begin(); ii != cell_map.end(); ++ii) {
/*DEBUG*/        out << "(" << ii->first << "->" << ii->second << ") ";
/*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

        // Initialize lists for ghost communication ------------------------- //
        ghost_list.reserve( cell_list.size() );

        // Initialize lists for vertex communication ------------------------ //
        vertex_list.reserve( 4*cell_list.size() );

    }

    // Create list of ghosts ================================================ //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), creating list of ghosts" << endl;
    {
        // Scope variables -------------------------------------------------- //
        int                                     n_vertex, n_neighs;
        int                                     k;
        vector<long>                            neighs;
        long                                    cell_idx;
        long                                    ghost_counter;
        unordered_map<long, bool>               is_ghost;
//         Cell                                    *cell;
        vector<long>::const_iterator            i, e;
        unordered_map<short, unordered_map<long, long> >::const_iterator   m;
        unordered_map<long, long>::const_iterator                          n;
        unordered_map<long, bool>::const_iterator                          l;

        // Extract ghosts --------------------------------------------------- //
        ghost_counter = 0;
        e = cell_list.cend();
        for (i = cell_list.cbegin(); i < e; ++i) {

//             cell_idx = *i;
//             cell = &cells[cell_idx];

            neighs = findCellNeighs(*i);
            n_neighs = neighs.size();
            for (k = 0; k < n_neighs; ++k) {
                is_ghost[neighs[k]] = true;
            } //next k
/*TO_BE_REMOVED*/{
/*TO_BE_REMOVED*///            n_vertex = cell->getVertexCount();
/*TO_BE_REMOVED*///            for ( j = 0; j < n_vertex; ++j ) {
/*TO_BE_REMOVED*///                ring_1 = Ring1(cell_idx, j);
/*TO_BE_REMOVED*///                n_neighs = ring_1.size();
/*TO_BE_REMOVED*///                for ( k = 0; k < n_neighs; ++k ) {
/*TO_BE_REMOVED*///                    is_ghost[ring_1[k]] = true;
/*TO_BE_REMOVED*///                } //next k
/*TO_BE_REMOVED*///            } //next j
/*TO_BE_REMOVED*/}
        } //next i
        for (i = cell_list.begin(); i < cell_list.end(); ++i) {
            cell_idx = *i;
            is_ghost[cell_idx] = false;
        } //next i

        for (l = is_ghost.begin(); l != is_ghost.end(); ++l) {
            if (l->second) {
                cell_idx = l->first;
                if ( ( find(m_ghost2id[rcv_rank].begin(),
                            m_ghost2id[rcv_rank].end(),
                            UnaryPredicate<const long, long>(cell_idx) ) == m_ghost2id[rcv_rank].end() )   // check #1: cell in not a ghost for receiver
                && ( ghost_map.find(cell_idx) == ghost_map.end() )                                     // check #2: cell has not already inserted in the ghost list
                && ( cell_map.find(cell_idx) == cell_map.end() ) ) {                                   // check #3: cell has not already inserted in the cell list
                 
/*TMP*/             {
/*TMP*/                 cout << "cells[cell_idx].isInterior(): " << m_cells[cell_idx].isInterior() << endl;
/*TMP*/                 cout << "cell_map.find( neigh_idx ) == cell_map.end(): " << bool( cell_map.find(cell_idx) == cell_map.end() ) << endl;
/*TMP*/                 cout << "ghost_map.find( neigh_idx ) == ghost_map.end(): " << bool( ghost_map.find(cell_idx) == ghost_map.end() ) << endl;
/*TMP*/             }

                    if ( m_cells[cell_idx].isInterior() ) {
                        ghost_list.push_back( pair<long, pair<long, short> >(cell_idx, pair<long, short>(cell_idx, snd_rank) ) );
                        ghost_map[cell_idx] = ghost_counter + n_cells;
                        ++ghost_counter;
                    }
                    else {
/*TMP*/                 cout << " we should never enter here" << endl;
                        for (m = m_ghost2id.begin(); m != m_ghost2id.end(); ++m) {
                            n = find( m->second.begin(), m->second.end(), UnaryPredicate<const long, long>(cell_idx) );
                            if ( n != m->second.end() ) {
                                ghost_list.push_back( pair<long, pair<long, short> >(cell_idx, pair<long, short>(n->first, m->first) ) );
                                ghost_map[cell_idx] = ghost_counter + n_cells;
                                ++ghost_counter;
                            }
                        } //next m
                    }
                }
            }
        }
        n_ghosts = ghost_counter;
/*DEBUG*/{
/*DEBUG*/    out << "    sending " << n_ghosts << " ghost(s): ";
/*DEBUG*/    vector< pair< long, pair<long, short> > >::iterator    kk;
/*DEBUG*/    for (kk = ghost_list.begin(); kk != ghost_list.end(); ++kk) {
/*DEBUG*/        out << "[" << kk->first << ", (" << kk->second.first << ", " << kk->second.second << ")], ";
/*DEBUG*/    } //next kk
/*DEBUG*/    out << endl;
/*DEBUG*/    out << "    ghost_map: ";
/*DEBUG*/    unordered_map<long, long>::const_iterator ii;
/*DEBUG*/    for (ii = ghost_map.begin(); ii != ghost_map.end(); ++ii) {
/*DEBUG*/           out << "(" << ii->first << "->" << ii->second << ") ";
/*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}
    }

    // Create list of vertices ============================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), creating list of vertices" << endl;
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
/*DEBUG*/    out << "    sending " << vertex_counter << " vertices " << vertex_list << endl;
/*DEBUG*/    out << "    vertex_map: ";
/*DEBUG*/    unordered_map<long, long>::const_iterator jj;
/*DEBUG*/    for (jj = vertex_map.begin(); jj != vertex_map.end(); ++jj) {
/*DEBUG*/           out << "(" << jj->first << "->" << jj->second << ") ";
/*DEBUG*/    }
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

    }

    // Communicate vertices ================================================= //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating vertices" << endl;
    {
        // Scope variables -------------------------------------------------- //
        int                                     j;
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

    // Send cells =========================================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating cells" << endl;
    {
        // Scope variables -------------------------------------------------- //
        int                                     j, k, l;
        int                                     n_vertices, n_faces, n_adj;
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
/*DEBUG*/   {
/*DEBUG*/       out << "    sending cell: " << endl;
/*DEBUG*/       cell_->display(out, 6);
/*DEBUG*/   }

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
/*DEUBG*/   {
/*DEUBG*/       out << "    re-arranged as:" << endl;
/*DEUBG*/       cell_->display(out, 6);
/*DEUBG*/   }

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
                for ( k = n_adj; k < interfs_[j].size(); ++k ) {
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
            m_ghost2id[rcv_rank][neigh_idx] = *i;
            notification << *i << neigh_idx;
        } //next i
/*DEBUG*/{
/*DEBUG*/    out << "  receiving new IDs, ";
/*DEBUG*/    out << "  ghost2idx[" << rcv_rank << "] = {";
/*DEBUG*/    for (m = m_ghost2id[rcv_rank].begin(); m != m_ghost2id[rcv_rank].end(); ++m) {
/*DEBUG*/        out << " (" << m->first << ", " << m->second << ")";
/*DEBUG*/    } //next m
/*DEBUG*/    out << " }" << endl << endl;
/*DEBUG*/}

        // Notify other processes ------------------------------------------- //
        int                             nproc;
        MPI_Comm_size(m_communicator, &nproc);
        for (j = 0; j < nproc; ++j) {
            if ( (j != snd_rank) && (j != rcv_rank) && (m_ghost2id.find(j) != m_ghost2id.end()) ) {
                MPI_Send(&notification_size, 1, MPI_LONG, j, 8+j, m_communicator);
                MPI_Send(notification.get_buffer(), notification_size, MPI_CHAR, j, 8+nproc+j, m_communicator);
            }
        } //next j
    }

    // Send ghost cells ===================================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), communicating ghosts" << endl;
    {
        // Scope variables -------------------------------------------------- //
        int                                     j, k;
        int                                     n_vertices, n_faces, n_adj;
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
                cell_idx = cell_->get_id();
                cell_->set_id( i->second.first);
            }
/*DEBUG*/   {
/*DEBUG*/       out << "    sending ghost:" << endl;
/*DEBUG*/       cell_->display(out, 6);
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
/*DEBUG*/       out << "    remapped as:" << endl;
/*DEBUG*/       cell_->display(out, 6);
/*DEBUG*/   }

            // Restore original id
            if ( i->second.second != snd_rank ) {
                cell_->set_id( cell_idx );
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
                for ( k = n_adj; k < interfs_[j].size(); ++k ) {
                    cell_->pushAdjacency( j, interfs_[j][k] );
                }
            } //next j
        } //next i

        // Communicate ghosts ----------------------------------------------- //
        MPI_Send(&buff_size, 1, MPI_LONG, rcv_rank, 4, m_communicator);
        MPI_Send(com_buff.get_buffer(), buff_size, MPI_CHAR, rcv_rank, 5, m_communicator);
/*DEBUG*/out << "    " << buff_size << " bytes sent" << endl << endl;
    }

    // Update ghost lists =================================================== //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), updating ghosts lists" << endl;
    {
        // Scope variables -------------------------------------------------- //
        bool                                            flag_delete;
        int                                             j, k;
        int                                             rank_id;
        int                                             n_neighs, n_adj;
        long                                            ghost_idx, neigh_idx;
        vector<long>                                    neighs;
        Cell                                           *cell_;
        unordered_map<short, unordered_map<long, long> >::const_iterator        ii, ee;
        unordered_map<long, long>::const_iterator       i, e;
        vector<long>::const_iterator                    m, n;

        // Update ghost list ------------------------------------------------ //
        n = cell_list.cend();
        for ( m = cell_list.cbegin(); m != n; ++m) {
/*DEBUG*/out << "    moving internal cell " << *m << " to ghosts" << endl;
            moveInternal2Ghost(*m);
        } //next m
            
        // Loop over ghost -------------------------------------------------- //
        ee = m_ghost2id.cend();
        for (ii = m_ghost2id.cbegin(); ii != ee; ++ii) {
            rank_id = ii->first;
            i = m_ghost2id[rank_id].cbegin();
            e = m_ghost2id[rank_id].cend();
            while ( i != e ) {
/*DEBUG*/       out << "    processing cell: " << i->second << endl;
                flag_delete = true;
                ghost_idx = i->second;
                cell_ = &m_cells[ghost_idx];
                neighs = findCellNeighs(ghost_idx);
                n_neighs = neighs.size();
/*DEBUG*/       out << "    n_neighs = " << n_neighs << endl;
                for (j = 0; j < n_neighs; ++j) {
/*DEBUG*/           out << "      neigh: " << m_cells[neighs[j]].isInterior() << endl;
                    flag_delete &= ( !m_cells[neighs[j]].isInterior() );
                }
//TO BE REMOVED                 n_vertices = cell_->getVertexCount();
//TO BE REMOVED                 for ( j = 0; j < n_vertices; ++j ) {
//TO BE REMOVED                     ring_1 = Ring1(ghost_idx, j);
//TO BE REMOVED                     n_adj = ring_1.size();
//TO BE REMOVED                     for ( k = 0; k < n_adj; ++k) {
//TO BE REMOVED                         neigh_idx = ring_1[k];
//TO BE REMOVED                         flag_delete &= ( !cells[neigh_idx].isInterior() );
//TO BE REMOVED                     } //next k
//TO BE REMOVED                 } //next j
                if ( flag_delete )  {
/*DEBUG*/           out << "      deleting ghost: " << ghost_idx << endl;

                    // Ghost has to be deleted
                    deleteCell(ghost_idx);
                    i = m_ghost2id[rank_id].erase( i );
                }
                else ++i;
            } //next i
        } //next ii

/*DEBUG*/{        
/*DEBUG*/    out << "  ghost2idx = {";
/*DEBUG*/    for (i = m_ghost2id[rcv_rank].begin(); i != m_ghost2id[snd_rank].end(); ++i) {
/*DEBUG*/        out << " (" << i->first << ", " << i->second << ")";
/*DEBUG*/    }
/*DEBUG*/    out << " }" << endl << endl;
/*DEBUG*/    out << "  nCells: " << m_nInternals << endl;
/*DEBUG*/    out << "  nGhosts: " << m_nGhosts << endl;
/*DEBUG*/}
    }

    // Remove isolated vertices ============================================= //
/*DEBUG*/out << "* sender (rank#" << snd_rank << "), removing isolated vertices" << endl << endl;
    {
        // Scope variables -------------------------------------------------- //
        // none

        // Remove Isolated vertices ----------------------------------------- //
        deleteOrphanVertices();
    }

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
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving vertices" << endl;
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
            v_local_mapping[i] = it->get_id();
        } //next i
/*DEBUG*/out << "    received " << n_vertex << " vertices " << v_local_mapping << endl;
        
    }

    // Receive cells ======================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving cells" << endl;
    {
        // Scope variables -------------------------------------------------- //
        int                                             j;
        int                                             n_vertices;
        long                                            i, cell_idx, ghost_idx, neigh_idx;
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
        cell.setInterior( true );
        reserveCells( m_nInternals + m_nGhosts + n_cells );
        for ( i = 0; i < n_cells; ++i ) {

            // Stream cell data
            com_buff >> cell;
            cell_idx = cell.get_id();

/*DEBUG*/   {
/*DEBUG*/       out << "    received cell:" << endl;
/*DEBUG*/       cell.display(out, 6);
/*DEBUG*/   }

            // Update connectivity
            n_vertices = cell.getVertexCount();
            for (j = 0; j < n_vertices; ++j) {
                cell.setVertex( j, v_local_mapping[ cell.getVertex( j ) ] );
            } //next j

            // Insert cell in cells list & update ghost lists
            if ( m_ghost2id[snd_rank].find(cell_idx) != m_ghost2id[snd_rank].end() ) {
                ghost_idx = m_ghost2id[snd_rank][cell_idx];
                cell.set_id( ghost_idx );
                feedback << ghost_idx;
/*DEBUG*/       {
/*DEBUG*/           out << "    (already exists as ghost), remapped as:" << endl;
/*DEBUG*/           cell.display(out, 6);
/*DEBUG*/       }
                m_cells[ghost_idx] = move( cell );
                moveGhost2Internal(ghost_idx);
                m_ghost2id[snd_rank].erase( cell_idx );
                c_local_mapping[ cell_count ] = ghost_idx;
                ++cell_count;
            }
            else {

                // Add cell before last internal cell
/*DEBUG*/       out << "    (adding new cell), remapped as:" << endl;
                it = addCell( move(cell), generateCellId() );
                feedback << long( it->get_id() );
/*DEBUG*/       it->display(out, 6);
                c_local_mapping[ cell_count ] = it->get_id();
                ++cell_count;
            }
        } //next i

/*DEBUG*/{
/*DEBUG*/    out << "    received " << cell_count << " cells: ";
/*DEBUG*/    unordered_map<long, long>::const_iterator    ii;
/*DEBUG*/    for (ii = c_local_mapping.begin(); ii != c_local_mapping.end(); ++ii) {
/*DEBUG*/           out << "(" << ii->first << ", " << ii->second << "), ";
/*DEBUG*/    } //next i
/*DEBUG*/    out << endl;
/*DEBUG*/}

        // Communicate new IDs to sender ------------------------------------ //
/*DEBUG*/{
/*DEBUG*/    out << "    sending new IDs: ";
/*DEBUG*/    unordered_map<long, long>::const_iterator      m;
/*DEBUG*/    for (m = c_local_mapping.begin(); m != c_local_mapping.end(); ++m) {
/*DEBUG*/        out << " (" << m->first << ", " << m->second << ") ";
/*DEBUG*/    } //next m
/*DEBUG*/    out << endl << endl;
/*DEBUG*/}

        // Communicate buffer size
        MPI_Send( &feedback_size, 1, MPI_LONG, snd_rank, 6, m_communicator );

        // Communicate buffer
        MPI_Send( feedback.get_buffer(), feedback_size, MPI_CHAR, snd_rank, 7, m_communicator );

    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;

    // Receive ghosts ======================================================= //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), receiving ghosts" << endl;
    {

        // Scope variables -------------------------------------------------- //
        short                                   ghost_rank;
        int                                     j;
        int                                     n_vertices;
        long                                    i;
        long                                    cell_idx, ghost_idx;
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
            cell_idx = cell.get_id();
/*DEBUG*/   {
/*DEBUG*/       out << "    receiving ghost (existing on rank: " << ghost_rank << ")" << endl;
/*DEBUG*/       cell.display(out, 6);
/*DEBUG*/   }

            // Update connectivity
            n_vertices = cell.getVertexCount();
            for ( j = 0; j < n_vertices; ++j ) {
                cell.setVertex( j, v_local_mapping[ cell.getVertex( j ) ] );
            } //next j

            // Update cell list
            if ( m_ghost2id[ghost_rank].find(cell_idx) == m_ghost2id[ghost_rank].end() ) {

/*DEBUG*/       out << "    (non-existent ghost cell)" << endl;

                // Cell has to be added after the last internal cell
                cell.setInterior( false );
                it = addCell( move(cell), generateCellId() );

                c_local_mapping[ n_cells + ghost_count ] = it->get_id();
                m_ghost2id[ghost_rank][ cell_idx ] = it->get_id();
                ++ghost_count;

/*DEBUG*/       {
/*DEBUG*/           out << "    remapped as:" << endl;
/*DEBUG*/           it->display(out, 6);
/*DEBUG*/       }
            }
            else {
/*DEBUG*/       out << "    (cell already exists as ghost)" << endl;
                cell.set_id(m_ghost2id[ghost_rank][cell_idx]);
                cell.setInterior( false );
/*DEBUG*/       {
/*DEBUG*/           out << "    remapped as:" << endl;
/*DEBUG*/           cell.display(out, 6);
/*DEBUG*/       }

                // Cell has to be overwritten
                m_cells[m_ghost2id[ghost_rank][cell_idx]] = move(cell);

                c_local_mapping[ n_cells + ghost_count ] = m_ghost2id[ghost_rank][cell_idx];
                ++ghost_count;
            }

        } //next i

/*DEBUG*/{
/*DEBUG*/    unordered_map<long, long>::const_iterator      m;
/*DEBUG*/    out << "    c_local_mapping is now: {";
/*DEBUG*/    for ( m = c_local_mapping.begin(); m != c_local_mapping.end(); ++m) {
/*DEBUG*/        out << " (" << m->first << ", " << m->second << ")";
/*DEBUG*/    } //next m
/*DEBUG*/    out << " }" << endl << endl;
/*DEBUG*/}
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;

    // Update adjacencies =================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), updating adjacencies" << endl;
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
/*DEBUG*/   out << "    updating adjacencies for cell ID " << cell_idx << endl;
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
/*DEBUG*/   out << "    cell ID " << cell_idx << " is now" << endl;
/*DEBUG*/   cell_->display(out, 6);
        } //next m

/*DEBUG*/out << "  ghost2idx = {";
/*DEBUG*/for (m = m_ghost2id[snd_rank].begin(); m != m_ghost2id[snd_rank].end(); ++m) {
/*DEBUG*/    out << " (" << m->first << ", " << m->second << ")";
/*DEBUG*/}
/*DEBUG*/out << " }" << endl << endl;

    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;

    // Remove duplicated vertices =========================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), removing duplicated vertices" << endl;
    {
        // Scope variables -------------------------------------------------- //
        // none

        // Remove duplicated vertices --------------------------------------- //
        out << "    n. vertices is: " << m_nVertices << endl;
        deleteCoincidentVertex();
        deleteOrphanVertices();
        out << "    (after cleaning), n.vertices is: " << m_nVertices << endl << endl;
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;

    // Update adjacencies =================================================== //
/*DEBUG*/out << "* receiver (rank#" << rcv_rank << "), updating adjacencies" << endl;
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
/*DEBUG*/out << "    updating adjacencies of simplicies: " << cell_list << endl << endl;
        updateAdjacencies(cell_list);
/*DEBUG*/{
/*DEBUG*/    for (int ii = 0; ii < cell_list.size(); ++ii) {
/*DEBUG*/           out << "    cell ID: " << cell_list[ii] << " is now: " << endl;
/*DEBUG*/           m_cells[cell_list[ii]].display(out, 4);
/*DEBUG*/    } //next ii
/*DEBUG*/}
    }
/*DEBUG*/out << "  nInternals: " << m_nInternals << endl;
/*DEBUG*/out << "  nGhosts: " << m_nGhosts << endl;

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
}

/*DEBUG*/{
/*DEBUG*/    out << "* display mesh infos" << endl;
/*DEBUG*/    displayCells(out);
/*DEBUG*/    out.close();
/*DEBUG*/}

return; }

/*!
	@}
*/

}

#endif
