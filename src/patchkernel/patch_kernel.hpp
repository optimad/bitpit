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

#ifndef __BITPIT_PATCH_KERNEL_HPP__
#define __BITPIT_PATCH_KERNEL_HPP__

#include <cstddef>
#include <deque>
#include <iostream>
#include <memory>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif
#include <set>
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_IO.hpp"
#if BITPIT_ENABLE_MPI==1
#	include "bitpit_communications.hpp"
#endif

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class PatchKernel : public VTKBaseStreamer {

friend class PatchInfo;
#if BITPIT_ENABLE_MPI==1
friend class PatchGlobalInfo;
#endif
friend class PatchManager;

public:
	typedef PiercedVector<Vertex>::iterator VertexIterator;
	typedef PiercedVector<Cell>::iterator CellIterator;
	typedef PiercedVector<Interface>::iterator InterfaceIterator;

	typedef PiercedVector<Vertex>::const_iterator VertexConstIterator;
	typedef PiercedVector<Cell>::const_iterator CellConstIterator;
	typedef PiercedVector<Interface>::const_iterator InterfaceConstIterator;

	typedef PiercedVector<Vertex>::range VertexRange;
	typedef PiercedVector<Cell>::range CellRange;
	typedef PiercedVector<Interface>::range InterfaceRange;

	typedef PiercedVector<Vertex>::const_range VertexConstRange;
	typedef PiercedVector<Cell>::const_range CellConstRange;
	typedef PiercedVector<Interface>::const_range InterfaceConstRange;

	enum WriteTarget {
		WRITE_TARGET_CELLS_ALL
#if BITPIT_ENABLE_MPI
		, WRITE_TARGET_CELLS_INTERNAL
#endif
	};

	/*!
		Functional for compare the position of two cells
	*/
	struct CellPositionLess
	{
		CellPositionLess(PatchKernel &patch, bool native = true)
			: m_patch(patch), m_native(native)
		{
		}

		bool operator()(const long &id_1, const long &id_2) const
		{
			std::array<double, 3> centroid_1;
			std::array<double, 3> centroid_2;
			if (m_native) {
				centroid_1 = m_patch.evalCellCentroid(id_1);
				centroid_2 = m_patch.evalCellCentroid(id_2);
			} else {
				centroid_1 = m_patch.PatchKernel::evalCellCentroid(id_1);
				centroid_2 = m_patch.PatchKernel::evalCellCentroid(id_2);
			}

			for (int k = 0; k < 3; ++k) {
				if (std::abs(centroid_1[k] - centroid_2[k]) <= m_patch.getTol()) {
					continue;
				}

				return centroid_1[k] < centroid_2[k];
			}

			// If we are here the two cell centroids coincide. It's not
			// possible to define an order for the two cells.
			std::ostringstream stream;
			stream << "It was not possible to define an order for cells " << id_1 << " and " << id_2 << ". ";
			stream << "The two cells have the same centroid.";
			throw std::runtime_error (stream.str());
		}

		PatchKernel &m_patch;
		bool m_native;
	};

	/*!
		Functional for compare the position of two cells
	*/
	struct CellPositionGreater
	{
		CellPositionGreater(PatchKernel &patch, bool native = true)
			: m_patch(patch), m_native(native)
		{
		}

		bool operator()(const long &id_1, const long &id_2) const
		{
			std::array<double, 3> centroid_1;
			std::array<double, 3> centroid_2;
			if (m_native) {
				centroid_1 = m_patch.evalCellCentroid(id_1);
				centroid_2 = m_patch.evalCellCentroid(id_2);
			} else {
				centroid_1 = m_patch.PatchKernel::evalCellCentroid(id_1);
				centroid_2 = m_patch.PatchKernel::evalCellCentroid(id_2);
			}

			for (int k = 0; k < 3; ++k) {
				if (std::abs(centroid_1[k] - centroid_2[k]) <= m_patch.getTol()) {
					continue;
				}

				return centroid_1[k] > centroid_2[k];
			}

			// If we are here the two cell centroids coincide. It's not
			// possible to define an order for the two cells.
			std::ostringstream stream;
			stream << "It was not possible to define an order for cells " << id_1 << " and " << id_2 << ". ";
			stream << "The two cells have the same centroid.";
			throw std::runtime_error (stream.str());
		}

		PatchKernel &m_patch;
		bool m_native;
	};

	/*!
		Spawn status
	*/
	enum SpawnStatus {
		SPAWN_UNNEEDED = -1,
		SPAWN_NEEDED,
		SPAWN_DONE
	};

	/*!
		Adaption status
	*/
	enum AdaptionStatus {
		ADAPTION_UNSUPPORTED = -1,
		ADAPTION_CLEAN,
		ADAPTION_DIRTY,
		ADAPTION_PREPARED,
		ADAPTION_ALTERED
	};

	virtual ~PatchKernel();

	virtual void reset();
	void resetVertices();
	void resetCells();
	void resetInterfaces();

	bool reserveVertices(size_t nVertices);
	bool reserveCells(size_t nCells);
	bool reserveInterfaces(size_t nInterfaces);

	std::vector<adaption::Info> update(bool trackAdaption = true, bool squeezeStorage = false);

	SpawnStatus getSpawnStatus() const;
	std::vector<adaption::Info> spawn(bool trackSpawn);

	AdaptionStatus getAdaptionStatus(bool global = false) const;
	std::vector<adaption::Info> adaption(bool trackAdaption = true, bool squeezeStorage = false);
	std::vector<adaption::Info> adaptionPrepare(bool trackAdaption = true);
	std::vector<adaption::Info> adaptionAlter(bool trackAdaption = true, bool squeezeStorage = false);
	void adaptionCleanup();

	void markCellForRefinement(const long &id);
	void markCellForCoarsening(const long &id);
	void enableCellBalancing(const long &id, bool enabled);

	bool isDirty(bool global = false) const;
	bool isExpert() const;

	int getId() const;
	int getDimension() const;
	virtual void setDimension(int dimension);
	bool isThreeDimensional() const;

	virtual long getVertexCount() const;
	PiercedVector<Vertex> &getVertices();
	const PiercedVector<Vertex> &getVertices() const;
	Vertex &getVertex(const long &id);
	const Vertex & getVertex(const long &id) const;
	const std::array<double, 3> & getVertexCoords(const long &id) const;
	long generateVertexId();
	VertexIterator addVertex(const std::array<double, 3> &coords, const long &id = Vertex::NULL_ID);
	VertexIterator addVertex(const Vertex &source, long id = Vertex::NULL_ID);
	VertexIterator addVertex(Vertex &&source, long id = Vertex::NULL_ID);
	long countFreeVertices() const;
	long countOrphanVertices() const;
	std::vector<long> findOrphanVertices();
	bool deleteOrphanVertices();
	std::vector<long> collapseCoincidentVertices(int nBins = 128);
	bool deleteCoincidentVertices(int nBins = 128);

	VertexIterator getVertexIterator(const long &id);
	VertexIterator vertexBegin();
	VertexIterator vertexEnd();

	VertexConstIterator getVertexConstIterator(const long &id) const;
	VertexConstIterator vertexConstBegin() const;
	VertexConstIterator vertexConstEnd() const;

	virtual long getCellCount() const;
	long getInternalCount() const;
	long getGhostCount() const;
	PiercedVector<Cell> &getCells();
	const PiercedVector<Cell> &getCells() const;
	Cell &getCell(const long &id);
	const Cell &getCell(const long &id) const;
	virtual ElementInfo::Type getCellType(const long &id) const;
	Cell &getLastInternal();
	const Cell &getLastInternal() const;
	Cell &getFirstGhost();
	const Cell &getFirstGhost() const;
	long generateCellId();
	CellIterator addCell(ElementInfo::Type type, const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, std::unique_ptr<long[]> &&connect, const long &id = Element::NULL_ID);
	CellIterator addCell(ElementInfo::Type type, bool interior, const std::vector<long> &connect, const long &id = Element::NULL_ID);
	CellIterator addCell(const Cell &source, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, long id = Element::NULL_ID);
	bool deleteCell(const long &id, bool updateNeighs = true, bool delayed = false);
	bool deleteCells(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
	bool setCellInternal(const long &id, bool isInternal);
	CellIterator moveGhost2Internal(const long &id);
	CellIterator moveInternal2Ghost(const long &id);
	virtual double evalCellSize(const long &id) const = 0;
	long countFreeCells() const;
	long countOrphanCells() const;
	virtual std::array<double, 3> evalCellCentroid(const long &id) const;
	std::vector<long> findCellNeighs(const long &id) const;
	std::vector<long> findCellNeighs(const long &id, int codimension, bool complete = true) const;
	std::vector<long> findCellFaceNeighs(const long &id) const;
	std::vector<long> findCellFaceNeighs(const long &id, const int &face) const;
	std::vector<long> findCellEdgeNeighs(const long &id, bool complete = true) const;
	std::vector<long> findCellEdgeNeighs(const long &id, const int &edge) const;
	std::vector<long> findCellVertexNeighs(const long &id, bool complete = true) const;
	std::vector<long> findCellVertexNeighs(const long &id, const int &vertex) const;
	std::vector<long> findCellVertexOneRing(const long &id, const int &vertex) const;
    void findFaceNeighCell(const long &cell_idx, const long &neigh_idx, int &face_loc_idx, int &intf_loc_idx);
	std::set<int> getInternalPIDs();
	std::vector<long> getInternalsByPID(int pid);

	CellIterator getCellIterator(const long &id);
	CellIterator cellBegin();
	CellIterator cellEnd();
	CellIterator internalBegin();
	CellIterator internalEnd();
	CellIterator ghostBegin();
	CellIterator ghostEnd();

	CellConstIterator getCellConstIterator(const long &id) const;
	CellConstIterator cellConstBegin() const;
	CellConstIterator cellConstEnd() const;
	CellConstIterator internalConstBegin() const;
	CellConstIterator internalConstEnd() const;
	CellConstIterator ghostConstBegin() const;
	CellConstIterator ghostConstEnd() const;

	virtual long getInterfaceCount() const;
	PiercedVector<Interface> &getInterfaces();
	const PiercedVector<Interface> & getInterfaces() const;
	Interface &getInterface(const long &id);
	const Interface &getInterface(const long &id) const;
	virtual ElementInfo::Type getInterfaceType(const long &id) const;
	long generateInterfaceId();
	InterfaceIterator addInterface(ElementInfo::Type type, const long &id = Element::NULL_ID);
	InterfaceIterator addInterface(const Interface &source, long id = Element::NULL_ID);
	InterfaceIterator addInterface(Interface &&source, long id = Element::NULL_ID);
	bool deleteInterface(const long &id, bool updateNeighs = true, bool delayed = false);
	bool deleteInterfaces(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
	long countFreeInterfaces() const;
	long countOrphanInterfaces() const;
	virtual std::array<double, 3> evalInterfaceCentroid(const long &id) const;

	InterfaceIterator getInterfaceIterator(const long &id);
	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

	InterfaceConstIterator getInterfaceConstIterator(const long &id) const;
	InterfaceConstIterator interfaceConstBegin() const;
	InterfaceConstIterator interfaceConstEnd() const;

	long countFaces() const;
	long countFreeFaces() const;

	bool sort();
	bool sortVertices();
	bool sortCells();
	bool sortInterfaces();

	bool squeeze();
	bool squeezeVertices();
	bool squeezeCells();
	bool squeezeInterfaces();

	long locatePoint(const double &x, const double &y, const double &z);
	virtual long locatePoint(const std::array<double, 3> &point) = 0;
        bool isSameFace(const long &, const int&, const long&, const int&);

	virtual void buildAdjacencies(bool resetAdjacencies = true);
	virtual void updateAdjacencies(const std::vector<long> &cellIds, bool resetAdjacencies = true);

	virtual void buildInterfaces(bool resetInterfaces = true);
	virtual void updateInterfaces(const std::vector<long> &cellIds, bool resetInterfaces = true);

	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);
	bool isBoundingBoxDirty(bool global = false) const;
	void updateBoundingBox(bool forcedUpdated = false);

	std::unordered_map<long, long> binSortVertex(PiercedVector<Vertex> vertices, int nBins = 128);
    std::unordered_map<long, long> binSortVertex(int nBins = 128);

	virtual void translate(std::array<double, 3> translation);
	void translate(double sx, double sy, double sz);
	virtual void scale(std::array<double, 3> scaling);
	void scale(double scaling);
	void scale(double sx, double sy, double sz);

	void setTol(double tolerance);
	double getTol() const;
	void resetTol();
	bool isTolCustomized() const;

	void extractEnvelope(PatchKernel &envelope) const;

	void displayTopologyStats(std::ostream &out, unsigned int padding = 0) const;
	void displayVertices(std::ostream &out, unsigned int padding = 0) const;
	void displayCells(std::ostream &out, unsigned int padding = 0) const;
	void displayInterfaces(std::ostream &out, unsigned int padding = 0) const;

	VTKUnstructuredGrid & getVTK();
	WriteTarget getVTKWriteTarget() const;
	void setVTKWriteTarget(WriteTarget targetCells);
	const PatchKernel::CellConstRange & getVTKCellWriteRange() const;
	void write(VTKWriteMode mode = VTKWriteMode::DEFAULT);
	void write(std::string name, VTKWriteMode mode = VTKWriteMode::DEFAULT);

	void flushData(std::fstream &stream, std::string name, VTKFormat format );

	int getDumpVersion() const;
	void dump(std::ostream &stream);
	void restore(std::istream &stream, bool reregister = false);

	void consecutiveRenumberVertices(long offset = 0);
	void consecutiveRenumberCells(long offset = 0);
	void consecutiveRenumberInterfaces(long offset = 0);
	void consecutiveRenumber(long offsetVertices, long offsetCells, long offsetInterfaces);

#if BITPIT_ENABLE_MPI==1
	virtual void setCommunicator(MPI_Comm communicator);
	void freeCommunicator();
	bool isCommunicatorSet() const;
	const MPI_Comm & getCommunicator() const;
	int getRank() const;
	int getProcessorCount() const;

	int getCellRank(const long &id) const;

	bool isRankNeighbour(int rank);
	std::vector<int> getNeighbourRanks();
	std::unordered_map<int, std::vector<long>> & getGhostExchangeTargets();
	const std::unordered_map<int, std::vector<long>> & getGhostExchangeTargets() const;
	std::vector<long> & getGhostExchangeTargets(int rank);
	const std::vector<long> & getGhostExchangeTargets(int rank) const;

	std::unordered_map<int, std::vector<long>> & getGhostExchangeSources();
	const std::unordered_map<int, std::vector<long>> & getGhostExchangeSources() const;
	std::vector<long> & getGhostExchangeSources(int rank);
	const std::vector<long> & getGhostExchangeSources(int rank) const;

	std::vector<adaption::Info> partition(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackChanges, bool squeezeStorage = false);
	std::vector<adaption::Info> partition(const std::vector<int> &cellRanks, bool trackChanges, bool squeezeStorage = false);
	std::vector<adaption::Info> partition(MPI_Comm communicator, bool trackChanges, bool squeezeStorage = false);
	std::vector<adaption::Info> partition(bool trackChanges, bool squeezeStorage = false);
	std::vector<adaption::Info> balancePartition(bool trackChanges, bool squeezeStorage = false);
	bool isPartitioned() const;

	adaption::Info sendCells(const int &sendRank, const int &recvRank, const std::vector<long> &cellsToSend);
#endif

protected:
	PiercedVector<Vertex> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	IndexGenerator m_vertexIdGenerator;
	IndexGenerator m_interfaceIdGenerator;
	IndexGenerator m_cellIdGenerator;

	long m_nInternals;
	long m_nGhosts;

	long m_lastInternalId;
	long m_firstGhostId;

	VTKUnstructuredGrid m_vtk ;
	WriteTarget m_vtkWriteTarget;
	CellConstRange m_vtkCellRange;
	std::unordered_map<long, long> m_vtkVertexMap;

	PatchKernel(bool expert);
	PatchKernel(const int &dimension, bool expert);
	PatchKernel(const int &id, const int &dimension, bool expert);

	void clearBoundingBox();
	bool isBoundingBoxFrozen() const;
	void setBoundingBoxFrozen(bool frozen);
	void setBoundingBoxDirty(bool dirty);
	void setBoundingBox(const std::array<double, 3> &minPoint, const std::array<double, 3> &maxPoint);

	bool deleteVertex(const long &id, bool delayed = false);
	bool deleteVertices(const std::vector<long> &ids, bool delayed = false);

	void dumpInterfaces(std::ostream &stream);
	void restoreInterfaces(std::istream &stream);

	void setSpawnStatus(SpawnStatus status);
	virtual std::vector<adaption::Info> _spawn(bool trackAdaption);

	void setAdaptionStatus(AdaptionStatus status);
	virtual std::vector<adaption::Info> _adaptionPrepare(bool trackAdaption);
	virtual std::vector<adaption::Info> _adaptionAlter(bool trackAdaption);
	virtual void _adaptionCleanup();
	virtual bool _markCellForRefinement(const long &id);
	virtual bool _markCellForCoarsening(const long &id);
	virtual bool _enableCellBalancing(const long &id, bool enabled);

	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	virtual int _getDumpVersion() const = 0;
	virtual void _dump(std::ostream &stream) = 0;
	virtual void _restore(std::istream &stream) = 0;

	virtual long _getCellNativeIndex(long id) const;

	virtual std::vector<long> _findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList = std::vector<long>()) const;
	virtual std::vector<long> _findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;
	virtual std::vector<long> _findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;

	void setExpert(bool expert);

	void addPointToBoundingBox(const std::array<double, 3> &point);
	void removePointFromBoundingBox(const std::array<double, 3> &point, bool delayedBoxUpdate = false);
#if BITPIT_ENABLE_MPI==1
	virtual std::vector<adaption::Info> _balancePartition(bool trackChanges);

	void setPartitioned(bool partitioned);

	void setGhostOwner(int id, int rank, bool updateExchangeData = false);
	void unsetGhostOwner(int id, bool updateExchangeData = false);
	void clearGhostOwners(bool updateExchangeData = false);

	void deleteGhostExchangeData();
	void deleteGhostExchangeData(int rank);

	void buildGhostExchangeData();
	void buildGhostExchangeData(int rank);
	void buildGhostExchangeData(const std::vector<int> &rank);

	void addGhostsToExchangeTargets(const std::vector<long> &ghostIds);
	void addGhostToExchangeTargets(const long ghostId);

	void removeGhostsFromExchangeTargets(const std::vector<long> &ghostIds);
	void removeGhostFromExchangeTargets(const long ghostId);
#endif

	template<typename item_t, typename id_t = long>
	std::unordered_map<id_t, id_t> consecutiveItemRenumbering(PiercedVector<item_t, id_t> &container, long offset);

	template<typename item_t, typename id_t = long>
	void mappedItemRenumbering(PiercedVector<item_t, id_t> &container, const std::unordered_map<id_t, id_t> &renumberMap);

private:
	double DEFAULT_TOLERANCE = 1e-14;

	bool m_boxFrozen;
	bool m_boxDirty;
	std::array<double, 3> m_boxMinPoint;
	std::array<double, 3> m_boxMaxPoint;
	std::array<int, 3> m_boxMinCounter;
	std::array<int, 3> m_boxMaxCounter;

	SpawnStatus m_spawnStatus;

	AdaptionStatus m_adaptionStatus;

	bool m_expert;

	int m_id;
	int m_dimension;

	bool m_hasCustomTolerance;
	double m_tolerance;

	int m_rank;
	int m_nProcessors;
#if BITPIT_ENABLE_MPI==1
	MPI_Comm m_communicator;
	bool m_partitioned;

	std::unordered_map<long, int> m_ghostOwners;
	std::unordered_map<int, std::vector<long>> m_ghostExchangeTargets;
	std::unordered_map<int, std::vector<long>> m_ghostExchangeSources;

	void addExchangeSources(const std::vector<long> &ghostIds);

    adaption::Info sendCells_sender(const int &recvRank, const std::vector<long> &cellsToSend);
    adaption::Info sendCells_receiver(const int &sendRank);
    adaption::Info sendCells_notified(const int &sendRank, const int &recvRank);
#endif

	void initialize();

	void beginAlteration();
	void endAlteration(bool squeezeStorage = false);

	void buildCellInterface(Cell *cell_1, int face_1, Cell *cell_2, int face_2, long interfaceId = Element::NULL_ID);

	VertexIterator createVertex(const std::array<double, 3> &coords, long id = Vertex::NULL_ID);
	InterfaceIterator createInterface(ElementInfo::Type type, long id = Element::NULL_ID);
	CellIterator createCell(ElementInfo::Type type, bool interior, long id = Element::NULL_ID);

	int findAdjoinNeighFace(const long &cellId, const long &neighId) const;

	void setId(int id);

	std::array<double, 3> evalElementCentroid(const Element &element) const;

	void mergeAdaptionInfo(std::vector<adaption::Info> &&source, std::vector<adaption::Info> &destination);
};

}

// Template implementation
#include "patch_kernel.tpp"

#endif
