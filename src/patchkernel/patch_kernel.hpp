/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "bitpit_containers.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "compiler.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class PatchKernel : public VTKBaseStreamer {

friend class PatchInfo;
friend class PatchNumberingInfo;
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
		Functional for comparing the position of two points.

		To identify the order of the two points, their coordinates will be
		compared. Before comparing the coordinates, it is necessary check
		if the two points are coincident within the specified tolerance.
	*/
	struct PointPositionLess
	{
		PointPositionLess(double tolerance)
			: m_tolerance(tolerance)
		{
		}

		virtual ~PointPositionLess() = default;

		bool operator()(std::array<double, 3> point_1, std::array<double, 3> point_2) const
		{
			// Check if the points are coincident
			//
			// To check if two points are coincident, the specified tolerance
			// will be used.
			bool coincident = true;
			for (int k = 0; k < 3; ++k) {
				if (!utils::DoubleFloatingEqual()(point_1[k], point_2[k], m_tolerance)) {
					coincident = false;
					break;
				}
			}

			if (coincident) {
				return false;
			}

			// Compare the position of the points
			//
			// If the points are not coincident, we loop over the coordinates
			// and we compare the first non-coincident coordinate.
			//
			// In order to guarantee a consistent ordering, the test to check
			// if two coordinates are coincident should be performed with no
			// tolerance.
			for (int k = 0; k < 3; ++k) {
				if (point_1[k] == point_2[k]) {
					continue;
				}

				return (point_1[k] < point_2[k]);
			}

			return false;
		}

		const double m_tolerance;
	};

	/*!
		Functional for comparing the position of two vertices.

		The comparison is made with respect to the vertex coordinates.
	*/
	struct VertexPositionLess : private PointPositionLess
	{
		VertexPositionLess(const PatchKernel &patch)
			: PointPositionLess(patch.getTol()), m_patch(patch)
		{
		}

		virtual ~VertexPositionLess() = default;

		bool operator()(long id_1, long id_2) const
		{
			if (id_1 == id_2) {
				return false;
			}

			const std::array<double, 3> &coords_1 = m_patch.getVertexCoords(id_1);
			const std::array<double, 3> &coords_2 = m_patch.getVertexCoords(id_2);

			return PointPositionLess::operator()(coords_1, coords_2);
		}

		const PatchKernel &m_patch;
	};

	/*!
		Functional for comparing the position of two vertices.

		The comparison is made with respect to the vertex coordinates.
	*/
	struct VertexPositionGreater : private VertexPositionLess
	{
		VertexPositionGreater(const PatchKernel &patch)
			: VertexPositionLess(patch)
		{
		}

		bool operator()(long id_1, long id_2) const
		{
			return !VertexPositionLess::operator()(id_1, id_2);
		}
	};

	/*!
		Functional for comparing the position of two cells.

		The comparison is made with respect to the cell centroid.
	*/
	struct CellPositionLess : private PointPositionLess
	{
		CellPositionLess(const PatchKernel &patch, bool native = true)
			: PointPositionLess(patch.getTol()), m_patch(patch), m_native(native)
		{
		}

		virtual ~CellPositionLess() = default;

		bool operator()(long id_1, long id_2) const
		{
			if (id_1 == id_2) {
				return false;
			}

			std::array<double, 3> centroid_1;
			std::array<double, 3> centroid_2;
			if (m_native) {
				centroid_1 = m_patch.evalCellCentroid(id_1);
				centroid_2 = m_patch.evalCellCentroid(id_2);
			} else {
				centroid_1 = m_patch.PatchKernel::evalCellCentroid(id_1);
				centroid_2 = m_patch.PatchKernel::evalCellCentroid(id_2);
			}

			return PointPositionLess::operator()(centroid_1, centroid_2);
		}

		const PatchKernel &m_patch;
		bool m_native;
	};

	/*!
		Functional for comparing the position of two cells.

		The comparison is made with respect to the cell centroid.
	*/
	struct CellPositionGreater : private CellPositionLess
	{
		CellPositionGreater(const PatchKernel &patch, bool native = true)
			: CellPositionLess(patch, native)
		{
		}

		bool operator()(long id_1, long id_2) const
		{
			return !CellPositionLess::operator()(id_1, id_2);
		}
	};

	/*!
		Functional for comparing the position of two cells.

		WARNING: the function is faster than the comparison based on cell
		centroid but its result is not indepenedent of the vertex order.

		The comparison is made with respect to the position of cell vertices.
	*/
	struct CellFuzzyPositionLess : private PointPositionLess
	{
		CellFuzzyPositionLess(PatchKernel &patch)
			: PointPositionLess(patch.getTol()), m_patch(patch)
		{
		}

		virtual ~CellFuzzyPositionLess() = default;

		bool operator()(long id_1, long id_2) const
		{
			if (id_1 == id_2) {
				return false;
			}

			// Select the first vertex of the first cell
			long vertexId_1 = m_patch.getCell(id_1).getVertexId(0);

			// The vertex of the second cell is choosen as the first vertex on
			// that cell not equal to the selected vertex of the first cell.
			long vertexId_2 = Vertex::NULL_ID;
			for (long candidateVertexId_2 : m_patch.getCell(id_2).getVertexIds()) {
				if (vertexId_1 != candidateVertexId_2) {
					vertexId_2 = candidateVertexId_2;
					break;
				}
			}

			if (vertexId_2 == Vertex::NULL_ID) {
				return false;
			}

			// Compare the two vertices
			const std::array<double, 3> &vertexCoords_1 = m_patch.getVertex(vertexId_1).getCoords();
			const std::array<double, 3> &vertexCoords_2 = m_patch.getVertex(vertexId_2).getCoords();

			return PointPositionLess::operator()(vertexCoords_1, vertexCoords_2);
		}

		PatchKernel &m_patch;
	};

	/*!
		Functional for comparing the position of two cells.

		WARNING: the function is faster than the comparison based on cell
		centroid but its result is not indepenedent of the vertex order.

		The comparison is made with respect to the position of cell vertices.
	*/
	struct CellFuzzyPositionGreater : private CellFuzzyPositionLess
	{
		CellFuzzyPositionGreater(PatchKernel &patch)
			: CellFuzzyPositionLess(patch)
		{
		}

		bool operator()(long id_1, long id_2) const
		{
			return !CellFuzzyPositionLess::operator()(id_1, id_2);
		}
	};

	/*!
		Adjacencies build strategy
	*/
	enum AdjacenciesBuildStrategy {
		ADJACENCIES_NONE = -1,
		ADJACENCIES_AUTOMATIC
	};

	/*!
		Interfaces build strategy
	*/
	enum InterfacesBuildStrategy {
		INTERFACES_NONE = -1,
		INTERFACES_AUTOMATIC
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
		ADAPTION_UNSUPPORTED = -1, //!< Adaptation is not supported by the patch
		ADAPTION_CLEAN,            //!< There are no adaptation operations in progress
		ADAPTION_DIRTY,            //!< There are cells marked for adaptation, but the changes
		                           //!< has not been performed
		ADAPTION_PREPARED,         //!< Partitioning prepare has been called
		ADAPTION_ALTERED           //!< Partitioning alters has been called
	};

	/*!
		Partitioning status
	*/
	enum PartitioningStatus {
		PARTITIONING_UNSUPPORTED = -1, //!< Partitioning is not supported by the patch
		PARTITIONING_CLEAN,            //!< There are no partitioning operations in progress
		PARTITIONING_PREPARED,         //!< Partitioning prepare has been called
		PARTITIONING_ALTERED           //!< Partitioning alter has been called
	};

	PatchKernel(PatchKernel &&other);
	PatchKernel & operator=(PatchKernel &&other);

	virtual ~PatchKernel();

	template<typename patch_t>
	static std::unique_ptr<patch_t> clone(const patch_t *original);

	virtual std::unique_ptr<PatchKernel> clone() const = 0;

	void setId(int id);

	virtual void reset();
	virtual void resetVertices();
	virtual void resetCells();
	virtual std::vector<adaption::Info> resetInterfaces(bool trackAdaption = false);

	bool reserveVertices(size_t nVertices);
	bool reserveCells(size_t nCells);
	bool reserveInterfaces(size_t nInterfaces);

	std::vector<adaption::Info> update(bool trackAdaption = true, bool squeezeStorage = false);

	virtual void simulateCellUpdate(const long id, adaption::Marker marker, std::vector<Cell> *virtualCells, PiercedVector<Vertex, long> *virtualVertices) const;

	SpawnStatus getSpawnStatus() const;
	std::vector<adaption::Info> spawn(bool trackSpawn);

	bool isAdaptionSupported() const;
	AdaptionStatus getAdaptionStatus(bool global = false) const;
	std::vector<adaption::Info> adaption(bool trackAdaption = true, bool squeezeStorage = false);
	std::vector<adaption::Info> adaptionPrepare(bool trackAdaption = true);
	std::vector<adaption::Info> adaptionAlter(bool trackAdaption = true, bool squeezeStorage = false);
	void adaptionCleanup();

	virtual void settleAdaptionMarkers();

	void markCellForRefinement(long id);
	void markCellForCoarsening(long id);
	void resetCellAdaptionMarker(long id);
    adaption::Marker getCellAdaptionMarker(long id);
	void enableCellBalancing(long id, bool enabled);

	bool isDirty(bool global = false) const;
	bool isExpert() const;

	int getId() const;
	int getDimension() const;
	virtual void setDimension(int dimension);
	bool isThreeDimensional() const;

	virtual int getVolumeCodimension() const = 0;
	virtual int getSurfaceCodimension() const = 0;
	virtual int getLineCodimension() const = 0;
	virtual int getPointCodimension() const = 0;

	bool empty(bool global = true) const;

	bool isVertexAutoIndexingEnabled() const;
	void setVertexAutoIndexing(bool enabled);

	virtual long getVertexCount() const;
	long getInternalVertexCount() const;
#if BITPIT_ENABLE_MPI==1
	long getGhostVertexCount() const;
#endif
	PiercedVector<Vertex> &getVertices();
	const PiercedVector<Vertex> &getVertices() const;
	Vertex &getVertex(long id);
	const Vertex & getVertex(long id) const;
	Vertex &getLastInternalVertex();
	const Vertex &getLastInternalVertex() const;
#if BITPIT_ENABLE_MPI==1
	Vertex &getFirstGhostVertex();
	const Vertex &getFirstGhostVertex() const;
#endif
	const std::array<double, 3> & getVertexCoords(long id) const;
	void getVertexCoords(std::size_t nVertices, const long *ids, std::unique_ptr<std::array<double, 3>[]> *coordinates) const;
	void getVertexCoords(std::size_t nVertices, const long *ids, std::array<double, 3> *coordinates) const;
	VertexIterator addVertex(const Vertex &source, long id = Vertex::NULL_ID);
	VertexIterator addVertex(Vertex &&source, long id = Vertex::NULL_ID);
	VertexIterator addVertex(const std::array<double, 3> &coords, long id = Vertex::NULL_ID);
	BITPIT_DEPRECATED(long countFreeVertices() const);
	long countBorderVertices() const;
	long countOrphanVertices() const;
	std::vector<long> findOrphanVertices();
	bool deleteOrphanVertices();
	std::vector<long> collapseCoincidentVertices();
	bool deleteCoincidentVertices();

	VertexIterator getVertexIterator(long id);
	VertexIterator vertexBegin();
	VertexIterator vertexEnd();
	VertexIterator internalVertexBegin();
	VertexIterator internalVertexEnd();
#if BITPIT_ENABLE_MPI==1
	VertexIterator ghostVertexBegin();
	VertexIterator ghostVertexEnd();
#endif

	VertexConstIterator getVertexConstIterator(long id) const;
	VertexConstIterator vertexConstBegin() const;
	VertexConstIterator vertexConstEnd() const;
	VertexConstIterator internalVertexConstBegin() const;
	VertexConstIterator internalVertexConstEnd() const;
#if BITPIT_ENABLE_MPI==1
	VertexConstIterator ghostVertexConstBegin() const;
	VertexConstIterator ghostVertexConstEnd() const;
#endif

	bool isCellAutoIndexingEnabled() const;
	void setCellAutoIndexing(bool enabled);

	virtual long getCellCount() const;
	long getInternalCellCount() const;
	BITPIT_DEPRECATED(long getInternalCount() const);
#if BITPIT_ENABLE_MPI==1
	long getGhostCellCount() const;
    BITPIT_DEPRECATED(long getGhostCount() const);
#endif
	PiercedVector<Cell> &getCells();
	const PiercedVector<Cell> &getCells() const;
	Cell &getCell(long id);
	const Cell &getCell(long id) const;
	virtual ElementType getCellType(long id) const;
	Cell &getLastInternalCell();
	const Cell &getLastInternalCell() const;
#if BITPIT_ENABLE_MPI==1
	Cell & getFirstGhostCell();
	BITPIT_DEPRECATED(Cell & getFirstGhost());
	const Cell & getFirstGhostCell() const;
	BITPIT_DEPRECATED(const Cell & getFirstGhost() const);
#endif
	CellIterator addCell(const Cell &source, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, const std::vector<long> &connectivity, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id = Element::NULL_ID);
#if BITPIT_ENABLE_MPI==1
	CellIterator addCell(const Cell &source, int owner, long id = Element::NULL_ID);
	CellIterator addCell(const Cell &source, int owner, int haloLayer, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, int owner, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, int owner, int haloLayer, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, int owner, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, int owner, int haloLayer, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, const std::vector<long> &connectivity, int owner, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, const std::vector<long> &connectivity, int owner, int haloLayer, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int owner, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int owner, int haloLayer, long id = Element::NULL_ID);
#endif
	bool deleteCell(long id);
	template<typename IdStorage>
	bool deleteCells(const IdStorage &ids);
#if BITPIT_ENABLE_MPI==1
	CellIterator ghostCell2InternalCell(long id);
	CellIterator internalCell2GhostCell(long id, int owner, int haloLayer);
#endif
	virtual double evalCellSize(long id) const = 0;
	BITPIT_DEPRECATED(long countFreeCells() const);
	long countBorderCells() const;
	long countOrphanCells() const;
	std::vector<long> findOrphanCells() const;
	long countDuplicateCells() const;
	std::vector<long> findDuplicateCells() const;

	virtual std::array<double, 3> evalCellCentroid(long id) const;
	virtual void evalCellBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;
	BITPIT_DEPRECATED(ConstProxyVector<std::array<double BITPIT_COMMA 3>> getCellVertexCoordinates(long id) const);
	void getCellVertexCoordinates(long id, std::unique_ptr<std::array<double, 3>[]> *coordinates) const;
	void getCellVertexCoordinates(long id, std::array<double, 3> *coordinates) const;
	std::vector<long> findCellNeighs(long id) const;
	void findCellNeighs(long id, std::vector<long> *neighs) const;
	std::vector<long> findCellNeighs(long id, int codimension, bool complete = true) const;
	void findCellNeighs(long id, int codimension, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellFaceNeighs(long id) const;
	void findCellFaceNeighs(long id, std::vector<long> *neighs) const;
	std::vector<long> findCellFaceNeighs(long id, int face) const;
	void findCellFaceNeighs(long id, int face, std::vector<long> *neighs) const;
	std::vector<long> findCellEdgeNeighs(long id, bool complete = true) const;
	void findCellEdgeNeighs(long id, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellEdgeNeighs(long id, int edge) const;
	void findCellEdgeNeighs(long id, int edge, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexNeighs(long id, bool complete = true) const;
	void findCellVertexNeighs(long id, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexNeighs(long id, int vertex) const;
	void findCellVertexNeighs(long id, int vertex, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexOneRing(long id, int vertex) const;
	void findCellVertexOneRing(long id, int vertex, std::vector<long> *neighs) const;
	bool findFaceNeighCell(long cellId, long neighId, int *cellFace, int *cellAdjacencyId) const;

	std::set<int> getInternalCellPIDs();
	std::vector<long> getInternalCellsByPID(int pid);

	std::vector<long> findVertexOneRing(long vertexId) const;
	void findVertexOneRing(long vertexId, std::vector<long> *ring) const;

	CellIterator getCellIterator(long id);
	CellIterator cellBegin();
	CellIterator cellEnd();
	CellIterator internalCellBegin();
	BITPIT_DEPRECATED(CellIterator internalBegin());
	CellIterator internalCellEnd();
	BITPIT_DEPRECATED(CellIterator internalEnd());
#if BITPIT_ENABLE_MPI==1
	CellIterator ghostCellBegin();
	BITPIT_DEPRECATED(CellIterator ghostBegin());
	CellIterator ghostCellEnd();
	BITPIT_DEPRECATED(CellIterator ghostEnd());
#endif

	CellConstIterator getCellConstIterator(long id) const;
	CellConstIterator cellConstBegin() const;
	CellConstIterator cellConstEnd() const;
	CellConstIterator internalCellConstBegin() const;
	BITPIT_DEPRECATED(CellConstIterator internalConstBegin() const);
	CellConstIterator internalCellConstEnd() const;
	BITPIT_DEPRECATED(CellConstIterator internalConstEnd() const);
#if BITPIT_ENABLE_MPI==1
	CellConstIterator ghostCellConstBegin() const;
	BITPIT_DEPRECATED(CellConstIterator ghostConstBegin() const);
	CellConstIterator ghostCellConstEnd() const;
	BITPIT_DEPRECATED(CellConstIterator ghostConstEnd() const);
#endif

	bool isInterfaceAutoIndexingEnabled() const;
	void setInterfaceAutoIndexing(bool enabled);

	virtual long getInterfaceCount() const;
	PiercedVector<Interface> &getInterfaces();
	const PiercedVector<Interface> & getInterfaces() const;
	Interface &getInterface(long id);
	const Interface &getInterface(long id) const;
	virtual ElementType getInterfaceType(long id) const;
	InterfaceIterator addInterface(const Interface &source, long id = Element::NULL_ID);
	InterfaceIterator addInterface(Interface &&source, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, const std::vector<long> &connectivity, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id = Element::NULL_ID);
	bool deleteInterface(long id);
	template<typename IdStorage>
	bool deleteInterfaces(const IdStorage &ids);
	BITPIT_DEPRECATED(long countFreeInterfaces() const);
	long countBorderInterfaces() const;
	long countOrphanInterfaces() const;
	std::vector<long> findOrphanInterfaces() const;
	bool deleteOrphanInterfaces();
	bool isInterfaceOrphan(long id) const;
	virtual std::array<double, 3> evalInterfaceCentroid(long id) const;
	virtual void evalInterfaceBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;
	BITPIT_DEPRECATED(ConstProxyVector<std::array<double BITPIT_COMMA 3>> getInterfaceVertexCoordinates(long id) const);
	void getInterfaceVertexCoordinates(long id, std::unique_ptr<std::array<double, 3>[]> *coordinates) const;
	void getInterfaceVertexCoordinates(long id, std::array<double, 3> *coordinates) const;

	InterfaceIterator getInterfaceIterator(long id);
	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

	InterfaceConstIterator getInterfaceConstIterator(long id) const;
	InterfaceConstIterator interfaceConstBegin() const;
	InterfaceConstIterator interfaceConstEnd() const;

	long countFaces() const;
	long countBorderFaces() const;
	BITPIT_DEPRECATED(long countFreeFaces() const);

	bool sort();
	bool sortVertices();
	bool sortCells();
	bool sortInterfaces();

	bool squeeze();
	bool squeezeVertices();
	bool squeezeCells();
	bool squeezeInterfaces();

	long locatePoint(double x, double y, double z) const;
	virtual long locatePoint(const std::array<double, 3> &point) const = 0;

	AdjacenciesBuildStrategy getAdjacenciesBuildStrategy() const;
	bool areAdjacenciesDirty(bool global = false) const;
	BITPIT_DEPRECATED(void buildAdjacencies());
	void initializeAdjacencies(AdjacenciesBuildStrategy strategy = ADJACENCIES_AUTOMATIC);
	void updateAdjacencies(bool forcedUpdated = false);
	void destroyAdjacencies();

	InterfacesBuildStrategy getInterfacesBuildStrategy() const;
	bool areInterfacesDirty(bool global = false) const;
	BITPIT_DEPRECATED(void buildInterfaces());
	std::vector<adaption::Info> initializeInterfaces(InterfacesBuildStrategy strategy = INTERFACES_AUTOMATIC, bool trackAdaption = false);
	std::vector<adaption::Info> updateInterfaces(bool forcedUpdated = false, bool trackAdaption = false);
	std::vector<adaption::Info> destroyInterfaces(bool trackAdaption = false);

	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const;
	void getBoundingBox(bool global, std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const;
	bool isBoundingBoxDirty(bool global = false) const;
	void updateBoundingBox(bool forcedUpdated = false);

	virtual void translate(const std::array<double, 3> &translation);
	void translate(double sx, double sy, double sz);
	virtual void rotate(const std::array<double, 3> &n0, const std::array<double, 3> &n1, double angle);
	void rotate(double n0x, double n0y, double n0z, double n1x, double n1y, double n1z, double angle);
	void scale(const std::array<double, 3> &scaling);
	virtual void scale(const std::array<double, 3> &scaling, const std::array<double, 3> &origin);
	void scale(double scaling);
	void scale(double scaling, const std::array<double, 3> &origin);
	void scale(double sx, double sy, double sz);
	void scale(double sx, double sy, double sz, const std::array<double, 3> &origin);

	void setTol(double tolerance);
	double getTol() const;
	void resetTol();
	bool isTolCustomized() const;

	void displayTopologyStats(std::ostream &out, unsigned int padding = 0) const;
	void displayVertices(std::ostream &out, unsigned int padding = 0) const;
	void displayCells(std::ostream &out, unsigned int padding = 0) const;
	void displayInterfaces(std::ostream &out, unsigned int padding = 0) const;

	VTKUnstructuredGrid & getVTK();
	WriteTarget getVTKWriteTarget() const;
	void setVTKWriteTarget(WriteTarget targetCells);
	const CellConstRange getVTKCellWriteRange() const;
	void write(VTKWriteMode mode = VTKWriteMode::DEFAULT);
	void write(const std::string &name, VTKWriteMode mode = VTKWriteMode::DEFAULT);

	void flushData(std::fstream &stream, const std::string &name, VTKFormat format) override;

	int getDumpVersion() const;
	bool dump(std::ostream &stream);
	bool dump(std::ostream &stream) const;
	void restore(std::istream &stream, bool reregister = false);

	void consecutiveRenumberVertices(long offset = 0);
	void consecutiveRenumberCells(long offset = 0);
	void consecutiveRenumberInterfaces(long offset = 0);
	void consecutiveRenumber(long offsetVertices, long offsetCells, long offsetInterfaces);

#if BITPIT_ENABLE_MPI==1
	const MPI_Comm & getCommunicator() const;
	int getRank() const;
	int getProcessorCount() const;

	bool isDistributed(bool allowDirty = false) const;
	int getOwner(bool allowDirty = false) const;

	void setHaloSize(std::size_t haloSize);
	std::size_t getHaloSize() const;

	int getCellOwner(long id) const;
	BITPIT_DEPRECATED(int getCellRank(long id) const);
	int getCellHaloLayer(long id) const;

	int getVertexOwner(long id) const;
	BITPIT_DEPRECATED(int getVertexRank(long id) const);

	bool isRankNeighbour(int rank);
	std::vector<int> getNeighbourRanks();

	const std::unordered_map<int, std::vector<long>> & getGhostVertexExchangeTargets() const;
	const std::vector<long> & getGhostVertexExchangeTargets(int rank) const;
	const std::unordered_map<int, std::vector<long>> & getGhostVertexExchangeSources() const;
	const std::vector<long> & getGhostVertexExchangeSources(int rank) const;

	const std::unordered_map<int, std::vector<long>> & getGhostCellExchangeTargets() const;
	BITPIT_DEPRECATED(const std::unordered_map<int BITPIT_COMMA std::vector<long>> & getGhostExchangeTargets() const);
	const std::vector<long> & getGhostCellExchangeTargets(int rank) const;
	BITPIT_DEPRECATED(const std::vector<long> & getGhostExchangeTargets(int rank) const);
	const std::unordered_map<int, std::vector<long>> & getGhostCellExchangeSources() const;
	BITPIT_DEPRECATED(const std::unordered_map<int BITPIT_COMMA std::vector<long>> & getGhostExchangeSources() const);
	const std::vector<long> & getGhostCellExchangeSources(int rank) const;
	BITPIT_DEPRECATED(const std::vector<long> & getGhostExchangeSources(int rank) const);

	bool isPartitioned() const;
	bool isPartitioningSupported() const;
	bool arePartitioningInfoDirty(bool global = true) const;
	PartitioningStatus getPartitioningStatus(bool global = false) const;
	double evalPartitioningUnbalance() const;
	double evalPartitioningUnbalance(const std::unordered_map<long, double> &cellWeights) const;
	BITPIT_DEPRECATED(std::vector<adaption::Info> partition(MPI_Comm communicator, const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, bool squeezeStorage = false, std::size_t haloSize = 1));
	std::vector<adaption::Info> partition(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, bool squeezeStorage = false);
	BITPIT_DEPRECATED(std::vector<adaption::Info> partition(MPI_Comm communicator, const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, bool squeezeStorage = false, std::size_t haloSize = 1));
	std::vector<adaption::Info> partition(const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, bool squeezeStorage = false);
	BITPIT_DEPRECATED(std::vector<adaption::Info> partition(MPI_Comm communicator, bool trackPartitioning, bool squeezeStorage = false, std::size_t haloSize = 1));
	std::vector<adaption::Info> partition(bool trackPartitioning, bool squeezeStorage = false);
	BITPIT_DEPRECATED(std::vector<adaption::Info> partitioningPrepare(MPI_Comm communicator, const std::unordered_map<long, int> &cellRanks, bool trackPartitioning, std::size_t haloSize = 1));
	std::vector<adaption::Info> partitioningPrepare(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning);
	BITPIT_DEPRECATED(std::vector<adaption::Info> partitioningPrepare(MPI_Comm communicator, const std::unordered_map<long, double> &cellWeights, bool trackPartitioning, std::size_t haloSize = 1));
	std::vector<adaption::Info> partitioningPrepare(const std::unordered_map<long, double> &cellWeights, bool trackPartitioning);
	BITPIT_DEPRECATED(std::vector<adaption::Info> partitioningPrepare(MPI_Comm communicator, bool trackPartitioning, std::size_t haloSize = 1));
	std::vector<adaption::Info> partitioningPrepare(bool trackPartitioning);
	std::vector<adaption::Info> partitioningAlter(bool trackPartitioning = true, bool squeezeStorage = false);
	void partitioningCleanup();
#endif

	std::array<double, 3> evalElementCentroid(const Element &element) const;
	void evalElementBoundingBox(const Element &element, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;
	BITPIT_DEPRECATED(ConstProxyVector<std::array<double BITPIT_COMMA 3>> getElementVertexCoordinates(const Element &element) const);
	void getElementVertexCoordinates(const Element &element, std::unique_ptr<std::array<double, 3>[]> *coordinates) const;
	void getElementVertexCoordinates(const Element &element, std::array<double, 3> *coordinates) const;

protected:
	typedef uint16_t AlterationFlags;
	typedef std::unordered_map<long, AlterationFlags> AlterationFlagsStorage;

#if BITPIT_ENABLE_MPI==1
	const static int DEFAULT_PARTITIONING_WEIGTH;
#endif

	const static AlterationFlags FLAG_NONE              = 0x0;
	const static AlterationFlags FLAG_DELETED           = (1u << 0);
	const static AlterationFlags FLAG_ADJACENCIES_DIRTY = (1u << 1);
	const static AlterationFlags FLAG_INTERFACES_DIRTY  = (1u << 2);
	const static AlterationFlags FLAG_DANGLING          = (1u << 3);

	PiercedVector<Vertex> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	AlterationFlagsStorage m_alteredCells;
	AlterationFlagsStorage m_alteredInterfaces;

#if BITPIT_ENABLE_MPI==1
	PatchKernel(MPI_Comm communicator, std::size_t haloSize, bool expert);
	PatchKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert);
	PatchKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert);
#else
	PatchKernel(bool expert);
	PatchKernel(int dimension, bool expert);
	PatchKernel(int id, int dimension, bool expert);
#endif
	PatchKernel(const PatchKernel &other);
	PatchKernel & operator=(const PatchKernel &other) = delete;

	void clearBoundingBox();
	bool isBoundingBoxFrozen() const;
	void setBoundingBoxFrozen(bool frozen);
	void setBoundingBoxDirty(bool dirty);
	void setBoundingBox(const std::array<double, 3> &minPoint, const std::array<double, 3> &maxPoint);

#if BITPIT_ENABLE_MPI==1
	bool isCommunicatorSet() const;
	virtual void setCommunicator(MPI_Comm communicator);
#endif

#if BITPIT_ENABLE_MPI==1
	CellIterator restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int owner, int haloLayer, long id);
#else
	CellIterator restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);
#endif

	InterfaceIterator restoreInterface(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);

#if BITPIT_ENABLE_MPI==1
	VertexIterator restoreVertex(const std::array<double, 3> &coords, int owner, long id);
#else
	VertexIterator restoreVertex(const std::array<double, 3> &coords, long id);
#endif

	bool deleteVertex(long id);
	template<typename IdStorage>
	bool deleteVertices(const IdStorage &ids);
#if BITPIT_ENABLE_MPI==1
	VertexIterator ghostVertex2InternalVertex(long id);
	VertexIterator internalVertex2GhostVertex(long id, int owner);
#endif

	void dumpVertices(std::ostream &stream) const;
	void restoreVertices(std::istream &stream);

	void dumpCells(std::ostream &stream) const;
	void restoreCells(std::istream &stream);

	void dumpInterfaces(std::ostream &stream) const;
	void restoreInterfaces(std::istream &stream);

	void updateLastInternalVertexId();
#if BITPIT_ENABLE_MPI==1
	void updateFirstGhostVertexId();
#endif

	void updateLastInternalCellId();
#if BITPIT_ENABLE_MPI==1
	void updateFirstGhostCellId();
#endif

	std::unordered_map<long, std::vector<long>> binGroupVertices(const PiercedVector<Vertex> &vertices, int nBins);
	std::unordered_map<long, std::vector<long>> binGroupVertices(int nBins);

	void setAdjacenciesBuildStrategy(AdjacenciesBuildStrategy status);
	void resetAdjacencies();
	void pruneStaleAdjacencies();
	virtual void _resetAdjacencies(bool release);
	virtual void _updateAdjacencies();

	void setInterfacesBuildStrategy(InterfacesBuildStrategy status);
	std::vector<adaption::Info> pruneStaleInterfaces(bool trackAdaption);
	virtual std::vector<adaption::Info> _resetInterfaces(bool trackAdaption, bool release);
	virtual std::vector<adaption::Info> _updateInterfaces(bool trackAdaption);

	bool testCellAlterationFlags(long id, AlterationFlags flags) const;
	AlterationFlags getCellAlterationFlags(long id) const;
	void resetCellAlterationFlags(long id, AlterationFlags flags = FLAG_NONE);
	void setCellAlterationFlags(AlterationFlags flags);
	void setCellAlterationFlags(long id, AlterationFlags flags);
	void unsetCellAlterationFlags(AlterationFlags flags);
	void unsetCellAlterationFlags(long id, AlterationFlags flags);

	bool testInterfaceAlterationFlags(long id, AlterationFlags flags) const;
	AlterationFlags getInterfaceAlterationFlags(long id) const;
	void resetInterfaceAlterationFlags(long id, AlterationFlags flags = FLAG_NONE);
	void setInterfaceAlterationFlags(AlterationFlags flags);
	void setInterfaceAlterationFlags(long id, AlterationFlags flags);
	void unsetInterfaceAlterationFlags(AlterationFlags flags);
	void unsetInterfaceAlterationFlags(long id, AlterationFlags flags);

	bool testAlterationFlags(AlterationFlags availableFlags, AlterationFlags requestedFlags) const;

	void setSpawnStatus(SpawnStatus status);
	virtual std::vector<adaption::Info> _spawn(bool trackAdaption);

	void setAdaptionStatus(AdaptionStatus status);
	virtual std::vector<adaption::Info> _adaptionPrepare(bool trackAdaption);
	virtual std::vector<adaption::Info> _adaptionAlter(bool trackAdaption);
	virtual void _adaptionCleanup();
	virtual bool _markCellForRefinement(long id);
	virtual bool _markCellForCoarsening(long id);
	virtual bool _resetCellAdaptionMarker(long id);
	virtual adaption::Marker _getCellAdaptionMarker(long id);
	virtual bool _enableCellBalancing(long id, bool enabled);

	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	virtual int _getDumpVersion() const = 0;
	virtual void _dump(std::ostream &stream) const = 0;
	virtual void _restore(std::istream &stream) = 0;

	virtual long _getCellNativeIndex(long id) const;

	virtual void _findCellNeighs(long id, const std::vector<long> *blackList, std::vector<long> *neighs) const;
	virtual void _findCellFaceNeighs(long id, int face, const std::vector<long> *blackList, std::vector<long> *neighs) const;
	virtual void _findCellEdgeNeighs(long id, int edge, const std::vector<long> *blackList, std::vector<long> *neighs) const;
	virtual void _findCellVertexNeighs(long id, int vertex, const std::vector<long> *blackList, std::vector<long> *neighs) const;

	void setExpert(bool expert);

	void extractEnvelope(PatchKernel &envelope) const;

	void addPointToBoundingBox(const std::array<double, 3> &point);
	void removePointFromBoundingBox(const std::array<double, 3> &point);
#if BITPIT_ENABLE_MPI==1
	virtual std::size_t _getMaxHaloSize();
	virtual void _setHaloSize(std::size_t haloSize);

	void setPartitioned(bool partitioned);
	void setPartitioningStatus(PartitioningStatus status);
	virtual std::vector<adaption::Info> _partitioningPrepare(const std::unordered_map<long, double> &cellWeights, double defaultWeight, bool trackPartitioning);
	virtual std::vector<adaption::Info> _partitioningPrepare(const std::unordered_map<long, int> &cellRanks, bool trackPartitioning);
	virtual std::vector<adaption::Info> _partitioningAlter(bool trackPartitioning);
	virtual void _partitioningCleanup();

	virtual std::vector<long> _findGhostCellExchangeSources(int rank);
#endif

	template<typename item_t, typename id_t = long>
	std::unordered_map<id_t, id_t> consecutiveItemRenumbering(PiercedVector<item_t, id_t> &container, long offset);

	template<typename item_t, typename id_t = long>
	void mappedItemRenumbering(PiercedVector<item_t, id_t> &container, const std::unordered_map<id_t, id_t> &renumberMap);

	virtual int findAdjoinNeighFace(const Cell &cell, int cellFace, const Cell &neigh) const;
	virtual bool isSameFace(const Cell &cell_A, int face_A, const Cell &cell_B, int face_B) const;

	std::vector<long> getOrderedCellsVertices(const std::vector<long> &cellIds, bool interior, bool ghost) const;
	std::vector<long> getOrderedCellsInterfaces(const std::vector<long> &cellIds) const;

private:
	struct GhostVertexInfo {
		int owner;
	};

	struct GhostCellInfo {
		int owner;
		int haloLayer;
	};

	std::unique_ptr<IndexGenerator<long>> m_vertexIdGenerator;
	std::unique_ptr<IndexGenerator<long>> m_interfaceIdGenerator;
	std::unique_ptr<IndexGenerator<long>> m_cellIdGenerator;

	long m_nInternalVertices;
#if BITPIT_ENABLE_MPI==1
	long m_nGhostVertices;
#endif

	long m_lastInternalVertexId;
#if BITPIT_ENABLE_MPI==1
	long m_firstGhostVertexId;
#endif

	long m_nInternalCells;
#if BITPIT_ENABLE_MPI==1
	long m_nGhostCells;
#endif

	long m_lastInternalCellId;
#if BITPIT_ENABLE_MPI==1
	long m_firstGhostCellId;
#endif

	VTKUnstructuredGrid m_vtk ;
	WriteTarget m_vtkWriteTarget;
	PiercedStorage<long, long> m_vtkVertexMap;

	bool m_boxFrozen;
	bool m_boxDirty;
	std::array<double, 3> m_boxMinPoint;
	std::array<double, 3> m_boxMaxPoint;
	std::array<int, 3> m_boxMinCounter;
	std::array<int, 3> m_boxMaxCounter;

	AdjacenciesBuildStrategy m_adjacenciesBuildStrategy;

	InterfacesBuildStrategy m_interfacesBuildStrategy;

	SpawnStatus m_spawnStatus;

	AdaptionStatus m_adaptionStatus;

	bool m_expert;

	int m_id;
	int m_dimension;

	bool m_toleranceCustom;
	double m_tolerance;

	int m_rank;
	int m_nProcessors;
#if BITPIT_ENABLE_MPI==1
	MPI_Comm m_communicator;
	PartitioningStatus m_partitioningStatus;

	int m_owner;

	std::size_t m_haloSize;

	int m_partitioningCellsTag;
	int m_partitioningVerticesTag;
	bool m_partitioningSerialization;
	std::unordered_map<long, int> m_partitioningOutgoings;
	std::vector<std::pair<int, int>> m_partitioningGlobalExchanges;

	bool m_partitioningInfoDirty;

	std::unordered_map<long, GhostVertexInfo> m_ghostVertexInfo;
	std::unordered_map<int, std::vector<long>> m_ghostVertexExchangeTargets;
	std::unordered_map<int, std::vector<long>> m_ghostVertexExchangeSources;

	std::unordered_map<long, GhostCellInfo> m_ghostCellInfo;
	std::unordered_map<int, std::vector<long>> m_ghostCellExchangeTargets;
	std::unordered_map<int, std::vector<long>> m_ghostCellExchangeSources;

	void setGhostVertexInfo(long id, int owner);
	void unsetGhostVertexInfo(long id);
	void clearGhostVerticesInfo();

	void setGhostCellInfo(long id, int owner, int haloLayer);
	void unsetGhostCellInfo(long id);
	void clearGhostCellsInfo();

	void computeCellHaloLayer(int id);

	std::vector<adaption::Info> _partitioningAlter_deleteGhosts(bool trackPartitioning);

	std::unordered_map<long, int> _partitioningAlter_evalGhostCellOwnershipChanges();
	void _partitioningAlter_applyGhostCellOwnershipChanges(int sendRank, std::unordered_map<long, int> *ghostCellOwnershipChanges);

	std::vector<adaption::Info> _partitioningAlter_sendCells(const std::unordered_set<int> &recvRanks, bool trackPartitioning, std::unordered_map<long, int> *ghostCellOwnershipChanges);
	std::vector<adaption::Info> _partitioningAlter_receiveCells(const std::unordered_set<int> &sendRanks, bool trackPartitioning, std::unordered_map<long, int> *ghostCellOwnershipChanges);

	void setPartitioningInfoDirty(bool dirty);

	void updatePartitioningInfo(bool forcedUpdated = false);

	void updateGhostCellExchangeInfo();

	void updateGhostVertexOwners();
	void updateGhostVertexExchangeInfo();

	void updateOwner();
	int evalOwner() const;

	std::unordered_map<long, int> evaluateExchangeVertexOwners() const;

	template<typename ExcludeList>
	bool confirmCellHaloLayer(const Cell &cell, int haloLayer, const ExcludeList &excludeList = ExcludeList()) const;
#endif

#if BITPIT_ENABLE_MPI==1
	void initialize(MPI_Comm communicator, std::size_t haloSize);
	void initializeHaloSize(std::size_t haloSize);
	void initializeCommunicator(MPI_Comm communicator);

	void freeCommunicator();
#else
	void initialize();
    void initializeSerialCommunicator();
#endif

	std::vector<adaption::Info> finalizeAlterations(bool trackAdaption, bool squeezeStorage = false);

	InterfaceIterator buildCellInterface(Cell *cell_1, int face_1, Cell *cell_2, int face_2, long interfaceId = Element::NULL_ID);

	void _setId(int id);

	bool testElementAlterationFlags(long id, AlterationFlags flags, const AlterationFlagsStorage &flagsStorage) const;
	AlterationFlags getElementAlterationFlags(long id, const AlterationFlagsStorage &flagsStorage) const;
	void resetElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const;
	void setElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const;
	void unsetElementAlterationFlags(AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const;
	void unsetElementAlterationFlags(long id, AlterationFlags flags, AlterationFlagsStorage *flagsStorage) const;

	void mergeAdaptionInfo(std::vector<adaption::Info> &&source, std::vector<adaption::Info> &destination);

	void setRestoredCellAlterationFlags(long id);
	void setAddedCellAlterationFlags(long id);
	void setDeletedCellAlterationFlags(long id);

	void setRestoredInterfaceAlterationFlags(long id);
	void setAddedInterfaceAlterationFlags(long id);
	void setDeletedInterfaceAlterationFlags(long id);

	void dumpVertexAutoIndexing(std::ostream &stream) const;
	void restoreVertexAutoIndexing(std::istream &stream);
	void createVertexIndexGenerator(bool populate);
	void importVertexIndexGenerator(const PatchKernel &source);

	void dumpCellAutoIndexing(std::ostream &stream) const;
	void restoreCellAutoIndexing(std::istream &stream);
	void createCellIndexGenerator(bool populate);
	void importCellIndexGenerator(const PatchKernel &source);

	void dumpInterfaceAutoIndexing(std::ostream &stream) const;
	void restoreInterfaceAutoIndexing(std::istream &stream);
	void createInterfaceIndexGenerator(bool populate);
	void importInterfaceIndexGenerator(const PatchKernel &source);

	VertexIterator _addInternalVertex(const std::array<double, 3> &coords, long id);

	void _restoreInternalVertex(const VertexIterator &iterator, const std::array<double, 3> &coords);
#if BITPIT_ENABLE_MPI==1
	void _restoreGhostVertex(const VertexIterator &iterator, const std::array<double, 3> &coords, int owner);
#endif

	void _deleteInternalVertex(long id);
#if BITPIT_ENABLE_MPI==1
	void _deleteGhostVertex(long id);
#endif

	CellIterator _addInternalCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);
#if BITPIT_ENABLE_MPI==1
	CellIterator _addGhostCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int owner, int haloLayer, long id);
#endif

	void _restoreInternalCell(const CellIterator &iterator, ElementType type, std::unique_ptr<long[]> &&connectStorage);
#if BITPIT_ENABLE_MPI==1
	void _restoreGhostCell(const CellIterator &iterator, ElementType type, std::unique_ptr<long[]> &&connectStorage, int owner, int haloLayer);
#endif

	void _deleteInternalCell(long id);
#if BITPIT_ENABLE_MPI==1
	void _deleteGhostCell(long id);
#endif

	InterfaceIterator _addInterface(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);

	void _restoreInterface(const InterfaceIterator &iterator, ElementType type, std::unique_ptr<long[]> &&connectStorage);

	void _deleteInterface(long id);

	void replaceVTKStreamer(const VTKBaseStreamer *original, VTKBaseStreamer *updated);

	template<typename Selector, typename Function, typename SeedContainer>
	void processCellsNeighbours(const SeedContainer &seeds, int nLayers, Selector isSelected, Function function);

};

}

// Template implementation
#include "patch_kernel.tpp"
#include "patch_kernel_parallel.tpp"

#endif
