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

#ifndef __BITPIT_PATCH_KERNEL_HPP__
#define __BITPIT_PATCH_KERNEL_HPP__

#include <cstddef>
#include <deque>
#include <memory>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_IO.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class IndexGenerator {

public:
	IndexGenerator();

	long generateId();
	long getLastId();
	void trashId(const long &id);
	void reset();

private:
	long m_id;
	std::deque<long> m_trash;

};

template<class T, class T1>
class UnaryPredicate;

template<class T, class T1>
bool operator==(const std::pair<T, T1> &pair_, const UnaryPredicate<T, T1> &pred_) { return ( pred_.value == pair_.second ); }

template<class T, class T1>
class UnaryPredicate {
    private:
    T                                   value;
    public:
    UnaryPredicate(
        T                               value_
    ) : value(value_) {}
    template<class U, class U1>
    friend bool operator==(const std::pair<U, U1>&, const UnaryPredicate<U, U1>&);
};

class PatchKernel : public VTKUnstructuredGrid {

public:
	typedef PiercedVector<Vertex>::iterator VertexIterator;
	typedef PiercedVector<Cell>::iterator CellIterator;
	typedef PiercedVector<Interface>::iterator InterfaceIterator;

	PatchKernel(const int &id, const int &dimension, bool epxert);

	virtual ~PatchKernel();

	void reset();
	void resetVertices();
	void resetCells();
	void resetInterfaces();

	bool reserveVertices(size_t nVertices);
	bool reserveCells(size_t nCells);
	bool reserveInterfaces(size_t nInterfaces);

	const std::vector<Adaption::Info> update(bool trackAdaption = true);

	void markCellForRefinement(const long &id);
	void markCellForCoarsening(const long &id);
	void enableCellBalancing(const long &id, bool enabled);

	bool isDirty() const;
	bool isExpert() const;

	int getId() const;
	int getDimension() const;
	bool isThreeDimensional() const;

	virtual long getVertexCount() const;
	PiercedVector<Vertex> &getVertices();
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

	virtual long getCellCount() const;
	long getInternalCount() const;
	long getGhostCount() const;
	PiercedVector<Cell> &getCells();
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
	virtual double evalCellSize(const long &id) = 0;
	long countFreeCells() const;
	long countOrphanCells() const;
	virtual std::array<double, 3> evalCellCentroid(const long &id);
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

	CellIterator getCellIterator(const long &id);
	CellIterator cellBegin();
	CellIterator cellEnd();
	CellIterator internalBegin();
	CellIterator internalEnd();
	CellIterator ghostBegin();
	CellIterator ghostEnd();

	virtual long getInterfaceCount() const;
	PiercedVector<Interface> &getInterfaces();
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
	virtual std::array<double, 3> evalInterfaceCentroid(const long &id);

	InterfaceIterator getInterfaceIterator(const long &id);
	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

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

	bool isPointInside(const double &x, const double &y, const double &z);
	virtual bool isPointInside(const std::array<double, 3> &point) = 0;
	long locatePoint(const double &x, const double &y, const double &z);
	virtual long locatePoint(const std::array<double, 3> &point) = 0;
        bool isSameFace(const long &, const int&, const long&, const int&);

        virtual void buildAdjacencies() = 0;
        virtual void updateAdjacencies(const std::vector<long>&) = 0;

	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint);
	bool isBoundingBoxDirty() const;
	void updateBoundingBox(bool forcedUpdated = false);
	void addPointToBoundingBox(const std::array<double, 3> &point);
	void removePointFromBoundingBox(const std::array<double, 3> &point, bool delayedBoxUpdate = false);

	std::unordered_map<long, long> binSortVertex(int nBins = 128);

	bool isAdaptionDirty() const;
	const std::vector<Adaption::Info> updateAdaption(bool trackAdaption = true);

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

	void write();
	void write(std::string name);

	const VTKFieldMetaData getMetaData(std::string name);
	void flushData(std::fstream &stream, VTKFormat format, std::string name);

#if BITPIT_ENABLE_MPI==1
	void setCommunicator(MPI_Comm communicator);
	bool isCommunicatorSet() const;
	const MPI_Comm & getCommunicator() const;
	int getRank() const;
	int getProcessorCount() const;

	std::unordered_map<short, std::unordered_map<long, long>> & getGhostExchangeData();
	const std::unordered_map<short, std::unordered_map<long, long>> & getGhostExchangeData() const;
	std::unordered_map<long, long> & getGhostExchangeData(short rank);
	const std::unordered_map<long, long> & getGhostExchangeData(short rank) const;

	const std::vector<Adaption::Info> partition(const std::vector<int> &cellRanks, bool trackChanges);
	Adaption::Info sendCells(const unsigned short &, const unsigned short &, const std::vector<long> &);
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

	void clearBoundingBox();
	bool isBoundingBoxFrozen() const;
	void setBoundingBoxFrozen(bool frozen);
	void setBoundingBoxDirty(bool dirty);
	void setBoundingBox(const std::array<double, 3> &minPoint, const std::array<double, 3> &maxPoint);

	bool deleteVertex(const long &id, bool delayed = false);
	bool deleteVertices(const std::vector<long> &ids, bool delayed = false);

	virtual const std::vector<Adaption::Info> _updateAdaption(bool trackAdaption) = 0;
	virtual bool _markCellForRefinement(const long &id) = 0;
	virtual bool _markCellForCoarsening(const long &id) = 0;
	virtual bool _enableCellBalancing(const long &id, bool enabled) = 0;
	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	virtual std::vector<long> _findCellFaceNeighs(const long &id, const int &face, const std::vector<long> &blackList = std::vector<long>()) const;
	virtual std::vector<long> _findCellEdgeNeighs(const long &id, const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;
	virtual std::vector<long> _findCellVertexNeighs(const long &id, const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;

	void setAdaptionDirty(bool dirty);
	void setExpert(bool expert);

#if BITPIT_ENABLE_MPI==1
	void resetGhostExchangeData();
	void resetGhostExchangeData(short rank);

	void setGhostExchangeData(const std::unordered_map<short, std::unordered_map<long, long>> &ghostInfo);
	void setGhostExchangeData(short rank, const std::unordered_map<long, long> &rankGhostInfo);
#endif

private:
	double DEFAULT_TOLERANCE = 1e-14;

	bool m_boxFrozen;
	bool m_boxDirty;
	std::array<double, 3> m_boxMinPoint;
	std::array<double, 3> m_boxMaxPoint;
	std::array<int, 3> m_boxMinCounter;
	std::array<int, 3> m_boxMaxCounter;

	bool m_adaptionDirty;

	bool m_expert;

	int m_id;
	int m_dimension;

	bool m_hasCustomTolerance;
	double m_tolerance;

	int m_rank;
	int m_nProcessors;
#if BITPIT_ENABLE_MPI==1
	MPI_Comm m_communicator;
        std::unordered_map<short, std::unordered_map<long, long> > m_ghost2id;
#endif

	VertexIterator createVertex(const std::array<double, 3> &coords, long id = Vertex::NULL_ID);
	InterfaceIterator createInterface(ElementInfo::Type type, long id = Element::NULL_ID);
	CellIterator createCell(ElementInfo::Type type, bool interior, long id = Element::NULL_ID);

	void setId(int id);
	void setDimension(int dimension);

	std::array<double, 3> evalElementCentroid(const Element &element);
};

}

#endif
