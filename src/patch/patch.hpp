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

#ifndef __PATCHMAN_PATCH_HPP__
#define __PATCHMAN_PATCH_HPP__

/*! \file */

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "patchPiercedVector.hpp"
#include "vertex.hpp"

#include <bitpit_IO.hpp>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace pman {

/*!
	\ingroup patch
	@{
*/

class Patch : public bitpit::VTKUnstructuredGrid {

public:
	Patch(const int &id, const int &dimension);

	virtual ~Patch();

	void reset();
	void reset_vertices();
	void reset_cells();
	void reset_interfaces();

	const std::vector<Adaption::Info> update(bool trackAdaption = true);

	void mark_cell_for_refinement(const long &id);
	void mark_cell_for_coarsening(const long &id);
	void enable_cell_balancing(const long &id, bool enabled);

	bool is_dirty() const;

	int get_id() const;
	int get_dimension() const;
	bool is_three_dimensional() const;

	std::string get_name() const;
	void set_name(std::string name);

	long get_vertex_count() const;
	bitpit::PiercedVector<Vertex> &vertices();
	Vertex &get_vertex(const long &id);
	const Vertex & get_vertex(const long &id) const;
	const std::array<double, 3> & get_vertex_coords(const long &id) const;

	long get_cell_count() const;
	bitpit::PiercedVector<Cell> &cells();
	Cell &get_cell(const long &id);
	const Cell &get_cell(const long &id) const;
	virtual double eval_cell_volume(const long &id) = 0;
	virtual double eval_cell_size(const long &id) = 0;
	virtual std::array<double, 3> eval_cell_centroid(const long &id);
	std::vector<long> extract_cell_neighs(const long &id) const;
	std::vector<long> extract_cell_neighs(const long &id, int codimension, bool complete = true) const;
	std::vector<long> extract_cell_face_neighs(const long &id) const;
	std::vector<long> extract_cell_face_neighs(const long &id, const int &face, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extract_cell_edge_neighs(const long &id, bool complete = true) const;
	std::vector<long> extract_cell_edge_neighs(const long &id, const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extract_cell_vertex_neighs(const long &id, bool complete = true) const;
	std::vector<long> extract_cell_vertex_neighs(const long &id, const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extract_cell_vertex_neighs(const long &id, const std::vector<int> &vertices, const std::vector<long> &blackList = std::vector<long>()) const;

	long get_interface_count() const;
	bitpit::PiercedVector<Interface> &interfaces();
	Interface &get_interface(const long &id);
	const Interface &get_interface(const long &id) const;
	virtual double eval_interface_area(const long &id) = 0;
	virtual std::array<double, 3> eval_interface_centroid(const long &id);
	virtual std::array<double, 3> eval_interface_normal(const long &id) = 0;

	void sort();
	void squeeze();

	void write_mesh();
	void write_mesh(std::string name);
	void write_field(std::string name, bitpit::VTKLocation location, const std::vector<double> &values);
	void write_field(std::string filename, std::string name, bitpit::VTKLocation location, const std::vector<double> &values);
	void write_cell_field(std::string name, const std::vector<double> &values);
	void write_cell_field(std::string filename, std::string name, const std::vector<double> &values);
	void write_vertex_field(std::string name, const std::vector<double> &values);
	void write_vertex_field(std::string filename, std::string name, const std::vector<double> &values);

	const bitpit::VTKFieldMetaData getMetaData(std::string name);
	void flushData(std::fstream &stream, bitpit::VTKFormat format, std::string name);

protected:
	bitpit::PiercedVector<Vertex> m_vertices;
	bitpit::PiercedVector<Cell> m_cells;
	bitpit::PiercedVector<Interface> m_interfaces;

	std::deque<long> m_unusedVertexIds;
	std::deque<long> m_unusedInterfaceIds;
	std::deque<long> m_unusedCellIds;

	long create_vertex();
	long create_vertex(const long &id);
	void delete_vertex(const long &id, bool delayed = false);

	long create_interface(ElementInfo::Type type = ElementInfo::UNDEFINED);
	long create_interface(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);
	void delete_interface(const long &id, bool delayed = false);

	long create_cell(bool internal = true, ElementInfo::Type type = ElementInfo::UNDEFINED);
	long create_cell(const long &id, bool internal = true, ElementInfo::Type type = ElementInfo::UNDEFINED);
	void delete_cell(const long &id, bool delayed = false);

	virtual const std::vector<Adaption::Info> _update(bool trackAdaption) = 0;
	virtual bool _mark_cell_for_refinement(const long &id) = 0;
	virtual bool _mark_cell_for_coarsening(const long &id) = 0;
	virtual bool _enable_cell_balancing(const long &id, bool enabled) = 0;

	void set_dirty(bool dirty);

private:
	bool m_dirty;

	int m_id;
	int m_dimension;
	std::string m_name;

	std::unordered_map<std::string, const std::vector<double> *> m_dataFields;
	std::unordered_map<std::string, bitpit::VTKLocation> m_dataLocations;
	std::unordered_map<std::string, bitpit::VTKFieldType> m_dataType;

	void set_id(int id);
	void set_dimension(int dimension);

	std::array<double, 3> eval_element_centroid(const Element &element);
};

/*!
	@}
*/

}

#endif
