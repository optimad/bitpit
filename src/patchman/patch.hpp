//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PATCH_HPP__
#define __PATCHMAN_PATCH_HPP__

/*! \file */

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "output_manager.hpp"
#include "patchman_piercedVector.hpp"
#include "vertex.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace pman {

class Patch {

public:
	Patch(const int &id, const int &dimension);

	virtual ~Patch();

	void reset();
	void reset_vertices();
	void reset_cells();
	void reset_interfaces();
	void reset_output();

	const std::vector<Adaption::Info> update(bool trackAdaption = true);

	void mark_cell_for_refinement(const long &id);
	void mark_cell_for_coarsening(const long &id);
	void enable_cell_balancing(const long &id, bool enabled);

	bool is_dirty() const;
	bool is_output_dirty() const;

	int get_id() const;
	int get_dimension() const;
	bool is_three_dimensional() const;

	std::string get_name() const;
	void set_name(std::string name);

	long get_vertex_count() const;
	PiercedVector<Vertex> &vertices();
	Vertex &get_vertex(const long &id);
	const Vertex & get_vertex(const long &id) const;
	const std::array<double, 3> & get_vertex_coords(const long &id) const;

	long get_cell_count() const;
	PiercedVector<Cell> &cells();
	Cell &get_cell(const long &id);
	const Cell &get_cell(const long &id) const;
	virtual double eval_cell_volume(const long &id) = 0;
	virtual double eval_cell_size(const long &id) = 0;
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
	PiercedVector<Interface> &interfaces();
	Interface &get_interface(const long &id);
	const Interface &get_interface(const long &id) const;
	virtual double eval_interface_area(const long &id) = 0;

	void sort();
	void squeeze();

	void write_mesh();
	void write_mesh(std::string name);
	void write_field(std::string name, int type, std::vector<double> values);
	void write_field(std::string filename, std::string name, int type, std::vector<double> values);
	void write_cell_field(std::string name, std::vector<double> values);
	void write_cell_field(std::string filename, std::string name, std::vector<double> values);
	void write_vertex_field(std::string name, std::vector<double> values);
	void write_vertex_field(std::string filename, std::string name, std::vector<double> values);
	OutputManager & get_output_manager();

protected:
	PiercedVector<Vertex> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

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

	void update_output_manager();

private:
	bool m_dirty;
	bool m_dirty_output;

	int m_id;
	int m_dimension;
	std::string m_name;

	vtkSmartPointer<OutputManager> m_output_manager;

	void set_id(int id);
	void set_dimension(int dimension);

};

}

#endif
