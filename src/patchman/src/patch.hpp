//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PATCH_HPP__
#define __PATCHMAN_PATCH_HPP__

/*! \file */

#include "cell.hpp"
#include "interface.hpp"
#include "output_manager.hpp"
#include "piercedVector.hpp"
#include "node.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace pman {

#define UNUSED(expr) do { (void)(expr); } while (0)

class Patch {

public:
	enum DataLocation {
		DATA_ON_CELL=0,
		DATA_ON_VERTEX
	};

	Patch(const int &id, const int &dimension);

	virtual ~Patch();

	void reset();
	void reset_vertices();
	void reset_cells();
	void reset_interfaces();
	void reset_output();

	bool update();
	bool update(std::vector<uint32_t> &cellMapping);

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
	PiercedVector<Node> &vertices();
	Node &get_vertex(const long &id);

	long get_cell_count() const;
	PiercedVector<Cell> &cells();
	Cell &get_cell(const long &id);

	long get_interface_count() const;
	PiercedVector<Interface> &interfaces();
	Interface &get_interface(const long &id);

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

	std::array<double, 3> & get_opposite_normal(std::array<double, 3> &normal);

protected:
	PiercedVector<Node> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	virtual std::array<double, 3> & _get_opposite_normal(std::array<double, 3> &normal) = 0;
	virtual bool _update(std::vector<uint32_t> &cellMapping) = 0;
	virtual bool _mark_cell_for_refinement(const long &id) = 0;
	virtual bool _mark_cell_for_coarsening(const long &id) = 0;
	virtual bool _enable_cell_balancing(const long &id, bool enabled) = 0;

	void set_dirty(bool dirty);

	void initialize_output();

private:
	bool m_dirty;
	int m_id;
	int m_dimension;
	std::string m_name;

	vtkSmartPointer<OutputManager> m_output_manager;

	void set_id(int id);
	void set_dimension(int dimension);

};

}

#endif
