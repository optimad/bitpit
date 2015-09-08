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

	~Patch();

	void reset();
	void reset_vertices();
	void reset_cells();
	void reset_interfaces();
	void reset_output();

	bool update();
	bool update(std::vector<uint32_t> &cellMapping);

	void mark_cell_for_refinement(const int &id);

	bool is_dirty() const;

	int get_id() const;
	int get_dimension() const;
	bool is_three_dimensional() const;

	std::string get_name() const;
	void set_name(std::string name);

	int get_vertex_count() const;
	PiercedVector<Node> &vertices();
	Node &get_vertex(const int &id);

	int get_cell_count() const;
	PiercedVector<Cell> &cells();
	Cell &get_cell(const int &id);

	int get_interface_count() const;
	PiercedVector<Interface> &interfaces();
	Interface &get_interface(const int &id);

	void output_write();
	void output_write(std::string name);
	OutputManager & get_output_manager();

	double * get_opposite_normal(double *normal);

protected:
	PiercedVector<Node> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	virtual double * _get_opposite_normal(double *normal) = 0;
	virtual bool _update(std::vector<uint32_t> &cellMapping) = 0;
	virtual bool _mark_cell_for_refinement(const int &id) = 0;

	void set_dirty(bool dirty);

	void output_initialize();

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
