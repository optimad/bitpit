//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PATCH_HPP__
#define __PATCHMAN_PATCH_HPP__

/*! \file */

#include "cell.hpp"
#include "interface.hpp"
#include "output_manager.hpp"
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

	void update();
	void update(const std::vector<uint32_t> &cellMapping);

	void mark_for_refinement(Cell &cell);
	void mark_for_refinement(const int &cellId);

	bool is_dirty() const;

	int get_id() const;
	int get_dimension() const;
	bool is_three_dimensional() const;

	std::string get_name() const;
	void set_name(std::string name);

	int get_vertex_count() const;
	int get_cell_count() const;
	int get_interface_count() const;

	void output_write();
	void output_write(std::string name);
	OutputManager & get_output_manager();

protected:
	std::vector<Node> m_vertices;
	std::vector<Cell> m_cells;
	std::vector<Interface> m_interfaces;

	virtual void _update(const std::vector<uint32_t> &cellMapping) = 0;
	virtual void _mark_for_refinement(Cell &cell) = 0;

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
