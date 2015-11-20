//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_CELL_HPP__
#define __PATCHMAN_CELL_HPP__

/*! \file */

#include "collapsedVector2D.hpp"
#include "element.hpp"

#include <memory>

namespace pman {

class Interface;
class Patch;

class Cell : public Element {

public:
	enum PositionType {
	    INTERNAL = 0,
	    GHOST
	};

	Cell();
	Cell(const long &id);
	Cell(const long &id, Patch *patch);

	Cell(Cell&& other) = default;
	Cell& operator=(Cell&& other) = default;

	void set_position_type(PositionType position);
	PositionType get_position_type() const;
	
	void set_volume(double * const volume);
	const double & get_volume() const;

	void initialize_interfaces(std::vector<std::vector<long>> &interfaces);
	void initialize_empty_interfaces(const int nInterfaces[]);
	void set_interface(const int &face, const int &index, const long &interface);
	void push_interface(const int &face, const long &interface);
	void delete_interface(const int &face, const int &i);
	void unset_interfaces();
	int get_interface_count() const;
	int get_interface_count(const int &face) const;
	long get_interface(const int &face, const int &index = 0) const;
	const long * get_interfaces() const;
	const long * get_interfaces(const int &face) const;

	std::vector<long> extract_neighs() const;
	std::vector<long> extract_neighs(int codimension, bool complete = true) const;

	std::vector<long> extract_face_neighs() const;
	std::vector<long> extract_face_neighs(const int &face, const std::vector<long> &blackList = std::vector<long>()) const;

	std::vector<long> extract_edge_neighs(bool complete = true) const;
	std::vector<long> extract_edge_neighs(const int &edge, const std::vector<long> &blackList = std::vector<long>()) const;

	std::vector<long> extract_vertex_neighs(bool complete = true) const;
	std::vector<long> extract_vertex_neighs(const int &vertex, const std::vector<long> &blackList = std::vector<long>()) const;
	std::vector<long> extract_vertex_neighs(const std::vector<int> &vertices, const std::vector<long> &blackList = std::vector<long>()) const;

protected:

private:
  	PositionType m_positionType;

	double *m_volume;

	CollapsedVector2D<long> m_interfaces;

	Cell(const Cell &other) = delete;
	Cell& operator = (const Cell &other) = delete;

};

}

#endif
