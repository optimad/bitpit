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

class Cell : public Element {

public:
	Cell();
	Cell(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Cell(Cell&& other) = default;
	Cell& operator=(Cell&& other) = default;

	void initialize(ElementInfo::Type type, int nInterfacesPerFace = 0);

	void set_interior(bool interior);
	bool is_interior() const;
	
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

protected:

private:
	bool m_interior;

	double m_volume;

	CollapsedVector2D<long> m_interfaces;

	Cell(const Cell &other) = delete;
	Cell& operator = (const Cell &other) = delete;

};

}

#endif
