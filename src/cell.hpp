//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_CELL_HPP__
#define __PATCHMAN_CELL_HPP__

/*! \file */

#include "element.hpp"

#include <memory>

namespace pman {

class Interface;

class CellData {};

class Cell : public Element {

public:
	enum PositionType {
	    INTERNAL = 0,
	    GHOST
	};

	Cell();
	Cell(const int &id);

	Cell(Cell&& other) = default;
	Cell& operator=(Cell&& other) = default;

	void set_position_type(PositionType position);
	PositionType get_position_type() const;
	
	void set_volume(double * const volume);
	double & get_volume() const;

	void set_centroid(double * const centroid);
	double * get_centroid() const;

	void initialize_interfaces(std::vector<std::vector<Interface *>> &interfaces);
	void initialize_empty_interfaces(const int nInterfaces[]);
	void set_interface(const int &face, const int &index, Interface *interface);
	void set_interfaces(const int &face, Interface **interface);
	void unset_interfaces();
	int get_interface_count(const int &face) const;
	Interface * get_interface(const int &face, const int &index = 0) const;
	Interface ** get_interfaces(const int &face) const;

	void set_data(std::unique_ptr<CellData> m_data);
	CellData * get_data() const;

protected:

private:
  	PositionType m_positionType;

	double *m_volume;
	double *m_centroid;

	std::unique_ptr<CollapsedArrayArray<Interface *> > m_interfaces;

	std::unique_ptr<CellData> m_data;

	Cell(const Cell &other) = delete;
	Cell& operator = (const Cell &other) = delete;

};

}

#endif
