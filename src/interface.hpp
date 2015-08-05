//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_INTERFACE_HPP__
#define __PATCHMAN_INTERFACE_HPP__

/*! \file */

#include <memory>

#include "element.hpp"

namespace pman {

class Cell;

class InterfaceData {};

class Interface : public Element {

public:
	enum PositionType {
	    INTERNAL = 0,
	    BOUNDARY,
	    GHOST
	};

	enum Side {
	    LEFT = 0,
	    RIGHT
	};

	static const int SIDE_COUNT = 2;

	Interface(const int &id = -1);

	void set_position_type(PositionType position);
	PositionType get_position_type() const;

	void set_normal(double normal[]);
	double * get_normal() const;

	void set_area(double *area);
	double get_area() const;

	void set_owner(Cell * neigh, const int &onwerFace);
	void unset_owner();
	Cell * get_owner() const;
	int get_owner_face() const;

	void set_neigh(Cell * neigh, const int &onwerFace);
	void unset_neigh();
	Cell * get_neigh() const;
	int get_neigh_face() const;

	void set_data(std::unique_ptr<InterfaceData> m_data);
	InterfaceData * get_data() const;

protected:

private:
  	PositionType m_positionType;

	double *m_area;
	double *m_normal;

	Cell *m_owner;
	int m_ownerFace;

	Cell *m_neigh;
	int m_neighFace;

	std::unique_ptr<InterfaceData> m_data;

};

}

#endif
