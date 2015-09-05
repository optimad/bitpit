//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_INTERFACE_HPP__
#define __PATCHMAN_INTERFACE_HPP__

/*! \file */

#include <array>
#include <memory>

#include "element.hpp"

namespace pman {

class Cell;
class Patch;

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

	Interface();
	Interface(const int &id);
	Interface(const int &id, Patch *patch);

	Interface(Interface&& other) = default;
	Interface& operator=(Interface&& other) = default;

	void set_position_type(PositionType position);
	PositionType get_position_type() const;

	void set_normal(double normal[]);
	double * get_normal() const;

	void set_area(double *area);
	double get_area() const;

	void set_owner(const int &owner, const int &onwerFace);
	void unset_owner();
	int get_owner() const;
	int get_owner_face() const;

	void set_neigh(const int &neigh, const int &onwerFace);
	void unset_neigh();
	int get_neigh() const;
	int get_neigh_face() const;

	std::array<int, 2> get_owner_neigh() const;
	void swap_owner_neigh();

	void set_data(std::unique_ptr<InterfaceData> m_data);
	InterfaceData * get_data() const;

protected:

private:
  	PositionType m_positionType;

	double *m_area;
	double *m_normal;

	int m_owner;
	int m_ownerFace;

	int m_neigh;
	int m_neighFace;

	std::unique_ptr<InterfaceData> m_data;

	Interface(const Interface &other) = delete;
	Interface& operator = (const Interface &other) = delete;

};

}

#endif
