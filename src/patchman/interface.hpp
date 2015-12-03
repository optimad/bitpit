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

class Interface : public Element {

public:
	enum Side {
	    LEFT = 0,
	    RIGHT
	};

	static const int SIDE_COUNT = 2;

	Interface();
	Interface(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Interface(Interface&& other) = default;
	Interface& operator=(Interface&& other) = default;

	bool is_border() const;

	void set_normal(const std::array<double, 3> &normal);
	const std::array<double, 3> & get_normal() const;

	std::array<std::array<double, 3>, 3> eval_rotation_from_cartesian();
	static std::array<std::array<double, 3>, 3> eval_rotation_from_cartesian(std::array<double, 3> &versor);
	std::array<std::array<double, 3>, 3> eval_rotation_to_cartesian();
	static std::array<std::array<double, 3>, 3> eval_rotation_to_cartesian(std::array<double, 3> &versor);
	static std::array<std::array<double, 3>, 3> eval_rotation_transpose(const std::array<std::array<double, 3>, 3> &R);

	void set_area(const double &area);
	const double & get_area() const;

	void set_owner(const long &owner, const int &onwerFace);
	void unset_owner();
	long get_owner() const;
	int get_owner_face() const;

	void set_neigh(const long &neigh, const int &onwerFace);
	void unset_neigh();
	long get_neigh() const;
	int get_neigh_face() const;

	std::array<long, 2> get_owner_neigh() const;

protected:

private:
	double m_area;
	std::array<double, 3> m_normal;

	long m_owner;
	int m_ownerFace;

	long m_neigh;
	int m_neighFace;

	Interface(const Interface &other) = delete;
	Interface& operator = (const Interface &other) = delete;

};

}

#endif
