//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_VERTEX_HPP__
#define __PATCHMAN_VERTEX_HPP__

/*! \file */

#include <memory>

#include "binary_stream.hpp"

namespace pman {

class Vertex;

obinarystream& operator<<(obinarystream &out, const Vertex &vertex);
ibinarystream& operator>>(ibinarystream &in, Vertex &vertex);

class Vertex {

friend obinarystream& operator<<(obinarystream&, const Vertex &);
friend ibinarystream& operator>>(ibinarystream&, Vertex &);

public:
	enum Coordinate {
		COORD_X =0,
		COORD_Y,
		COORD_Z,
	};

	Vertex();
	Vertex(const long &id);
	Vertex(const long &id, std::array<double, 3> &coords);

	Vertex(Vertex &&other) = default;
	Vertex& operator=(Vertex &&other) = default;

	double & operator[](int coord_id);
	const double & operator[](int coord_id) const;

	void set_id(const long &id);
	long get_id() const;

	void set_coords(std::array<double, 3> &coords);
	const std::array<double, 3> & get_coords() const;

	static const long NULL_VERTEX_ID;

private:
	long m_id;
	std::array<double, 3> m_coords;

};

}

#endif
