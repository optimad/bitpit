//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITP_MESH_VERTEX_HPP__
#define __BITP_MESH_VERTEX_HPP__

/*! \file */

#include <array>
#include <memory>

#include <bitpit_containers.hpp>

/*!
	\ingroup Common
	@{
*/

class Vertex;

bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &out, const Vertex &vertex);
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &in, Vertex &vertex);

class Vertex {

friend bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream&, const Vertex &);
friend bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream&, Vertex &);

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

	bool operator==(const Vertex &other);

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

/*!
	@}
*/

#endif
