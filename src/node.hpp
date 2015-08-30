//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_NODE_HPP__
#define __PATCHMAN_NODE_HPP__

/*! \file */

#include <memory>

namespace pman {

class Node {

public:
	enum Coordinate {
		COORD_X =0,
		COORD_Y,
		COORD_Z,
	};

	Node();
	Node(const int &id);

	int get_id() const;

	void set_coords(std::unique_ptr<double[]> coords);
	double * get_coords() const;
	
protected:
	static const int NULL_NODE_ID;

private:
	int m_id;
	std::unique_ptr<double[]> m_coords;

	void set_id(const int &id);

};

}

#endif
