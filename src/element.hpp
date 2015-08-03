//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_ELEMENT_HPP__
#define __PATCHMAN_ELEMENT_HPP__

/*! \file */

#include <cstddef>
#include <memory>

#include "node.hpp"
#include "collapsedArrayArray.hpp"

class Element {

public:
	enum Type {
	    POINT = 0,
	    LINE,
	    TRIANGLE,
	    QUADRANGLE,
	    POLYGON,
	    TETRAHEDRON,
	    HEXAHEDRON,
	    PYRAMID,
	    PRISM,
	    POLYHEDRON
	};

	Element();
	Element(const int &id = -1);

	int get_id() const;
	
	void set_local_id(int id);
	int get_local_id() const;
	
	void set_type(Element::Type type);
	Element::Type get_type() const;
	
	void set_connect(std::unique_ptr<Node*[]> connect);
	void unset_connect();
	Node ** get_connect() const;

	int get_face_count() const;
	static int get_face_count(Element::Type type);

	int get_vertex_count() const;
	static int get_vertex_count(Element::Type type);

protected:

private:
	int m_id;
	int m_local_id;

	Element::Type m_type;
	
	std::unique_ptr<Node*[]> m_connect;

	void set_id(const int &id);

};

#endif
