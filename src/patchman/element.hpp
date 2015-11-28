//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_ELEMENT_HPP__
#define __PATCHMAN_ELEMENT_HPP__

/*! \file */

#include <cstddef>
#include <memory>
#include <vector>

namespace pman {

class Patch;

class ElementInfo {

public:
	enum Type {
		UNDEFINED = -1,
		POINT,
		LINE,
		TRIANGLE,
		RECTANGLE,
		QUADRANGLE,
		POLYGON,
		TETRAHEDRON,
		BRICK,
		HEXAHEDRON,
		PYRAMID,
		PRISM,
		POLYHEDRON
	};

	Type type;
	int dimension;

	int nVertices;
	int nEdges;
	int nFaces;

	static const ElementInfo undefinedInfo;
	static const ElementInfo pointInfo;
	static const ElementInfo lineInfo;
	static const ElementInfo triangleInfo;
	static const ElementInfo rectangleInfo;
	static const ElementInfo quadrangleInfo;
	static const ElementInfo tetrahedronInfo;
	static const ElementInfo brickInfo;
	static const ElementInfo hexahedronInfo;
	static const ElementInfo pyramidInfo;
	static const ElementInfo prismInfo;

	std::vector<Type> face_type;
	std::vector<std::vector<int>> face_connect;

	std::vector<Type> edge_type;
	std::vector<std::vector<int>> edge_connect;

	ElementInfo();
	ElementInfo(ElementInfo::Type type);

	static const ElementInfo & get_element_info(ElementInfo::Type type);

private:
	void initializeUndefinedInfo();
	void initializePointInfo();
	void initializeLineInfo();
	void initializeTriangleInfo();
	void initializeRectangleInfo();
	void initializeQuadrangleInfo();
	void initializeTetrahedronInfo();
	void initializeBrickInfo();
	void initializeHexahedronInfo();
	void initializePyramidInfo();
	void initializePrismInfo();

};

class Element {

public:
	/*!
		Hasher for the ids.

		Since ids are unique, the hasher can be a function that
		takes an id and cast it to a size_t.

		The hasher is defined as a struct, because a struct can be
		passed as an object into metafunctions (meaning that the type
		deduction for the template paramenters can take place, and
		also meaning that inlining is easier for the compiler). A bare
		function would have to be passed as a function pointer.
		To transform a function template into a function pointer,
		the template would have to be manually instantiated (with a
		perhaps unknown type argument).

	*/
	struct IdHasher {
		/*!
			Function call operator that casts the specified
			value to a size_t.

			\tparam U type of the value
			\param value is the value to be casted
			\result Returns the value casted to a size_t.
		*/
		template<typename U>
		constexpr std::size_t operator()(U&& value) const noexcept
		{
			return static_cast<std::size_t>(std::forward<U>(value));
		}
	};

	Element();
	Element(const long &id);
	Element(const long &id, Patch *patch);

	Element(Element&& other) = default;
	Element& operator=(Element&& other) = default;

	Patch * get_patch() const;
	void set_patch(Patch *patch);
	int get_patch_dimension() const;
	bool is_patch_three_dimensional() const;

	void set_id(const long &id);
	long get_id() const;
	
	void set_type(ElementInfo::Type type);
	ElementInfo::Type get_type() const;

	int get_dimension() const;
	bool is_three_dimensional() const;
	
	void set_connect(std::unique_ptr<long[]> connect);
	void unset_connect();
	const long * get_connect() const;

	void set_centroid(std::array<double, 3> &centroid);
	const std::array<double, 3> & get_centroid() const;

	int get_face_count() const;
	ElementInfo::Type get_face_type(const int &face) const;
	std::vector<int> get_face_local_connect(const int &face) const;

	int get_edge_count() const;
	std::vector<int> get_edge_local_connect(const int &edge) const;

	int get_vertex_count() const;
	int get_vertex(const int &vertex) const;

	double eval_length() const;

	static const long NULL_ELEMENT_ID;

protected:
	static std::array<double, 3> cross_3D(std::array<double, 3> &x, std::array<double, 3> &y);
	static void normalize_3D(std::array<double, 3> &x);
	static void transpose_3D(std::array<std::array<double, 3>, 3> &A);

private:
	Patch *m_patch;

	long m_id;

	ElementInfo::Type m_type;

	std::array<double, 3> m_centroid;
	std::unique_ptr<long[]> m_connect;

	Element(const Element &other) = delete;
	Element& operator = (const Element &other) = delete;

};

}

#endif
