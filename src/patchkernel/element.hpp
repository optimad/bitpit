/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __BITPIT_ELEMENT_HPP__
#define __BITPIT_ELEMENT_HPP__

#include <cstddef>
#include <memory>
#include <vector>

#include "bitpit_containers.hpp"

namespace bitpit {
	class Element;
}

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, bitpit::Element& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const bitpit::Element& element);

namespace bitpit {

class ElementInfo {

public:
	enum Type {
		UNDEFINED  = -1,
		VERTEX     = 1,
		LINE       = 3,
		TRIANGLE   = 5,
		POLYGON    = 7,
		PIXEL      = 8,
		QUAD       = 9,
		TETRA      = 10,
		VOXEL      = 11,
		HEXAHEDRON = 12,
		WEDGE      = 13,
		PYRAMID    = 14,
		POLYHEDRON = 42
	};

	Type type;
	int dimension;

	int nVertices;
	int nEdges;
	int nFaces;

	static const ElementInfo undefinedInfo;
	static const ElementInfo vertexInfo;
	static const ElementInfo lineInfo;
	static const ElementInfo triangleInfo;
	static const ElementInfo pixelInfo;
	static const ElementInfo quadInfo;
	static const ElementInfo tetraInfo;
	static const ElementInfo voxelInfo;
	static const ElementInfo hexahedronInfo;
	static const ElementInfo pyramidInfo;
	static const ElementInfo wedgeInfo;

	std::vector<Type> face_type;
	std::vector<std::vector<int>> faceConnect;
	std::vector<std::vector<int>> faceEdges;

	std::vector<Type> edge_type;
	std::vector<std::vector<int>> edgeConnect;

	ElementInfo();
	ElementInfo(ElementInfo::Type type);

	static const ElementInfo & getElementInfo(ElementInfo::Type type);

private:
	void initializeUndefinedInfo();
	void initializeVertexInfo();
	void initializeLineInfo();
	void initializeTriangleInfo();
	void initializePixelInfo();
	void initializeQuadInfo();
	void initializeTetraInfo();
	void initializeVoxelInfo();
	void initializeHexahedronInfo();
	void initializePyramidInfo();
	void initializeWedgeInfo();

	void initializeFaceEdges(const std::vector<ElementInfo *> &facesInfo = std::vector<ElementInfo *>());
};

class Element {

friend bitpit::OBinaryStream& (::operator<<) (bitpit::OBinaryStream& buf, const Element& element);
friend bitpit::IBinaryStream& (::operator>>) (bitpit::IBinaryStream& buf, Element& element);

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
	Element(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Element(const Element &other);
	Element(Element&& other) = default;
	Element& operator = (const Element &other);
	Element& operator=(Element&& other) = default;

	void swap(Element &other) noexcept;

	void initialize(ElementInfo::Type type);

	const ElementInfo & getInfo() const;

	void setId(const long &id);
	long getId() const;
	
	void setType(ElementInfo::Type type);
	ElementInfo::Type getType() const;

	int getDimension() const;
	bool isThreeDimensional() const;
	
	void setConnect(std::unique_ptr<long[]> &&connect);
	void unsetConnect();
	const long * getConnect() const;
	long * getConnect();

	int getFaceCount() const;
	ElementInfo::Type getFaceType(const int &face) const;
	const std::vector<int> & getFaceLocalConnect(const int &face) const;
	std::vector<long> getFaceConnect(int face) const;

	int getEdgeCount() const;
	const std::vector<int> & getEdgeLocalConnect(const int &edge) const;
	std::vector<long> getEdgeConnect(int edge) const;

	void setVertex(const int &index, const long &vertex);
	int getVertexCount() const;
	long getVertex(const int &vertex) const;

	static const long NULL_ID;

	unsigned int getBinarySize();

private:
	long m_id;

	ElementInfo::Type m_type;

	std::unique_ptr<long[]> m_connect;

	void _initialize(ElementInfo::Type type);

};

extern template class PiercedVector<Element>;

}

#endif
