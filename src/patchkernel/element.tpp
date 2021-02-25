/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#ifndef __BITPIT_ELEMENT_TPP__
#define __BITPIT_ELEMENT_TPP__

namespace bitpit {

/*!
	\class ElementHalfItem
	\ingroup patchelements

	\brief The ElementHalfItem class is the base calss for defining element
	half-edge and half-faces.

	ElementHalfItem is the base calss for defining element half-edge and
	half-faces.
*/

/*!
	Constructor.

	\param element is a reference to the element the owns the item
	\param vertexIds are the ids of the vertices
	\param winding is the winding order of the vertices
*/
template<class DerivedElement>
ElementHalfItem<DerivedElement>::ElementHalfItem(DerivedElement &element, ConstProxyVector<long> &&vertexIds, Winding winding)
    : m_element(element), m_vertexIds(std::move(vertexIds)), m_winding(winding)
{
	// Find the vertex with the lowest id, this will be used as first vertex
	// when iterating over the ids.
	std::size_t nVertices = m_vertexIds.size();

	m_firstVertexId = 0;
	if (nVertices > 0) {
		long smallestVertexId = m_vertexIds[m_firstVertexId];
		for (std::size_t i = 1; i < nVertices; ++i) {
			long vertexId = m_vertexIds[i];
			if (vertexId < smallestVertexId) {
				m_firstVertexId = i;
				smallestVertexId = vertexId;
			}
		}
	}
}

/*!
	Get the element the half-item belongs to.

	\result Returns the element the half-item belongs to.
*/
template<class DerivedElement>
DerivedElement & ElementHalfItem<DerivedElement>::getElement() const
{
    return m_element;
}

/*!
	Get the list of vertex ids of the edge.

	\result Returns the list of vertex ids of the edge.
*/
template<class DerivedElement>
const ConstProxyVector<long> & ElementHalfItem<DerivedElement>::getVertexIds() const
{
    return m_vertexIds;
}

/*!
	Get vertex winding order.

	\result Vertex winding order.
*/
template<class DerivedElement>
typename ElementHalfItem<DerivedElement>::Winding ElementHalfItem<DerivedElement>::getWinding() const
{
	return m_winding;
}

/*!
	Set vertex winding order.

	\param winding is the vertex winding order
*/
template<class DerivedElement>
void ElementHalfItem<DerivedElement>::setWinding(Winding winding)
{
	m_winding = winding;
}

/*!
	Comparison operator for the ElementHalfItem class.

	\param other is the object to be compared with
	\result Returns true if the two half-items define the same half-item,
	false otherwise.
*/
template<class DerivedElement>
bool ElementHalfItem<DerivedElement>::operator==(const ElementHalfItem<DerivedElement> &other) const
{
	const bitpit::ConstProxyVector<long> &vertexIds_1 = m_vertexIds;
	const bitpit::ConstProxyVector<long> &vertexIds_2 = other.m_vertexIds;

	std::size_t nVertices = vertexIds_1.size();
	if (nVertices != vertexIds_2.size()) {
		return false;
	}

	int winding = m_winding * other.m_winding;
	std::size_t offset = other.m_firstVertexId - winding * m_firstVertexId + 2 * nVertices;
	for (std::size_t i = 0; i < nVertices; ++i) {
		std::size_t k1 = i;
		std::size_t k2 = (offset + winding * i) % nVertices;
		if (vertexIds_1[k1] != vertexIds_2[k2]) {
			return false;
		}
	}

	return true;
}

/*!
	Comparison operator for the ElementHalfItem class.

	\param other is the object to be compared with
	\result Returns true if the two half-items define different half-item,
	false otherwise.
*/
template<class DerivedElement>
bool ElementHalfItem<DerivedElement>::operator!=(const ElementHalfItem<DerivedElement> &other) const
{
	return !((*this) == other);
}

/*!
	\class ElementHalfItem::Hasher
	\ingroup patchelements

	\brief The ElementHalfItem::Hasher class allows to create hashes for the
	half-items.
*/

/*!
	Generete the hash for the specified half-item.

	\param item is the half-item
	\result The hash of the specified half-item.
*/
template<class DerivedElement>
std::size_t ElementHalfItem<DerivedElement>::Hasher::operator()(const ElementHalfItem &item) const
{
	const ConstProxyVector<long> &vertexIds = item.m_vertexIds;
	std::size_t nVertices = vertexIds.size();

	std::size_t hash = nVertices;
	if (item.m_winding == ElementHalfItem::WINDING_NATURAL) {
		for (std::size_t i = 0; i < nVertices; ++i) {
			std::size_t k = (item.m_firstVertexId + i) % nVertices;
			utils::hashing::hash_combine(hash, vertexIds[k]);
		}
	} else {
		// Reversing the winding of pixel elements needs special handling. In
		// pixel elements, vertices are ordered using a Z-order, wherase in
		// all other elements ordering of the vertices is anti-clockwise.
		if (item.m_element.getType() != ElementType::VOXEL) {
			for (std::size_t i = nVertices; i > 0; --i) {
				std::size_t k = (item.m_firstVertexId + i) % nVertices;
				utils::hashing::hash_combine(hash, vertexIds[k]);
			}
		} else {
			utils::hashing::hash_combine(hash, vertexIds[(item.m_firstVertexId + 1) % nVertices]);
			utils::hashing::hash_combine(hash, vertexIds[(item.m_firstVertexId + 0) % nVertices]);
			utils::hashing::hash_combine(hash, vertexIds[(item.m_firstVertexId + 3) % nVertices]);
			utils::hashing::hash_combine(hash, vertexIds[(item.m_firstVertexId + 2) % nVertices]);
		}
	}

	return hash;
}

/*!
	\class ElementHalfEdge
	\ingroup patchelements

	\brief The ElementHalfEdge class defines element half-edge items.

	ElementHalfEdge is the class that defines element half-edge items.
*/

/*!
	Constructor.

	\param element is a reference to the element the owns the edge
	\param edge if the local index of the edge
	\param winding is the winding order of the vertexIds
*/
template<class DerivedElement>
ElementHalfEdge<DerivedElement>::ElementHalfEdge(DerivedElement &element, int edge, Winding winding)
    : ElementHalfItem<DerivedElement>(element, element.getEdgeVertexIds(edge), winding),
      m_edge(edge)
{
}

/*!
	Get the local edge index.

	\result Returns the local edge index.
*/
template<class DerivedElement>
int ElementHalfEdge<DerivedElement>::getEdge() const
{
    return m_edge;
}

/*!
	\class ElementHalfFace
	\ingroup patchelements

	\brief The ElementHalfFace class defines element half-faces.

	ElementHalfFace is the class that defines element half-faces. Each face
	can be seen as two half-faces: one belonging to an element and the other
	belonging to the neighbouring element. A half-face is identify by its
	vertexIds and by the winding order of the vertexIds.
*/

/*!
	Constructor.

	\param element is a reference to the element the owns the face
	\param face if the local face of the element
	\param winding is the winding order of the vertexIds
*/
template<class DerivedElement>
ElementHalfFace<DerivedElement>::ElementHalfFace(DerivedElement &element, int face, Winding winding)
    : ElementHalfItem<DerivedElement>(element, element.getFaceVertexIds(face), winding),
      m_face(face)
{
}

/*!
	Get the local face index.

	\result Returns the local face index.
*/
template<class DerivedElement>
int ElementHalfFace<DerivedElement>::getFace() const
{
    return m_face;
}

}

#endif
