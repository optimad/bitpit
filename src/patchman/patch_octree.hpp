//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef DISABLE_OCTREE
#ifndef __PATCHMAN_PATCH_OCTREE_HPP__
#define __PATCHMAN_PATCH_OCTREE_HPP__

/*! \file */

#include "patch.hpp"

#include "Class_Para_Tree.hpp"

#include <assert.h>
#include <deque>
#include <cstddef>
#include <vector>
#include <unordered_set>

namespace pman {

struct OctreeLevelInfo{
    int    level;
    double h;
    double area;
    double volume;
};

class PatchOctree : public Patch {

public:
	struct OctantInfo {
		OctantInfo() : id(0), internal(true) {};
		OctantInfo(uint32_t _id, bool _internal) : id(_id), internal(_internal) {};

		uint32_t id;
		bool internal;
	};

	PatchOctree(const int &id, const int &dimension, std::array<double, 3> origin,
			double length, double dh);

	~PatchOctree();

	double eval_cell_volume(const long &id);

	OctantInfo get_cell_octant(const long &id) const;
	int get_cell_level(const long &id);

	long get_octant_id(const OctantInfo &octantInfo) const;
	const std::vector<uint32_t> & get_octant_connect(const OctantInfo &octantInfo);

	/*!
		\brief Gets the octree associated with the patch.

		\tparam dimension is the dimension of the octree. It MUST match the dimension
		of the patch.
		\result A pointer to the octree associated to the patch.
	*/
	template<int dimension>
	Class_Para_Tree<dimension> * get_tree()
	{
		assert(dimension == get_dimension());

		if (dimension == 2) {
			return dynamic_cast<Class_Para_Tree<dimension> *>(&m_tree_2D);
		} else {
			return dynamic_cast<Class_Para_Tree<dimension> *>(&m_tree_3D);
		}
	}

protected:
	std::array<double, 3> & _get_opposite_normal(std::array<double, 3> &normal);
	const std::vector<Adaption::Info> _update(bool trackAdaption);
	bool _mark_cell_for_refinement(const long &id);
	bool _mark_cell_for_coarsening(const long &id);
	bool _enable_cell_balancing(const long &id, bool enabled);

private:
	typedef std::bitset<72> OctantHash;

	struct FaceInfo {
		FaceInfo() : id(Element::NULL_ELEMENT_ID), face(-1) {};
		FaceInfo(long _id, int _face) : id(_id), face(_face) {};

		bool operator==(const FaceInfo &other) const
		{
			return (id == other.id && face == other.face);
		}

		long id;
		int face;
	};

	struct FaceInfoHasher
	{
		std::size_t operator()(const FaceInfo& k) const
		{
			using std::hash;
			using std::string;

			return ((hash<long>()(k.id) ^ (hash<int>()(k.face) << 1)) >> 1);
		}
	};

	typedef std::unordered_set<FaceInfo, FaceInfoHasher> FaceInfoSet;

	std::unordered_map<long, uint32_t, Element::IdHasher> m_cell_to_octant;
	std::unordered_map<long, uint32_t, Element::IdHasher> m_cell_to_ghost;
	std::unordered_map<uint32_t, long> m_octant_to_cell;
	std::unordered_map<uint32_t, long> m_ghost_to_cell;

	Class_Para_Tree<2> m_tree_2D;
	Class_Para_Tree<3> m_tree_3D;

	vector<double> m_tree_dh;
	vector<double> m_tree_area;
	vector<double> m_tree_volume;

	std::vector<std::array<double, 3> > m_normals;

	bool set_marker(const long &id, const int8_t &value);

	OctantHash evaluate_octant_hash(const OctantInfo &octantInfo);

	std::vector<unsigned long> import_octants(std::vector<OctantInfo> &octantTreeIds);
	std::vector<unsigned long> import_octants(std::vector<OctantInfo> &octantTreeIds, FaceInfoSet &danglingInfoSet);

	FaceInfoSet remove_cells(std::vector<long> &cellIds);

	long create_vertex(uint32_t treeId);

	long create_interface(uint32_t treeId,
                            std::unique_ptr<long[]> &vertices,
                            std::array<FaceInfo, 2> &faces);

	long create_cell(uint32_t treeId, bool interior,
	                 std::unique_ptr<long[]> &vertices,
	                 std::vector<std::vector<long>> &interfaces,
	                 std::vector<std::vector<bool>> &ownerFlags);
	void delete_cell(long id);
};

}

#endif
#endif