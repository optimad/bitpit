// =================================================================================== //
// PRE-COMPILATION INSTRUCTIONS                                                        //
// =================================================================================== //
#ifndef IOFUNCT_HPP_
#define IOFUNCT_HPP_

// includes
#if NOMPI == 0
#include <mpi.h>
#endif
#include <algorithm>
#include <string>
#include <functional>
#include <cctype>
#include <fstream>
#include <iomanip>
#include "Class_Para_Tree.hpp"
#include "logFunct.hpp"

// =================================================================================== //
// NAMESPACES                                                                          //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// PROTOTYPES                                                                          //
// =================================================================================== //

typedef vector<vector<uint32_t>	>	u32vector2D;
typedef vector<vector<uint64_t>	>	u64vector2D;
typedef vector<vector<double>	>	dvector2D;


// ----------------------------------------------------------------------------------- //
void writeLocalTree(const u32vector2D & nodes, const u32vector2D & connectivity,
		const u32vector2D & ghostNodes, const u32vector2D & ghostConnectivity,
		const Class_Para_Tree<2> & ptree, string filename);

// ----------------------------------------------------------------------------------- //

void writePhysicalTree(const dvector2D & nodes, const u32vector2D & connectivity,
		const dvector2D & ghostNodes, const u32vector2D & ghostConnectivity,
		const Class_Para_Tree<2> & ptree, string filename);

// ----------------------------------------------------------------------------------- //

void writeLocalTree(const u32vector2D & nodes, const u32vector2D & connectivity,
		const u32vector2D & ghostNodes, const u32vector2D & ghostConnectivity,
		const Class_Para_Tree<3> & ptree, string filename);

// ----------------------------------------------------------------------------------- //

void writePhysicalTree(const dvector2D & nodes, const u32vector2D & connectivity,
		const dvector2D & ghostNodes, const u32vector2D & ghostConnectivity,
		const Class_Para_Tree<3> & ptree, string filename);

// ----------------------------------------------------------------------------------- //

#endif
