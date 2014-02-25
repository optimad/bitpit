// =================================================================================== //
// Input/Output TOOLS                                                                  //
//                                                                                     //
// Functions for PABLO Input and Output.                                               //
//                                                                                     //
// LIST OF FUNCTIONS                                                                   //
// - writeLocalTree : write on .vtu file local octants and ghosts if filled            //
// =================================================================================== //
// INFO                                                                                //
// =================================================================================== //
// Author   : Edoardo Lombard                                                         //
// Company  : Optimad Engineering srl                                                  //
// Date     : Feb 21, 2014                                                              //
// Version  : v 1.0                                                                    //
//                                                                                     //
// All rights reserved                                                                 //
// =================================================================================== //

// =================================================================================== //
// PRE-COMPILATION INSTRUCTIONS                                                        //
// =================================================================================== //
#ifndef IOFUNCT_HPP_
#define IOFUNCT_HPP_

// includes
#include <algorithm>
#include <string>
#include <functional>
#include <cctype>
#include "Class_Para_Tree.hpp"

// =================================================================================== //
// NAMESPACES                                                                          //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// PROTOTYPES                                                                          //
// =================================================================================== //

typedef vector<vector<uint32_t>	>	u32vector2D;
typedef vector<vector<uint64_t>	>	u64vector2D;


// ----------------------------------------------------------------------------------- //
void writeLocalTree(const u32vector2D & nodes, const u32vector2D & connectivity,
		const u32vector2D & ghostNodes, const u32vector2D & ghostConnectivity,
		const Class_Para_Tree & ptree, string filename);

// ----------------------------------------------------------------------------------- //

#endif
