/*
 * Class_Local_Tree.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: edoardo
 */

#include "Class_Local_Tree.hpp"

Class_Local_Tree::Class_Local_Tree() {
	// TODO Auto-generated constructor stub
	Class_Octant oct0(0,0,0,0);
	Class_Octant octf(0,0,0,MAX_LEVEL);
	Class_Octant octl(maxlength-1,maxlength-1,maxlength-1,MAX_LEVEL);
	octants.resize(1);
	octants[0] = oct0;
	first_desc = octf;
	last_desc = octl;
	sizeghosts = 0;

}

Class_Local_Tree::~Class_Local_Tree() {
	// TODO Auto-generated destructor stub
}

//-------------------------------------------------------------------------------- //
// Other methods ----------------------------------------------------------------- //

void Class_Local_Tree::Refine() {


}
