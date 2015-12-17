// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "classIntersection.hpp"

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

classIntersection::classIntersection(){
	owners[0] = 0;
	owners[1] = 0;
	iface = 0;
	isnew = false;
	isghost = false;
	finer = 0;
	bound = pbound = false;
	dim = 2;
};

classIntersection::classIntersection(uint8_t dim_){
	owners[0] = 0;
	owners[1] = 0;
	iface = 0;
	isnew = false;
	isghost = false;
	finer = 0;
	bound = pbound = false;
	dim = dim_;
};

classIntersection::classIntersection(const classIntersection & intersection){
	owners[0] = intersection.owners[0];
	owners[1] = intersection.owners[1];
	iface = intersection.iface;
	isnew = intersection.isnew;
	isghost = intersection.isghost;
	finer = intersection.finer;
	bound = intersection.bound;
	pbound = intersection.pbound;
	dim = intersection.dim;
};

classIntersection& classIntersection::operator =(const classIntersection & intersection){
	owners[0] = intersection.owners[0];
	owners[1] = intersection.owners[1];
	iface = intersection.iface;
	isnew = intersection.isnew;
	isghost = intersection.isghost;
	finer = intersection.finer;
	bound = intersection.bound;
	pbound = intersection.pbound;
	dim = intersection.dim;
	return *this;
};

bool classIntersection::operator ==(const classIntersection & intersection){
	bool check = true;
	check = check && (owners[0] == intersection.owners[0]);
	check = check && (owners[1] == intersection.owners[1]);
	check = check && (iface == intersection.iface);
	check = check && (isnew == intersection.isnew);
	check = check && (isghost == intersection.isghost);
	check = check && (finer == intersection.finer);
	check = check && (bound == intersection.bound);
	check = check && (pbound == intersection.pbound);
	check = check && (dim == intersection.dim);
	return check;

};

// =================================================================================== //
// METHODS
// =================================================================================== //

// =================================================================================== //
// BASIC GET/SET METHODS
// =================================================================================== //


/*!Get the owner with exiting normal;
 */
uint32_t classIntersection::getOut(){
	return owners[finer];
};

/*!Get the owner with entering normal;
 */
uint32_t classIntersection::getIn(){
	return owners[!finer];
};

void classIntersection::getNormal(int8_t normal[3]){
	for (int i=0; i<dim; i++){
		normal[i] = CG::normals[iface][i];
	}
};

bool classIntersection::getBound(){
	return bound;
};

bool classIntersection::getIsGhost(){
	return isghost;
};

bool classIntersection::getPbound(){
	return pbound;
};

