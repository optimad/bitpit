// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ClassIntersection.hpp"

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

// =================================================================================== //
// CONSTRUCTORS AND OPERATORS
// =================================================================================== //

ClassIntersection::ClassIntersection(){
	m_owners[0] = 0;
	m_owners[1] = 0;
	m_iface = 0;
	m_isnew = false;
	m_isghost = false;
	m_finer = 0;
	m_bound = m_pbound = false;
	m_dim = 2;
};

ClassIntersection::ClassIntersection(uint8_t dim_){
	m_owners[0] = 0;
	m_owners[1] = 0;
	m_iface = 0;
	m_isnew = false;
	m_isghost = false;
	m_finer = 0;
	m_bound = m_pbound = false;
	m_dim = dim_;
};

ClassIntersection::ClassIntersection(const ClassIntersection & intersection){
	m_owners[0] = intersection.m_owners[0];
	m_owners[1] = intersection.m_owners[1];
	m_iface = intersection.m_iface;
	m_isnew = intersection.m_isnew;
	m_isghost = intersection.m_isghost;
	m_finer = intersection.m_finer;
	m_bound = intersection.m_bound;
	m_pbound = intersection.m_pbound;
	m_dim = intersection.m_dim;
};

ClassIntersection& ClassIntersection::operator =(const ClassIntersection & intersection){
	m_owners[0] = intersection.m_owners[0];
	m_owners[1] = intersection.m_owners[1];
	m_iface = intersection.m_iface;
	m_isnew = intersection.m_isnew;
	m_isghost = intersection.m_isghost;
	m_finer = intersection.m_finer;
	m_bound = intersection.m_bound;
	m_pbound = intersection.m_pbound;
	m_dim = intersection.m_dim;
	return *this;
};

bool ClassIntersection::operator ==(const ClassIntersection & intersection){
	bool check = true;
	check = check && (m_owners[0] == intersection.m_owners[0]);
	check = check && (m_owners[1] == intersection.m_owners[1]);
	check = check && (m_iface == intersection.m_iface);
	check = check && (m_isnew == intersection.m_isnew);
	check = check && (m_isghost == intersection.m_isghost);
	check = check && (m_finer == intersection.m_finer);
	check = check && (m_bound == intersection.m_bound);
	check = check && (m_pbound == intersection.m_pbound);
	check = check && (m_dim == intersection.m_dim);
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
uint32_t ClassIntersection::getOut(){
	return m_owners[m_finer];
};

/*!Get the owner with entering normal;
 */
uint32_t ClassIntersection::getIn(){
	return m_owners[!m_finer];
};

/*!Get the direction of the exiting normal;
 * \param[out] normal Components of the exiting normal.
 * \param[in] normals Basic matrix with components of the elementary normals.
 */
void ClassIntersection::getNormal(int8_t normal[3], int8_t normals[6][3]){
	for (int i=0; i<m_dim; i++){
		normal[i] = normals[m_iface][i];
	}
};

/*!Get the boundary condition of the intersection;
 * \return Boolean true/false if the intersection is/is not a boundary intersection
 */
bool ClassIntersection::getBound(){
	return m_bound;
};

/*!Get the ghost information about the intersection;
 * \return Boolean true/false if the intersection is/is not a ghost intersection
 */
bool ClassIntersection::getIsGhost(){
	return m_isghost;
};

/*!Get the partition boundary condition of the intersection;
 * \return Boolean true/false if the intersection is/is not a process boundary intersection
 */
bool ClassIntersection::getPbound(){
	return m_pbound;
};

