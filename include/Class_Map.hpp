#ifndef CLASS_MAP_HPP_
#define CLASS_MAP_HPP_

#include "Class_Global.hpp"
#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <map>
#include <iostream>
#include <vector>
#include <string.h>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //


/*!
 *	\date			23/apr/2014
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *
 *	\brief Transformation Mapper
 *
 *	Definition of the transformation from the logical domain to the physical domain.
 *	It contains a default (temporary) implementation of a scaling and translation mapper
 *	of logical octree.
 *	Class_Map has to be implemented and customized by the user for different applications.
 */
template<int dim>
class Class_Map{

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
public:
	double X0;						/**<Coordinate X of the origin of the octree in the physical domain*/
	double Y0;						/**<Coordinate Y of the origin of the octree in the physical domain*/
	double Z0;						/**<Coordinate Z of the origin of the octree in the physical domain*/
	double L;						/**<Side length of octree in the physical domain*/
	Class_Global<dim> globals;		/**<Global variables*/

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Class_Map();
	Class_Map(double & X, double & Y, double & Z, double & LL);

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	/*! Transformation of coordinate X.
	 * \param[in] X Coordinate X from logical domain.
	 * \return Coordinate X in physical domain.
	 */
	double mapX(uint32_t const & X);

	/*! Transformation of coordinate Y.
	 * \param[in] Y Coordinate Y from logical domain.
	 * \return Coordinate Y in physical domain.
	 */
	double mapY(uint32_t const & Y);

	/*! Transformation of coordinate Z.
	 * \param[in] Z Coordinate Z from logical domain.
	 * \return Coordinate Z in physical domain.
	 */
	double mapZ(uint32_t const & Z);

	/*! Transformation of coordinate X.
	 * \param[in] X Coordinate X from physical domain.
	 * \return Coordinate X in logical domain.
	 */
	uint32_t mapX(double const & X);

	/*! Transformation of coordinate Y.
	 * \param[in] Y Coordinate Y from physical domain.
	 * \return Coordinate Y in logical domain.
	 */
	uint32_t mapY(double const & Y);

	/*! Transformation of coordinate Z.
	 * \param[in] Z Coordinate Z from physical domain.
	 * \return Coordinate Z in logical domain.
	 */
	uint32_t mapZ(double const & Z);

	/*! Transformation of size of an octant.
	 * \param[in] size Size of octant from logical domain.
	 * \return Size of octant in physical domain.
	 */
	double mapSize(uint32_t const & size);

	/*! Transformation of area of an octant.
	 * \param[in] area Area of octant from logical domain.
	 * \return Area of octant in physical domain.
	 */
	double mapArea(uint32_t const & area);

	/*! Transformation of volume of an octant.
	 * \param[in] volume Volume of octant from logical domain.
	 * \return Coordinate Volume of octant in physical domain.
	 */
	double mapVolume(uint32_t const & volume);

	/*! Transformation of coordinates of center of an octant.
	 * \param[in] center Pointer to coordinates of center from logical domain.
	 * \param[out] mapcenter Coordinates of center in physical domain.
	 */
	void mapCenter(double* & center,
					vector<double> & mapcenter);

	/*! Transformation of coordinates of nodes of an octant.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodes(uint32_t (*nodes)[3],
					vector<vector<double> > & mapnodes);

	/*! Transformation of coordinates of nodes of an intersection.
	 * \param[in] nodes Pointer to coordinates of nodes from logical domain.
	 * \param[out] mapnodes Coordinates of nodes in physical domain.
	 */
	void mapNodesIntersection(uint32_t (*nodes)[3],
					vector<vector<double> > & mapnodes);

	/*! Transformation of components of normal of an intersection.
	 * \param[in] nodes Pointer to components of normal from logical domain.
	 * \param[out] mapnodes components of normal in physical domain.
	 */
	void mapNormals(vector<int8_t> normal,
					vector<double> & mapnormal);
	// ------------------------------------------------------------------------------- //

};

#include "Class_Map.tpp"

// end class_map

#endif /* CLASS_MAP_HPP_ */
