/*
 * octant.hpp
 *
 *  Created on: Feb 10, 2014
 *      Author: Edoardo Lombardi
 */

#ifndef OCTANT_HPP_
#define OCTANT_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "preprocessor_defines.dat"
#include "global.hpp"
#include "Log_Funct.hpp"
#include <vector>
#include <string.h>
#include "inlinedFunct.hpp"


// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;


// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

class Class_Octant{

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //

private:
	uint32_t  x, y, z;			// Coordinates
	uint8_t   level;			// Refinement level (0=root)
	int8_t    marker;			// Set for Refinement(m>0) or Coarsening(m<0) |m|-times
	bool      info[16];			// Info[0..5] : true if 0..5 face is a boundary face [bound];
								// Info[6..11]: true if 0..5 face is a process boundary face [pbound];
								// Info[12/13]: true if octant is new after refinement/coarsening;
								// Info[14]   : true if balancing is not required for this octant;
								// Info[15]   : true if octant is a scary ghost.


	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //

public:
	Class_Octant();
	Class_Octant(int8_t level, int32_t x, int32_t y, int32_t z);
	Class_Octant(const Class_Octant &octant);


	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //

	// Basic Get/Set methods --------------------------------------------------------- //

public:
	uint32_t  getx() const;
	uint32_t  gety() const;
	uint32_t  getz() const;
	uint8_t   getlevel() const;
	int8_t    getmarker() const;
	bool      getbound(uint8_t face) const;		// Get if face is boundary
	bool      getpbound(uint8_t face) const;	// Get if face is process boundary
	bool      getisnewR() const;				// Get if octant is new after refinement
	bool      getisnewC() const;				// Get if octant is new after coarsening
	bool      getbalance() const;				// Get if balancing-blocked octant
	bool      getisghost() const;				// For ghostbusters : get if octant is a ghost

	void      setmarker(int8_t marker);			// Set refinement/coarsening marker
	void      setbalance(bool balance);			// Set if balancing-blocked octant
	uint64_t  computeMorton();					// Compute Morton index of the octant (without level)
private:
	void      setlevel(uint8_t level);

	//-------------------------------------------------------------------------------- //
	// Other Get/Set methods --------------------------------------------------------- //

public:
	uint32_t  getsize() const;					// Get the size of octant
	uint32_t  getvolume() const;				// Get the volume of octant

private:

	//-------------------------------------------------------------------------------- //
	// Other methods ----------------------------------------------------------------- //

public:
	void  buildchildren(vector<Class_Octant> & children);	// Builds children of octant and stores them in vector children[nchildren]

private:

	// ------------------------------------------------------------------------------- //


};//end Class_Octant;




#endif /* OCTREE_HPP_ */
