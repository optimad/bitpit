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
	// FRIENDSHIPS ------------------------------------------------------------------- //
	friend class Class_Local_Tree;

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
	uint32_t	getX() const;
	uint32_t	getY() const;
	uint32_t	getZ() const;
	uint8_t		getLevel() const;
	int8_t		getMarker() const;
	bool		getBound(uint8_t face) const;		// Get if face is boundary
	bool		getPbound(uint8_t face) const;	// Get if face is process boundary
	bool		getIsNewR() const;				// Get if octant is new after refinement
	bool		getIsNewC() const;				// Get if octant is new after coarsening
	bool		getBalance() const;				// Get if balancing-blocked octant
	bool		getIsGhost() const;				// For ghostbusters : get if octant is a ghost

	void		setMarker(int8_t marker);			// Set refinement/coarsening marker
	void		setBalance(bool balance);			// Set if balancing-blocked octant

private:
	void		setLevel(uint8_t level);

	//-------------------------------------------------------------------------------- //
	// Other Get/Set methods --------------------------------------------------------- //

public:
	uint32_t	getSize() const;					// Get the size of octant
	uint32_t	getVolume() const;					// Get the volume of octant
	uint32_t	(*getNodes())[DIM];					// Get a pointer to the array (size [nnodes][3]) with the nodes of octant
	uint64_t	computeMorton() const;				// Compute Morton index of the octant (without level)
	uint64_t	computeMorton();

private:

	//-------------------------------------------------------------------------------- //
	// Other methods ----------------------------------------------------------------- //

public:
	Class_Octant*	buildChildren();							// Builds children of octant and return a pointer to an ordered array children[nchildren]
	uint64_t* 		computeHalfSizeMorton(uint8_t iface, 		// Compute Morton index (without level) of "n=sizehf" half-size (or same size if level=maxlevel)
										  uint8_t & sizehf);	// possible neighbours of octant throught face iface (sizehf=0 if boundary octant)

private:

	// ------------------------------------------------------------------------------- //


};//end Class_Octant;




#endif /* OCTREE_HPP_ */
