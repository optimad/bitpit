#ifndef ARRAY_HPP_
#define ARRAY_HPP_

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include <cstdint>

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*! @ingroup PABLO 
 * @{ 
 *
 *	\date		03/mar/2014
 *  \authors	Marco Cisternino
 *  \authors	Edoardo Lombardi
 *	\copyright	Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par		License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *  \brief Customized array definition
 *
 *  Array contains a pointer to an array of integer values.
 *  Implemented here for fast using in PABLO.
 *
 */
class Array {

	friend class ParaTree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
private:
	uint32_t m_arraySize;		/**< Size of array */
	int* m_array;				/**< Pointer to array of integers */

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Array();
	Array(const Array& other);
	~Array();
private:
	Array(uint32_t size, int value);

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	Array& operator=(const Array& rhs);
};

/* @} */

#endif /* ARRAY_HPP_ */
