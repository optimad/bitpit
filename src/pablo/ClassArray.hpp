#ifndef CLASSARRAY_HPP_
#define CLASSARRAY_HPP_

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
 *  ClassArray contains a pointer to an array of integer values.
 *  Implemented here for fast using in PABLO.
 *
 */
class ClassArray {

	friend class ClassParaTree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
	uint32_t m_arraySize;		/**< Size of array */
	int* m_array;				/**< Pointer to array of integers */

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	ClassArray();
	ClassArray(uint32_t size, int value);
	ClassArray(const ClassArray& other);
	~ClassArray();

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	ClassArray& operator=(const ClassArray& rhs);
};

/* @} */

#endif /* CLASS_ARRAY_HPP_ */
