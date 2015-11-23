//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_TYPE_HPP__
#define __PATCHMAN_TYPE_HPP__

/*! \file */

#include <vector>

namespace pman {

struct Adaption
{
	enum Type {
		TYPE_UNKNOWN = -1,
		TYPE_CREATION = 0,
		TYPE_DELETION,
		TYPE_REFINEMENT,
		TYPE_COARSENING,
		TYPE_RENUMBERING
	};

	enum Entity {
		ENTITY_UNKNOWN = -1,
		ENTITY_CELL,
		ENTITY_INTERFACE
	};

	struct Info
	{
		Info()
			: type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN)
		{
		}

		Type type;
		Entity entity;
		std::vector<unsigned long> previous;
		std::vector<unsigned long> current;
	};
};

}

#endif
