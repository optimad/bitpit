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

class Patch;

class FlatMapping
{

public:
	FlatMapping();
	FlatMapping(Patch *patch);

	virtual ~FlatMapping();

	virtual void update(const std::vector<Adaption::Info> adaptionData) = 0;

	const std::vector<long> & get_numbering() const;
	const std::vector<long> & get_mapping() const;

protected:
	Patch *m_patch;
	std::vector<long> m_numbering;
	std::vector<long> m_mapping;

};


class CellFlatMapping : public FlatMapping
{

public:
	CellFlatMapping();
	CellFlatMapping(Patch *patch);

	~CellFlatMapping();

	void update(const std::vector<Adaption::Info> adaptionData);

};

}

#endif
