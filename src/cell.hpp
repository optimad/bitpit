//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_CELL_HPP__
#define __PATCHMAN_CELL_HPP__

/*! \file */

#include "element.hpp"

#include <memory>

namespace pman {

class Interface;
class Patch;

class CellData {};

class Cell : public Element {

public:
	enum PositionType {
	    INTERNAL = 0,
	    GHOST
	};

	Cell();
	Cell(const int &id);
	Cell(const int &id, Patch *patch);

	Cell(Cell&& other) = default;
	Cell& operator=(Cell&& other) = default;

	void set_position_type(PositionType position);
	PositionType get_position_type() const;
	
	void set_volume(double * const volume);
	double & get_volume() const;

	void initialize_interfaces(std::vector<std::vector<int>> &interfaces);
	void initialize_empty_interfaces(const int nInterfaces[]);
	void set_interface(const int &face, const int &index, const int &interface);
	void set_interfaces(const int &face, int interface[]);
	void unset_interfaces();
	int get_interface_count() const;
	int get_interface_count(const int &face) const;
	int get_interface(const int &face, const int &index = 0) const;
	const int * get_interfaces() const;
	const int * get_interfaces(const int &face) const;

	void set_data(std::unique_ptr<CellData> m_data);
	CellData * get_data() const;

protected:

private:
  	PositionType m_positionType;

	double *m_volume;

	std::unique_ptr<CollapsedArray2D<int> > m_interfaces;

	std::unique_ptr<CellData> m_data;

	Cell(const Cell &other) = delete;
	Cell& operator = (const Cell &other) = delete;

};

}

#endif
