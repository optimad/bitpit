/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/*!
	\example patch_example_00001.cpp

	\brief Creation and manipulation a list of cells

	This example shows how to use the pierced vector to create and
	manipulate a list of cells.

	<b>To run</b>: ./patch_example_00001 \n
*/

#include "bitpit_patch.hpp"

using namespace bitpit;

void printCellIds(PiercedVector<Cell> &cells)
{
	std::cout << std::endl << "  List of cell ids:" << std::endl;

	if (cells.empty()) {
		std::cout << std::endl << "  Vector is empty!" << std::endl;
		return;
	}

	for (auto const &cell : cells) {
		std::cout << "    > id = " << cell.get_id() << std::endl;
	}
}

void fillCellList(int nCells, PiercedVector<Cell> &cells)
{
	for (int i = 0; i < nCells; i++) {
		// std::cout << "  Inserting (at the end) cell with id = " << i << std::endl;

		cells.emplace_back(i);
	}
}

int main(int argc, char *argv[])
{
	// Creating an emtpy list of cells
	std::cout << std::endl << "::: Creating an empty pierced vector :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Pierced vector created" << std::endl;

	PiercedVector<Cell> cells;

	// Filling the list of cells
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	int nCells = 10;
	fillCellList(nCells, cells);

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Deleting cells
	int id_erase;

	std::cout << std::endl << "::: Deleting cells :::" << std::endl;
	std::cout << std::endl;

	nCells = cells.size();
	for (int id_erase = 0; id_erase < nCells; id_erase += 2) {
		if (!cells.exists(id_erase)) {
			continue;
		}

		std::cout << "  Deleting element with id = " << id_erase << std::endl;
		cells.erase(id_erase);
	}

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < 10; ++i) {
		std::cout << "  Element id before which there are " << i << " elements: " <<  cells.get_size_marker(i) << std::endl;
	}

	// Resizing the cell container
	std::cout << std::endl << "::: Resizing the cell container :::" << std::endl;
	std::cout << std::endl;

	int size = 3;
	std::cout << "  Resizing the cell container to " << size << " elements" << std::endl;

	cells.resize(size);

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Inserting cells again
	int id_insert;

	std::cout << std::endl << "::: Inserting cells :::" << std::endl;
	std::cout << std::endl;

	id_insert = 10;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		cells.emplace(id_insert);
	}

	id_insert = 13;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		cells.emplace_back(id_insert);
	}

	id_insert = 15;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		cells.emplace(id_insert);
	}

	id_insert = 17;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		cells.emplace_back(id_insert);
	}

	id_insert = 45;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		cells.emplace(id_insert);
	}

	id_insert = 102;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		cells.emplace_back(id_insert);
	}

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	std::cout << std::endl << "::: Inserting cells after/before :::" << std::endl;
	std::cout << std::endl;

	long id_reference;

	id_insert    = 111;
	id_reference = 13;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_before(id_reference, id_insert);
	}

	id_insert    = 123;
	id_reference = 102;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_before(id_reference, id_insert);
	}


	id_insert    = 124;
	id_reference = 10;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_before(id_reference, id_insert);
	}

	id_insert    = 125;
	id_reference = 124;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_after(id_reference, id_insert);
	}

	id_insert    = 126;
	id_reference = 102;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_after(id_reference, id_insert);
	}

	id_insert    = 127;
	id_reference = 45;
	if (!cells.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		cells.emplace_after(id_reference, id_insert);
	}

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// List of ids
	std::cout << std::endl << "::: List of cells id :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  List of the ids (not ordered)"  << std::endl;
	for (auto const &id : cells.get_ids(false)) {
		std::cout << "    > id = " << id << std::endl;
	}

	std::cout << std::endl;
	std::cout << "  List of the ids (ordered)"  << std::endl;
	for (auto const &id : cells.get_ids(true)) {
		std::cout << "    > id = " << id << std::endl;
	}

	// Sort the vector
	std::cout << std::endl << "::: Sort cells by ids :::" << std::endl;

	cells.sort();

	printCellIds(cells);

	// Delete all the cells
	std::cout << std::endl << "::: Deleting all cells :::" << std::endl;

	while (!cells.empty()) {
		cells.pop_back();
	}

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillCellList(nCells, cells);

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Deleting cells
	std::cout << std::endl << "::: Deleting cells :::" << std::endl;
	std::cout << std::endl;

	id_erase = 0;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 1;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 2;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 3;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 4;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 8;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 7;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	// Inserting cells
	std::cout << std::endl << "::: Inserting cells :::" << std::endl;
	std::cout << std::endl;

	id_insert = 40;
	std::cout << "  Inserting element with id = " << id_insert << std::endl;
	cells.insert(id_insert);

	id_insert = 41;
	std::cout << "  Inserting element with id = " << id_insert << std::endl;
	cells.insert(id_insert);

	id_insert = 42;
	std::cout << "  Inserting element with id = " << id_insert << std::endl;
	cells.insert(id_insert);

	printCellIds(cells);

	// Moving elements
	std::cout << std::endl << "::: Moving cells :::" << std::endl;
	std::cout << std::endl;

	long id_move;

	id_move      = 42;
	id_reference = 40;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " before element with id = " << id_reference << std::endl;
	cells.move_before(id_reference, id_move);

	printCellIds(cells);

	id_move      = 40;
	id_reference = 41;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " after element with id = " << id_reference << std::endl;
	cells.move_after(id_reference, id_move);

	printCellIds(cells);

	id_move      = 5;
	id_reference = 9;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " after element with id = " << id_reference << std::endl;
	cells.move_after(id_reference, id_move);

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Squeezee the vector
	std::cout << std::endl << "::: Squeeze :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Capacity before squeeze = " << cells.capacity() << std::endl;

	cells.squeeze();

	std::cout << "  Capacity after squeeze  = " << cells.capacity() << std::endl;

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Clear the vector
	std::cout << std::endl << "::: Clear the vector :::" << std::endl;

	cells.clear();

	printCellIds(cells);

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillCellList(nCells, cells);

	printCellIds(cells);

	// Find marker
	std::cout << std::endl << "::: Marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  cells.get_size_marker(4) << std::endl;

	// Done
	std::cout << std::endl << "::: Done :::" << std::endl;
	std::cout << std::endl;

}
