//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

/*!
	\example manipulatePiercedVector.cpp

	\brief Creation and manipulation of a pierced vector

	This example shows how a pierced vector can be created and
	manipulated.

	<b>To run</b>: ./manipulatePiercedVector \n
*/

#include "BitP_Mesh_COMMON.hpp"
#include "BitP_Mesh_PATCHMAN.hpp"

void printCellIds(bitpit::PiercedVector<Cell> &cells)
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

void fillCellList(int nCells, bitpit::PiercedVector<Cell> &cells)
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

	bitpit::PiercedVector<Cell> cells;

	// Filling the list of cells
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	int nCells = 10;
	fillCellList(nCells, cells);

	printCellIds(cells);

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

	// Resizing the cell container
	std::cout << std::endl << "::: Resizing the cell container :::" << std::endl;
	std::cout << std::endl;

	int size = 3;
	std::cout << "  Resizing the cell container to " << size << " elements" << std::endl;

	cells.resize(size);

	printCellIds(cells);

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

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillCellList(nCells, cells);

	printCellIds(cells);

	// Deleting cells
	std::cout << std::endl << "::: Deleting cells :::" << std::endl;
	std::cout << std::endl;

	id_erase = 0;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 1;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 3;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 4;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 7;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 8;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	id_erase = 9;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	cells.erase(id_erase);

	printCellIds(cells);

	// Squeezee the vector
	std::cout << std::endl << "::: Squeeze :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Capacity before squeeze = " << cells.capacity() << std::endl;

	cells.squeeze();

	std::cout << "  Capacity after squeeze  = " << cells.capacity() << std::endl;

	printCellIds(cells);

	// Clear the vector
	std::cout << std::endl << "::: Clear the vector :::" << std::endl;

	cells.clear();

	printCellIds(cells);

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillCellList(nCells, cells);

	printCellIds(cells);

	// Done
	std::cout << std::endl << "::: Done :::" << std::endl;
	std::cout << std::endl;
	
}
