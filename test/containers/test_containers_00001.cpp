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

#include "bitpit_containers.hpp"

using namespace bitpit;

void printElements(PiercedVector<double> &container)
{
	std::cout << std::endl << "  List of elements:" << std::endl;

	if (container.empty()) {
		std::cout << "    > Container is empty!" << std::endl;
		return;
	}

	for (auto itr = container.cbegin(); itr != container.cend(); ++itr) {
		auto id     = itr.getId();
		auto flatId = container.evalFlatIndex(id);
		auto value  = *itr;

		std::cout << std::endl;
		std::cout << "    > id     : " << id << std::endl;
		std::cout << "      value  : " << value << std::endl;
		std::cout << "      flatId : " << flatId << std::endl;
	}
}

void fillContainer(int nElements, PiercedVector<double> &container)
{
	for (int i = 0; i < nElements; i++) {
		PiercedVector<double>::iterator iterator = container.emplace(i);
		*iterator = i;
	}
}

int main()
{
	// Creating an emtpy PiercedVector
	std::cout << std::endl << "::: Creating an empty PiercedVector :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Pierced vector created" << std::endl;

	PiercedVector<double> container;

	// Filling the PiercedVector
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	int nElements = 10;
	fillContainer(nElements, container);

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Deleting all elements of the container
	int id_erase;

	std::cout << std::endl << "::: Testing deletion of all elements of the container :::" << std::endl;
	std::cout << std::endl;

	nElements = container.size();
	for (int id_erase = 0; id_erase < nElements; id_erase += 2) {
		if (!container.exists(id_erase)) {
			continue;
		}

		std::cout << "  Deleting element with id = " << id_erase << std::endl;
		container.erase(id_erase);
	}

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < 10; ++i) {
		std::cout << "  Element id before which there are " << i << " elements: " <<  container.getSizeMarker(i) << std::endl;
	}

	// Resizing the container
	std::cout << std::endl << "::: Testing resize :::" << std::endl;
	std::cout << std::endl;

	int size = 3;
	std::cout << "  Resizing the container to " << size << " elements" << std::endl;

	container.resize(size);

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Inserting elements in the container
	int id_insert;

	std::cout << std::endl << "::: Testing insertion :::" << std::endl;
	std::cout << std::endl;

	id_insert = 10;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		container.emplace(id_insert, id_insert);
	}

	id_insert = 13;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		container.emplaceBack(id_insert, id_insert);
	}

	id_insert = 15;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		container.emplace(id_insert, id_insert);
	}

	id_insert = 17;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		container.emplaceBack(id_insert, id_insert);
	}

	id_insert = 45;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at first avilable position) element with id = " << id_insert << std::endl;
		container.emplace(id_insert, id_insert);
	}

	id_insert = 102;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (at the end) element with id = " << id_insert << std::endl;
		container.emplaceBack(id_insert, id_insert);
	}

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	std::cout << std::endl << "::: Testing insertion after/before :::" << std::endl;
	std::cout << std::endl;

	long id_reference;

	id_insert    = 111;
	id_reference = 13;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceBefore(id_reference, id_insert, id_insert);
	}

	id_insert    = 123;
	id_reference = 102;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceBefore(id_reference, id_insert, id_insert);
	}


	id_insert    = 124;
	id_reference = 10;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (before the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceBefore(id_reference, id_insert, id_insert);
	}

	id_insert    = 125;
	id_reference = 124;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceAfter(id_reference, id_insert, id_insert);
	}

	id_insert    = 126;
	id_reference = 102;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceAfter(id_reference, id_insert, id_insert);
	}

	id_insert    = 127;
	id_reference = 45;
	if (!container.exists(id_insert)) {
		std::cout << "  Inserting (after the element with id = " << id_reference << ") element with id = " << id_insert << std::endl;
		container.emplaceAfter(id_reference, id_insert, id_insert);
	}

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// List of ids
	std::cout << std::endl << "::: List of element ids :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  List of the ids (not ordered)"  << std::endl;
	for (auto const &id : container.getIds(false)) {
		std::cout << "    > id = " << id << std::endl;
	}

	std::cout << std::endl;
	std::cout << "  List of the ids (ordered)"  << std::endl;
	for (auto const &id : container.getIds(true)) {
		std::cout << "    > id = " << id << std::endl;
	}

	// Sort the vector
	std::cout << std::endl << "::: Sort elements by ids :::" << std::endl;

	container.sort();

	printElements(container);

	// Delete all the elements
	std::cout << std::endl << "::: Testing deletion of all elements of the container :::" << std::endl;

	while (!container.empty()) {
		container.popBack();
	}

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillContainer(nElements, container);

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Deleting elements
	std::cout << std::endl << "::: Testing element deletion :::" << std::endl;
	std::cout << std::endl;

	id_erase = 0;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 1;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 2;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 3;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 4;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 8;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	id_erase = 7;
	std::cout << "  Deleting element with id = " << id_erase << std::endl;
	container.erase(id_erase);

	// Inserting elements in the container
	std::cout << std::endl << "::: Testing insertion of element :::" << std::endl;
	std::cout << std::endl;

	id_insert = 40;
	std::cout << "  Emplacing element with id = " << id_insert << std::endl;
	container.emplace(id_insert, id_insert);

	id_insert = 41;
	std::cout << "  Emplacing element with id = " << id_insert << std::endl;
	container.emplace(id_insert, id_insert);

	id_insert = 42;
	std::cout << "  Emplacing element with id = " << id_insert << std::endl;
	container.emplace(id_insert, id_insert);

	printElements(container);

	// Moving elements
	std::cout << std::endl << "::: Testing move :::" << std::endl;
	std::cout << std::endl;

	long id_move;

	id_move      = 42;
	id_reference = 40;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " before element with id = " << id_reference << std::endl;
	container.moveBefore(id_reference, id_move);

	printElements(container);

	id_move      = 40;
	id_reference = 41;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " after element with id = " << id_reference << std::endl;
	container.moveAfter(id_reference, id_move);

	printElements(container);

	id_move      = 5;
	id_reference = 9;
	std::cout << std::endl;
	std::cout << "  Moving element with id = " << id_move << " after element with id = " << id_reference << std::endl;
	container.moveAfter(id_reference, id_move);

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Squeezee the vector
	std::cout << std::endl << "::: Testing squeeze :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Capacity before squeeze = " << container.capacity() << std::endl;

	container.squeeze();

	std::cout << "  Capacity after squeeze  = " << container.capacity() << std::endl;

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Clear the vector
	std::cout << std::endl << "::: Clear the vector :::" << std::endl;

	container.clear();

	printElements(container);

	// Filling the vector again
	std::cout << std::endl << "::: Filling the vector :::" << std::endl;

	fillContainer(nElements, container);

	printElements(container);

	// Find marker
	std::cout << std::endl << "::: Testing marker :::" << std::endl;
	std::cout << std::endl;

	std::cout << "  Element id before which there are 4 elements: " <<  container.getSizeMarker(4) << std::endl;

	// Done
	std::cout << std::endl << "::: Done :::" << std::endl;
	std::cout << std::endl;

}
