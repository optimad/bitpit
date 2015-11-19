//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#ifndef __PATCHMAN_OUTPUT_MANAGER_HPP__
#define __PATCHMAN_OUTPUT_MANAGER_HPP__

#include <string>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <unordered_map>

namespace pman {

class Cell;
class Node;

class OutputManager : public vtkUnstructuredGrid {

public:
	// Tipo di dati
	enum {
		DATA_ON_CELL = 0,
		DATA_ON_VERTEX
	};

	// Tipo di file
	enum {
		WRITER_XML = 0,
		WRITER_LEGACY
	};

	// Tipo di elemento
	enum {
		ELEM_POINT = 1,
		ELEM_LINE,
		ELEM_TRIA,
		ELEM_QUAD,
		ELEM_TETRA,
		ELEM_HEXA,
		ELEM_WEDGE,
		ELEM_PYRA
	};

	static OutputManager *New();

	void initialize(long nGuessCells, long nGuessVerts);
	void finalize();

	static int getVTKCellType(Element::Type type);
	static void mapConnectFromVTK(Element::Type elemType, int elemMapDegree, int *elemConnect);
	static void mapConnectToVTK(Element::Type elemType, int elemMapDegree, int *elemConnect);

	vtkTypeMacro(OutputManager, vtkUnstructuredGrid);

	vtkIdType InsertNextCell(const Cell &cell, std::unordered_map<long, vtkIdType> &vertexMap);
	vtkIdType InsertNextVertex(const Node &node);

	void resetFields();
	int addField(int fieldType, const char* fieldName);
	int getFieldIndex(int fieldType, const char* fieldName);
	void addFieldValues(int fieldType, int fieldIndex, double* field);
	void addNextFieldValue(int fieldType, int fieldIndex, double value);
	std::string write(std::string filePrefix);
	std::string writeTime(int timeStep, double time, std::string filePrefix, int clearWrapper);

protected:
	OutputManager();
	~OutputManager();

	void writeTimeWrapper(int timeStep, double time, int reset,
		std::string filePrefix, std::string dataFile);

private:
	// Numero massimo di lati
	const int MAX_SIDES = 4;

	vtkSmartPointer<vtkPoints> m_vertexes;

};

}

#endif
