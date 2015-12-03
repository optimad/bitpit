//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include "cell.hpp"
#include "node.hpp"
#include "output_manager.hpp"
#include "patch.hpp"

#include "vtkGenericCell.h"
#include "vtkObjectFactory.h"
#include "vtkSetGet.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkZLibDataCompressor.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace pman {

// Description:
// Crea un nuovo oggetto
vtkStandardNewMacro(OutputManager);


// Description:
// Costruttore del dataset non strutturato
OutputManager::OutputManager()
     : vtkUnstructuredGrid()
{

}


// Description:
// Distruttore del dataset non strutturato
OutputManager::~OutputManager()
{
}


// Description:
// Inizializza il dataset
void OutputManager::initialize(long nGuessCells, long nGuessVerts)
{
	m_vertexes = vtkSmartPointer<vtkPoints>::New();
	m_vertexes->SetNumberOfPoints(nGuessVerts);
	m_vertexes->Initialize();

	Initialize();
	Allocate(nGuessCells);
}


// Description:
// Finalizza il dataset
void OutputManager::finalize()
{
	m_vertexes->Squeeze();
	SetPoints(m_vertexes);

	Squeeze();
}


// Description:
// Restituisce il tipo di elemento VTK
int OutputManager::getVTKCellType(ElementInfo::Type type)
{
	switch (type)  {
		case ElementInfo::VERTEX:
			return VTK_VERTEX;
		case ElementInfo::LINE:
			return VTK_LINE;
		case ElementInfo::TRIANGLE:
			return VTK_TRIANGLE;
		case ElementInfo::PIXEL:
			return VTK_PIXEL;
		case ElementInfo::QUAD:
			return VTK_QUAD;
		case ElementInfo::TETRA:
			return VTK_TETRA;
		case ElementInfo::VOXEL:
			return VTK_VOXEL;
		case ElementInfo::HEXAHEDRON:
			return VTK_HEXAHEDRON;
		case ElementInfo::WEDGE:
			return VTK_WEDGE;
		case ElementInfo::PYRAMID:
			return VTK_PYRAMID;
		default:
			return VTK_EMPTY_CELL;
	}

}

// Description:
// Extract the connectivity in the VTK format and order.
std::vector<vtkIdType> OutputManager::extractVTKCellConnect(const Cell &cell, std::unordered_map<long, vtkIdType> &vertexMap)
{
	int nCellVertices = cell.get_vertex_count();

	const long *cellconnect = cell.get_connect();
	std::vector<vtkIdType> vtkCellConnect;
	for (long k = 0; k < nCellVertices; k++) {
		vtkIdType VTKVertexId = vertexMap.at(cellconnect[k]);
		vtkCellConnect.push_back(VTKVertexId);
	}

	return vtkCellConnect;
}

// Description:
// Aggiunge il successivo vertice della mesh
vtkIdType OutputManager::InsertNextVertex(const Node &vertex)
{
	const std::array<double, 3> &coords = vertex.get_coords();
	double x = coords[0];
	double y = coords[1];
	double z = coords[2];

	return m_vertexes->InsertNextPoint(x, y, z);
}

// Description:
// Aggiunge una cella alla mesh
vtkIdType OutputManager::InsertNextCell(const Cell &cell, std::unordered_map<long, vtkIdType> &vertexMap)
{
	vtkIdType nCellVerts = cell.get_vertex_count();
	int vtkCellType = getVTKCellType(cell.get_type());
	std::vector<vtkIdType> vtkCellConnect = extractVTKCellConnect(cell, vertexMap);

	return vtkUnstructuredGrid::InsertNextCell(vtkCellType, nCellVerts, vtkCellConnect.data());
}


// Description:
// Azzera i dati memorizzati.
void OutputManager::resetFields()
{
	GetPointData()->Initialize();
	GetCellData()->Initialize();
}


// Description:
// Aggiunge un campo.
int OutputManager::addField(int fieldType, const char* fieldName)
{
	vtkSmartPointer<vtkDoubleArray> fieldArray = vtkSmartPointer<vtkDoubleArray>::New();
	fieldArray->SetName(fieldName);

	int index;
	if (fieldType == OutputManager::DATA_ON_VERTEX) {
		index = GetPointData()->AddArray(fieldArray);
	} else if (fieldType == OutputManager::DATA_ON_CELL) {
		index = GetCellData()->AddArray(fieldArray);
	}

	return index;
}


// Description:
// Determina l'index associato ad un array.
int OutputManager::getFieldIndex(int fieldType, const char* fieldName)
{
	int fieldIndex;
	if (fieldType == OutputManager::DATA_ON_VERTEX) {
		static_cast<vtkDoubleArray *>(GetPointData()->GetArray(fieldName, fieldIndex));
	} else if (fieldType == OutputManager::DATA_ON_CELL) {
		static_cast<vtkDoubleArray *>(GetCellData()->GetArray(fieldName, fieldIndex));
	}

	return fieldIndex;
}


// Description:
// Aggiunge i valori di un campo.
void OutputManager::addFieldValues(int fieldType, int fieldIndex, double* field)
{
	int dimField;
	vtkDoubleArray *fieldArray;
	if (fieldType == DATA_ON_VERTEX) {
		dimField   = GetNumberOfPoints();
		fieldArray = static_cast<vtkDoubleArray *>(GetPointData()->GetArray(fieldIndex));
	} else if (fieldType == DATA_ON_CELL) {
		dimField   = GetNumberOfCells();
		fieldArray = static_cast<vtkDoubleArray *>(GetCellData()->GetArray(fieldIndex));
	}

	for (int i = 0; i < dimField; i++)
		fieldArray->InsertNextValue(field[i]);

	fieldArray->Squeeze();
}


// Description:
// Aggiunge il prossimo valore di un campo.
void OutputManager::addNextFieldValue(int fieldType, int fieldIndex, double value)
{
	vtkDoubleArray *fieldArray;
	if (fieldType == DATA_ON_VERTEX) {
		fieldArray = static_cast<vtkDoubleArray *>(GetPointData()->GetArray(fieldIndex));
	} else if (fieldType == DATA_ON_CELL) {
		fieldArray = static_cast<vtkDoubleArray *>(GetCellData()->GetArray(fieldIndex));
	}

	fieldArray->InsertNextValue(value);
}

// Description:
// Scrive su file la griglia non strutturata e i relativi dati.
std::string OutputManager::write(std::string filePrefix)
{
	// Le opzioni possibili sono: WRITER_XML, WRITER_LEGACY.
	int writerType = WRITER_XML;

	// Opzioni per i file XML
	//
	// dataMode: le scelte sono  vtkXMLWriter::Ascii, vtkXMLWriter::Binary
	//           o vtkXMLWriter::Appended.
	// encodeAppendData: valido solo per la modialità Appended,
	//                   determina se codificare i dati in base64
	// compressData: determina se comprimere i dati
	int dataMode         = vtkXMLWriter::Binary;
	int encodeAppendData = 0;
	int compressData     = 1;

	// Opzioni per i file LEGACY
	//
	// fieldType: VTK_ASCII,VTK_BINARY
	int fileType = VTK_ASCII;

	// Scrittura del file
	std::string fileName(filePrefix);

	if (writerType == WRITER_LEGACY) {
		fileName = fileName + ".vtk";
		vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		writer->SetHeader(filePrefix.c_str());
		writer->SetFileType(fileType);
		writer->SetInputData(this);
		writer->SetFileName(fileName.c_str());
		writer->Write();
	} else if (writerType == WRITER_XML) {
		fileName = fileName + ".vtu";
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

		writer->SetDataMode(dataMode);
		if (dataMode == vtkXMLWriter::Appended)
			writer->SetEncodeAppendedData(encodeAppendData);

		if (compressData == 1) {
			vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
			writer->SetCompressor(compressor);
		} else
			writer->SetCompressor(NULL);

		writer->SetInputData(this);
		writer->SetFileName(fileName.c_str());
		writer->Write();
	}

	return fileName;
}


// Description:
// Scrive su file la griglia non strutturata e i relativi dati per
// lo step temporale "timeStep".
std::string OutputManager::writeTime(int timeStep, double time, std::string filePrefix, int clearWrapper)
{

	// Scrittura file con i dati
	char strTimeStep[1024];
	sprintf(strTimeStep, "%06i", timeStep);

	std::string stepPrefix = filePrefix + "_" + std::string(strTimeStep);

	std::string fileName = write(stepPrefix);

	// Scrittura wrapper
	writeTimeWrapper(timeStep, time, clearWrapper, filePrefix.c_str(), fileName.c_str());

	return fileName;
}

// Description:
// Scrive il file wrapper che raccoglie i dataset.
void OutputManager::writeTimeWrapper(int timeStep, double time, int reset,
		          std::string filePrefix, std::string dataFile)
{

	std::ostringstream strStream;
	strStream << time;
	std::string strTime = strStream.str();

	std::string dataFileNoPath = dataFile;
	int endPath = dataFileNoPath.find_last_of('/');
	dataFileNoPath = dataFileNoPath.substr(endPath + 1);

	std::string strDataSet = "<DataSet";
	strDataSet += " timestep=\"" + strTime  + "\""	;
	strDataSet += " file=\"" + dataFileNoPath + "\"";
	strDataSet += "/>";

	std::string wrapperName = filePrefix + ".pvd";

	std::vector<std::string> vtkCollection;

	// Determino se devo creare un nuovo wrapper o modificar quello
	// già esistente.
	bool createNewWrapper = (timeStep == 0) || (reset == 1);
	if (!createNewWrapper) {
		ifstream oldWrapperFile;
		oldWrapperFile.open(wrapperName, ifstream::in);
		oldWrapperFile.close();
		if (oldWrapperFile.fail())
			createNewWrapper = true;
	}

	// Se devo modificare il wrapper esistente carico i timestep
	// contenuti nel wrapper. Leggo solo quelli precedenti al
	// timestep che voglio aggiungere
	if (!createNewWrapper) {
		std::string line;
		ifstream oldWrapperFile;
		oldWrapperFile.open(wrapperName);
		while(!oldWrapperFile.eof()) {
			getline(oldWrapperFile, line);

			int timestepInit = line.find("timestep=\"", 0);
			if (timestepInit <= 0)
				continue;
			timestepInit += 11;

			int timestepEnd  = line.find("\"", timestepInit);

			std::string fileTimeStr = line.substr(timestepInit, timestepEnd - timestepInit + 1);
			std::istringstream i(fileTimeStr);
   			double fileTime;
   			i >> fileTime;

			// Aggiungo solo i timestep precedenti
			if (time <= fileTime)
				break;

			std::size_t startpos = line.find_first_not_of(" \t");
    			if (std::string::npos != startpos)
        			line = line.substr(startpos);

			vtkCollection.push_back(line);
		}
		oldWrapperFile.close();
	}

	vtkCollection.push_back(strDataSet);

	// Scrivo il file
	ofstream wrapperFile(wrapperName, std::ios_base::trunc);

	wrapperFile << "<?xml version=\"1.0\"?>\n";
	wrapperFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	wrapperFile << "\t<Collection>\n";

	for (uint i = 0; i < vtkCollection.size(); i++)
		wrapperFile << "\t\t" << vtkCollection[i] << "\n";

	wrapperFile << "\t</Collection>\n";
	wrapperFile << "</VTKFile>\n";

	wrapperFile.close();
}

}
