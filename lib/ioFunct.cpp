/*
 * ioFunct.cpp
 *
 *  Created on: 21/feb/2014
 *      Author: Marco Cisternino
 */

#include <ioFunct.hpp>
#include <fstream>
#include <iomanip>
#include <string>
#include "Class_Para_Tree.hpp"
#include "logFunct.hpp"

using namespace std;

void writeLocalTree(const u32vector2D& nodes, const u32vector2D& connectivity,
		const u32vector2D& ghostNodes, const u32vector2D& ghostConnectivity,
		const Class_Para_Tree & ptree, string filename) {

	stringstream name;
	name << "s" << std::setfill('0') << std::setw(4) << ptree.nproc << "-p" << std::setfill('0') << std::setw(4) << ptree.rank << "-" << filename << ".vtu";

	ofstream out(name.str().c_str());
	if(!out.is_open()){
		stringstream ss;
		ss << filename << "*.vtu cannot be opened and it won't be written.";
		writeLog(ss.str());
		return;
	}
	int nofNodes = nodes.size();
	int nofGhostNodes = ghostNodes.size();
	int nofOctants = connectivity.size();
	int nofGhosts = ghostConnectivity.size();
	int nofAll = nofGhosts + nofOctants;
	out << "<?xml version=\"1.0\"?>" << endl
		<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
		<< "  <UnstructuredGrid>" << endl
		<< "    <Piece NumberOfCells=\"" << connectivity.size() + ghostConnectivity.size() << "\" NumberOfPoints=\"" << nodes.size() + ghostNodes.size() << "\">" << endl;
    out << "      <Points>" << endl
    	<< "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< DIM <<"\" format=\"ascii\">" << endl
    	<< "          " << std::fixed;
    for(int i = 0; i < nofNodes; i++)
    {
    	for(int j = 0; j < DIM; ++j)
    		out << std::setprecision(6) << nodes[i][j] << " ";
    	if((i+1)%4==0 && i!=nofNodes-1)
    		out << endl << "          ";
    }
    for(int i = 0; i < nofGhostNodes; i++)
    {
    	for(int j = 0; j < DIM; ++j)
    		out << std::setprecision(6) << ghostNodes[i][j] << " ";
    	if((i+1)%4==0 && i!=nofNodes-1)
    		out << endl << "          ";
    }
    out << endl << "        </DataArray>" << endl
    	<< "      </Points>" << endl
        << "      <Cells>" << endl
        << "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
        << "          ";
    for(int i = 0; i < nofOctants; i++)
      {
        for(int j = 0; j < nnodes; j++)
          {
            out << connectivity[i][j] << " ";
          }
        if((i+1)%3==0 && i!=nofOctants-1)
          out << endl << "          ";
      }
    for(int i = 0; i < nofGhosts; i++)
      {
        for(int j = 0; j < nnodes; j++)
          {
            out << ghostConnectivity[i][j] + nofNodes << " ";
          }
        if((i+1)%3==0 && i!=nofGhosts-1)
          out << endl << "          ";
      }
    out << endl << "        </DataArray>" << endl
        << "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
        << "          ";
    for(int i = 0; i < nofAll; i++)
      {
        out << (i+1)*nnodes << " ";
        if((i+1)%12==0 && i!=nofAll-1)
          out << endl << "          ";
      }
    out << endl << "        </DataArray>" << endl
        << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
        << "          ";
    for(int i = 0; i < nofAll; i++)
      {
        int type;
        if(DIM==3)
          type = 11;
        else if(DIM==2)
          type = 9;
        out << type << " ";
        if((i+1)%12==0 && i!=nofAll-1)
          out << endl << "          ";
      }
    out << endl << "        </DataArray>" << endl
        << "      </Cells>" << endl
        << "    </Piece>" << endl
        << "  </UnstructuredGrid>" << endl
        << "</VTKFile>" << endl;


    if(ptree.rank == 0){
        name.str("");
        name << "s" << std::setfill('0') << std::setw(4) << ptree.nproc << "-" << filename << ".pvtu";
        ofstream pout(name.str().c_str());
        if(!pout.is_open()){
        	stringstream ss;
    		ss << filename << "*.pvtu cannot be opened and it won't be written.";
    		writeLog(ss.str());
    		return;
        }

        pout << "<?xml version=\"1.0\"?>" << endl
        		<< "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
        		<< "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
        		<< "    <PPointData>" << endl
        		<< "    </PPointData>" << endl
        		<< "    <PCellData Scalars=\"\">" << endl;
        pout << "    </PCellData>" << endl
        		<< "    <PPoints>" << endl
        		<< "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
        		<< "    </PPoints>" << endl;
        for(int i = 0; i < ptree.nproc; i++)
        	pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << ptree.nproc << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
        pout << "  </PUnstructuredGrid>" << endl
        		<< "</VTKFile>";

        pout.close();

    }
    MPI_Barrier(MPI_COMM_WORLD);


}
