/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

// ========================================================================== //
//                              DGF IO FUNCTIONS                              //
//                                                                            //
// I/O functions for .dgf file formats                                        //
// ========================================================================== //
// INTFO                                                                      //
// ========================================================================== //
// Author  : Alessandro Alaia                                                 //
// Version : v3.0                                                             //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

# include "DGF.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //
using namespace std;

namespace bitpit{

/*!
    \class DGFObj
    \ingroup DuneGridFormat
    \brief Interface to DGF I/O function

    This class has been designed to allow an easy interface between end-user
    and DGF I/O functions.
*/

// class DGFObj methods ==================================================== //

// Constructors ------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class DGF obj.

    Initialize an empty interface.
*/
DGFObj::DGFObj(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET MEMBERS TO DEFAULT VALUE                                               //
// ========================================================================== //

// General info
dgf_name = "";

// Error flags
err = 0;

// dgf content
data.nV = -1;
data.nS = -1;


return; }

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class DGFObj.

    Initialize an interface associated to a dgf file with specified name.

    \param[in] filename dgf file name
*/
DGFObj::DGFObj(
    std::string                          filename
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET MEMBERS TO DEFAULT VALUE                                               //
// ========================================================================== //

// General info
dgf_name = utils::string::trim(filename);

// Error flags
err = 0;

// dgf content
data.nV = -1;
data.nS = -1;

return; }

// Destructors -------------------------------------------------------------- //
// none

// Public methods ----------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Open stream to dgf file associated with the interface.

    \param[in] mode opening mode ("in": input, "out": output, "app": append mode)
*/
void DGFObj::open(
    std::string                          mode
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN STREAM                                                                //
// ========================================================================== //
if (mode.compare("in") == 0) {
    if (!ifile_handle.is_open()) {

        // Close open output stream
        close("out");

        // Open input stream
        ifile_handle.open(dgf_name, ifstream::in | ifstream::binary);
        if (!ifile_handle.good()) { err = 1; }

    }
}
else if (mode.compare("out") == 0) {
    if (!ofile_handle.is_open()) {

        // Close open input stream
        close("in");

        // Open output stream
        ofile_handle.open(dgf_name, ifstream::out);
        if (!ofile_handle.good()) { err = 1; }        

    }
}
else if (mode.compare("app") == 0) {
    if (!ofile_handle.is_open()) {

        // Close open input stream
        close("in");

        // Open output stream in "append" mode
        ofile_handle.open(dgf_name, ifstream::out | ifstream::app);
        if (!ofile_handle.good()) { err = 1; }

    }
}

return; };

// -------------------------------------------------------------------------- //
/*!
    Close stream from/to the dgf file associated with the interface.

    \param[in] mode opening mode used for the stream ("in": input, "out": output,
    "app": append mode)
*/
void DGFObj::close(
    std::string                          mode
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CLOSE STREAM                                                               //
// ========================================================================== //
if (mode.compare("in") == 0) {
    if (ifile_handle.is_open()) {
        ifile_handle.close();
    }
}
else if (mode.compare("out") == 0) {
    if (ofile_handle.is_open()) {
        ofile_handle.close();
    }
}
else if (mode.compare("app") == 0) {
    if (ofile_handle.is_open()) {
        ofile_handle.close();
    }
}
else if (mode.compare("inout") == 0) {
    if (ifile_handle.is_open()) {
        ifile_handle.close();
    }
    if (ofile_handle.is_open()) {
        ofile_handle.close();
    }
}



return; };

// -------------------------------------------------------------------------- //
/*!
    Clear info gathered from the associated dgf file,
*/
void DGFObj::clear(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// RESET MEMBERS VALUES TO DEFAULT                                                       //
// ========================================================================== //

// Error flags
err = 0;
dgf_error.resize(0);

// dgf content
data.nV = -1;
data.nS = -1;
data.nV_data.resize(0);
data.sV_data.resize(0);
data.nS_data.resize(0);
data.sS_data.resize(0);

// Close open stream
close("inout");

return; };

// -------------------------------------------------------------------------- //
/*!
    Display content and info to output stream

    \param[in,out] out output stream
*/
void DGFObj::display(
    std::ostream                        &out
) {

// ========================================================================== //
// void DGFObj::display(                                                     //
//     ostream     &out)                                                      //
//                                                                            //
// Display dgf content and infos.                                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - out            : ostream, output stream                                  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int         n, i;

// ========================================================================== //
// DISPLAY INFOS                                                              //
// ========================================================================== //

// General infos ------------------------------------------------------------ //
out << "dgf object:" << endl;
out << "  dgf name              : '" << dgf_name << "'" << endl;
out << "  input stream status   : ";
if (ifile_handle.is_open()) { out << "open"; }
else                        { out << "closed"; }
out << endl;
out << "  output stream status  : ";
if (ofile_handle.is_open()) { out << "open"; }
else                        { out << "closed"; }
out << endl;

// Error status ------------------------------------------------------------- //
if (err == 1) {
    out << "  **ERROR** dgf file is missing!!" << endl;
}

// Infos -------------------------------------------------------------------- //
if (data.nV >= 0) {
    out << "  # of vertices         : " << data.nV << endl;
}
if (data.nS >= 0) {
    out << "  # of simplicies       : " << data.nS << endl;
}
if (data.nV_data.size() > 0) {
    out << "  vertex dataset names  : " << data.sV_data << endl;
    out << "  # of vertex data      : " << data.nV_data << endl;
}
if (data.nS_data.size() > 0) {
    out << "  simplex dataset names : " << data.sS_data << endl;
    out << "  # of simplex data     : " << data.nS_data << endl;
}

// Error on data sets ------------------------------------------------------- //
if (dgf_error.size() > 0) {

    out << "ERROR report: " << endl;

    // Vertex
    out << "  vertex: " << endl;
    switch (dgf_error[0][0]) {
        case 0: {out << "    no errors" << endl; break; }
        case 1: {out << "    **WARNING** unterminated data block!!" << endl; break; }
        case 2: {out << "    **WARNING** missing data" << endl; break; }
    };

    // Simplex
    out << "  simplex: " << endl;
    switch (dgf_error[1][0]) {
        case 0: {out << "    no errors" << endl; break; }
        case 1: {out << "    **WARNING** unterminated data block!!" << endl; break; }
        case 2: {out << "    **WARNING** missing data" << endl; break; }
    };

    // Vertex data sets
    n = data.nV_data.size();
    for (i = 0; i < n; i++) {
        out << "  vertex dataset '" << data.sV_data[i] << "': " << endl;
        switch (dgf_error[2][i]) {
            case 0: {out << "    no errors" << endl; break; }
            case 1: {out << "    **WARNING** unterminated data block!!" << endl; break; }
            case 2: {out << "    **WARNING** missing data" << endl; break; }
        };
    } //next i

    // Simplex data sets
    n = data.nS_data.size();
    for (i = 0; i < n; i++) {
        out << "  simplex dataset '" << data.sS_data[i] << "': " << endl;
        switch (dgf_error[3][i]) {
            case 0: {out << "    no errors" << endl; break; }
            case 1: {out << "    **WARNING** unterminated data block!!" << endl; break; }
            case 2: {out << "    **WARNING** missing data" << endl; break; }
        };
    } //next i
    
}

return; }

// -------------------------------------------------------------------------- //
/*!
    Scan dgf file associated to the interface and gather infos.
*/
void DGFObj::scan(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// SCAN DGF FILE                                                              //
// ========================================================================== //
err = dgf::scan(ifile_handle,
               data.nV,
               data.nS,
               data.sV_data,
               data.sS_data,
               data.nV_data,
               data.nS_data);

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; }

// -------------------------------------------------------------------------- //
/*!
    Scan dgf file associated with the interface and perform check on format
    error(s).
*/
void DGFObj::check(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SCAN DGF FILE                                                              //
// ========================================================================== //
scan();

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// CHECK DATA COHERENCY                                                       //
// ========================================================================== //
err = dgf::check(ifile_handle, dgf_error);

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; }

// -------------------------------------------------------------------------- //
/*!
    Load mesh data from DGF file.

    \param[in,out] nV on input stores the number of vertices already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh vertices loaded from the dgf file.
    \param[in,out] nS on input stores the number of cells already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh cells loaded from the dgf file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of
    the vertices already acquired from the dgf file. New coordinates are appended
    at the end of V.
    \param[in,out] S cell->vertex connectivity data. On output stores the connectivity
    entries for cells acquired from the dgf file. New connectivity entries are
    appended at the and of S.
*/
void DGFObj::load(
    int                                 &nV,
    int                                 &nS,
    std::vector<std::vector<double> >   &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// READ MESH DATA FROM DGF FILE                                               //
// ========================================================================== //
err = dgf::readMesh(ifile_handle, nV, nS, V, S);

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Load mesh data from DGF file. Overloading of member function DGFObj::load()
    for container vector<array<double,3>>

    \param[in,out] nV on input stores the number of vertices already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh vertices loaded from the dgf file.
    \param[in,out] nS on input stores the number of cells already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh cells loaded from the dgf file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of
    the vertices already acquired from the dgf file. New coordinates are appended
    at the end of V.
    \param[in,out] S cell->vertex connectivity data. On output stores the connectivity
    entries for cells acquired from the dgf file. New connectivity entries are
    appended at the and of S.
*/
void DGFObj::load(
    int                                 &nV,
    int                                 &nS,
    std::vector<std::array<double,3> >  &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// READ MESH DATA FROM DGF FILE                                               //
// ========================================================================== //
err = dgf::readMesh(ifile_handle, nV, nS, V, S);

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Load mesh data from DGF file. Overloading of member function DGFObj::load()
    for container vector<array<double,3>>

    \param[in,out] nV on input stores the number of vertices already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh vertices loaded from the dgf file.
    \param[in,out] nS on input stores the number of cells already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh cells loaded from the dgf file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of
    the vertices already acquired from the dgf file. New coordinates are appended
    at the end of V.
    \param[in,out] S cell->vertex connectivity data. On output stores the connectivity
    entries for cells acquired from the dgf file. New connectivity entries are
    appended at the and of S.
    \param[out] PID part identifier list for each simplex.
    \param[in] pidName name of the SIMPLEXDATA to  be used as PID.
*/
void DGFObj::load(
    int                                 &nV,
    int                                 &nS,
    std::vector<std::array<double,3> >  &V,
    std::vector<std::vector<int> >      &S,
    std::vector<int>                    &PID,
    std::string                         pidName
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// READ MESH DATA FROM DGF FILE                                               //
// ========================================================================== //
err = dgf::readMesh(ifile_handle, nV, nS, V, S);

// ========================================================================== //
// READ PID                                                                   //
// ========================================================================== //
int nData = PID.size() ;
err = dgf::readSimplexData(ifile_handle, nData, PID, pidName ) ;
if (nData == (int) PID.size()) {
    PID.resize(nS - nData, 0);
}

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Save mesh data from DGF file.

    \param[in] nV number of vertices in the mesh
    \param[in] nS number of cells in the mesh
    \param[in] V vertex coordinate list.
    \param[in] S cell->vertex connectivity data.
*/
void DGFObj::save(
    int                                 &nV,
    int                                 &nS,
    std::vector<std::vector<double> >   &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// void DGFObj::save(                                                        //
//     int         &nV,                                                       //
//     int         &nS,                                                       //
//     dvector2D   &V,                                                        //
//     ivector2D   &S)                                                        //
//                                                                            //
// SAve mesh data into a dgf file.                                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - nV          : int, number of mesh vertices                               //
// - nS          : int, number of simplicies                                  //
// - V           : dvector2D, vertex coordinate list. V[i][0], V[i][1], ...   //
//                 are the x, y, ... coordinates of the i-th vertex           //
// - S           : ivector2D, simplex-vertex connectivity S[i][0], S[i][1],   //
//                 ... are the global indices of vertices of the i-th simplex //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("out");

// ========================================================================== //
// READ MESH DATA FROM DGF FILE                                               //
// ========================================================================== //
err = dgf::writeMesh(ofile_handle, nV, nS, V, S);

// ========================================================================== //
// CLOSE OUTPUT STREAM                                                        //
// ========================================================================== //
close("out");

return; }

// -------------------------------------------------------------------------- //
/*!
    Save mesh data from DGF file. Overloading of member function DGFObj::save()
    for container vector<array<double,3>>

    \param[in] nV number of vertices in the mesh
    \param[in] nS number of cells in the mesh
    \param[in] V vertex coordinate list.
    \param[in] S cell->vertex connectivity data.
*/
void DGFObj::save(
    int                                 &nV,
    int                                 &nS,
    std::vector<std::array<double,3> >  &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("out");

// ========================================================================== //
// READ MESH DATA FROM DGF FILE                                               //
// ========================================================================== //
err = dgf::writeMesh(ofile_handle, nV, nS, V, S);

// ========================================================================== //
// CLOSE OUTPUT STREAM                                                        //
// ========================================================================== //
close("out");

return; }

// Private methods ---------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Dummy function to end recursive calls to DGFObj::loadVData().
*/
void DGFObj::loadVData(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; }

// -------------------------------------------------------------------------- //
/*!
    Dummy function to end recursive calls to DGFObj::loadSData().
*/
void DGFObj::loadSData(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; }

// -------------------------------------------------------------------------- //
/*!
    Dummy function to end recursive calls to DGFObj::appendVData().
*/
void DGFObj::appendVData(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CLOSE OUTPUT STREAM                                                        //
// ========================================================================== //
close("app");

return; }

// -------------------------------------------------------------------------- //
/*!
    Dummy function to end recursive calls to DGFObj::appendSData().
*/
void DGFObj::appendSData(
    void
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// CLOSE OUTPUT STREAM                                                        //
// ========================================================================== //
close("app");

return; }

// Scanning routines ======================================================== //

// -------------------------------------------------------------------------- //
/*!
    Scan dgf data block and returns the number of entries in the block.

    \param[in,out] file_handle input stream from dgf file
    \param[in,out] n on output stores the number of entries found in the data
    block

    \result error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan dgf file.
*/
unsigned int dgf::scanData(
    std::ifstream                       &file_handle,
    int                                 &n
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int        current_pos;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN DATA SET                                                              //
// ========================================================================== //

// Get current line
current_pos = file_handle.tellg();
word = "begin";

// Loop until end of data block
while (!file_handle.eof()
    && ((word.compare("#") != 0)
    &&  (word.compare("VERTEX") != 0)
    &&  (word.compare("SIMPLEX") != 0)
    &&  (word.compare("VERTEXDATA") != 0)
    &&  (word.compare("SIMPLEXDATA") != 0))) {

    // Get next line
    current_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Check if data exist
    if (sline >> word) {
        n++;
    }
    
} //next line
if ((word.compare("#") == 0)
 || (word.compare("VERTEX") == 0)
 || (word.compare("SIMPLEX") == 0)
 || (word.compare("VERTEXDATA") == 0)
 || (word.compare("SIMPLEXDATA") == 0)){
    file_handle.clear();
    file_handle.seekg(current_pos);
    n--;
}

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Scan dgf file and gather infos

    \param[in,out] file_handle input stream from dgf file.
    \param[in,out] nV on output stores the number of mesh vertices
    \param[in,out] nS on output stores the number of mesh cells
    \param[in,out] sV_data name associated to each vertex data block
    \param[in,out] sS_data name associated to each cell data block
    \param[in,out] nV_data number of data entries for each vertex data block
    \param[in,out] nS_data number of data entries for each cell data block

    \result error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan dgf file.
*/
unsigned int dgf::scan(
    std::ifstream                       &file_handle,
    int                                 &nV,
    int                                 &nS,
    std::vector<std::string>            &sV_data,
    std::vector<std::string>            &sS_data,
    std::vector<int>                    &nV_data,
    std::vector<int>                    &nS_data
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int    err = 0;
long int        start_pos;
string          line, word;
stringstream    sline;

// Counters
int             n;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// RESET INPUT DATA                                                           //
// ========================================================================== //
nV = 0;
nS = 0;
sV_data.resize(0);
sS_data.resize(0);
nV_data.resize(0);
nS_data.resize(0);

// ========================================================================== //
// RESET CURSOR POSITION TO FILE BEGIN                                        //
// ========================================================================== //
start_pos = file_handle.tellg();
file_handle.clear();
file_handle.seekg(0);

// ========================================================================== //
// SCAN DGF FILE                                                              //
// ========================================================================== //
word = "begin";
while (!file_handle.eof()) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for relevant keyword
    if (sline >> word) {
        if (word.compare("VERTEX") == 0) {

            // Get number of vertices
            err = dgf::scanData(file_handle, nV);
        }
        else if (word.compare("SIMPLEX") == 0) {

            // Get number of simplicies
            err = dgf::scanData(file_handle, nS);
        }
        else if (word.compare("VERTEXDATA") == 0) {
        
            // Get data set name
            if (sline >> word)  { sV_data.push_back(word);  }
            else                { sV_data.push_back("");    }

            // Get number of data in the dataset
            n = 0;
            err = dgf::scanData(file_handle, n);
            nV_data.push_back(n);
        }
        else if (word.compare("SIMPLEXDATA") == 0) {

            // Get data set name
            if (sline >> word)  { sS_data.push_back(word);  }
            else                { sS_data.push_back("");    }

            // Get number of data in the dataset
            n = 0;
            err = dgf::scanData(file_handle, n);
            nS_data.push_back(n);
        }
    }
} //next line

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(err); }

// Checking routines ======================================================== //

// -------------------------------------------------------------------------- //
/*!
    Scan dgf data block and check for format errors.

    \param[in,out] file_handle input stream from dgf file.
    \param[in,out] err_code error code for format error founds in the data block
        err_code = 0: no error(s) found
        err_code = 1: unterminated data block
        err_code = 2: unreadable data

    \result error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan dgf file.
*/
unsigned int dgf::checkData(
    std::ifstream                       &file_handle,
    int                                 &err_code
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int        current_pos;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); };

// ========================================================================== //
// RESET INPUT VARIABLES                                                      //
// ========================================================================== //
err_code = 0;

// ========================================================================== //
// CHECK DATA BLOCK                                                           //
// ========================================================================== //
word = "begin";
while (!file_handle.eof()
    && ((word.compare("#") != 0)
    &&  (word.compare("VERTEX") != 0)
    &&  (word.compare("SIMPLEX") != 0)
    &&  (word.compare("VERTEXDATA") != 0)
    &&  (word.compare("SIMPLEXDATA") != 0))) {

    // Get current line
    current_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for relevant keywords
    if (!(sline >> word)) {
        err_code = 2;
    }
} //next line
if (word.compare("#") != 0) {
    file_handle.clear();
    file_handle.seekg(current_pos);
    err_code = 1;
}

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Scan dgf file and check for format error(s) within the file.

    \param[in,out] file_handle input stream from dgf file.
    \param[in,out] err_code on output stores the error code for format error(s)
    found in each data section:
        err_code[0][0] stores the list of error codes found in the vertex section
        err_code[1][0] stores the list of error codes found in the cell section
        err_code[2][i] stores the list of error codes found in the i-th vertex data section
        err_code[3][i] stores the list of error codes found in the i-th cell data section
    Error codes are listed below:
        err_code[i][j] = 0: no error(s) found
        err_code[i][j] = 1: unterminated data block
        err_code[i][j] = 2: unreadable data

    \result error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan dgf file.
*/
unsigned int dgf::check(
    std::ifstream                       &file_handle,
    std::vector<std::vector<int> >      &err_code
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int    err = 0;
long int        start_pos;
int             loc_err;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// CHECK INPUT STREAM                                                         //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// RESET INPUT VARIABLES                                                      //
// ========================================================================== //
err_code.resize(4);

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// SCAN DGF FILE                                                              //
// ========================================================================== //
word = "begin";
while ((!file_handle.eof()) && (word.compare("#") != 0)) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for keywords
    if (sline >> word) {
        if (word.compare("VERTEX") == 0) {
            err_code[0].resize(1);
            err = dgf::checkData(file_handle, err_code[0][0]);
        }
        else if (word.compare("SIMPLEX") == 0) {
            err_code[1].resize(1);
            err = dgf::checkData(file_handle, err_code[1][0]);
        }
        else if (word.compare("VERTEXDATA") == 0) {
            err = dgf::checkData(file_handle, loc_err);
            err_code[2].push_back(loc_err);
        }
        else if (word.compare("SIMPLEXDATA") == 0) {
            err = dgf::checkData(file_handle, loc_err);
            err_code[3].push_back(loc_err);
        }
    }
} //next line

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(err); }

// Input functions ========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Load mesh from dgf file.

    \param[in] file_handle input file stream
    \param[in,out] nV on input stores the number of vertices already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh vertices loaded from the dgf file.
    \param[in,out] nS on input stores the number of cells already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh cells loaded from the dgf file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of
    the vertices already acquired from the dgf file. New coordinates are appended
    at the end of V.
    \param[in,out] S cell->vertex connectivity data. On output stores the connectivity
    entries for cells acquired from the dgf file. New connectivity entries are
    appended at the and of S.

    \result error flag for I/O errors.
        err = 0: no error(s) occurred
        err = 1: failed to import mesh data from dgf file.
*/
unsigned int dgf::readMesh(
    std::ifstream                       &file_handle,
    int                                 &nV,
    int                                 &nS,
    std::vector<std::vector<double> >   &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int    err;
long int        start_pos;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// LOAD MESH DATA FROM FILE                                                   //
// ========================================================================== //
word = "begin";
while (!file_handle.eof()
    && (word.compare("#") != 0)) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for keywords
    if (sline >> word) {
        if (word.compare("VERTEX") == 0) {
            err = readData(file_handle, nV, V);
        }
        else if (word.compare("SIMPLEX") == 0) {
            err = readData(file_handle, nS, S);
        }
    }
} //next line

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(err); };

// -------------------------------------------------------------------------- //
/*!
    Load mesh from dgf file. Overloading of function dgf::readMesh()
    for container vector<array<double,3>>

    \param[in] file_handle input file stream
    \param[in,out] nV on input stores the number of vertices already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh vertices loaded from the dgf file.
    \param[in,out] nS on input stores the number of cells already acquired
    from the dgf file. On output stores the input value incremented by 
    the number of mesh cells loaded from the dgf file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of
    the vertices already acquired from the dgf file. New coordinates are appended
    at the end of V.
    \param[in,out] S cell->vertex connectivity data. On output stores the connectivity
    entries for cells acquired from the dgf file. New connectivity entries are
    appended at the and of S.

    \result error flag for I/O errors.
        err = 0: no error(s) occurred
        err = 1: failed to import mesh data from dgf file.
*/
unsigned int dgf::readMesh(
    std::ifstream                       &file_handle,
    int                                 &nV,
    int                                 &nS,
    std::vector<std::array<double,3> >  &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int    err;
long int        start_pos;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// LOAD MESH DATA FROM FILE                                                   //
// ========================================================================== //
word = "begin";
while (!file_handle.eof()
    && (word.compare("#") != 0)) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for keywords
    if (sline >> word) {
        if (word.compare("VERTEX") == 0) {
            err = dgf::readData(file_handle, nV, V);
        }
        else if (word.compare("SIMPLEX") == 0) {
            err = dgf::readData(file_handle, nS, S);
        }
    }
} //next line

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(err); };

// Output functions ========================================================= //

// -------------------------------------------------------------------------- //
/*!
    Export mesh to dgf file.

    \param[in,out] file_handle output stream to dgf file.
    \param[in] nV number of mesh vertices
    \param[in] nS number of mesh cells
    \param[in] V vertex coordinate list
    \param[in] S cell->vertex connectivity

    \result error flag for I/O errors:
    err = 0: no error(s) encountered
    err = 1: failed to write data to dgf file.
*/
unsigned int dgf::writeMesh(
    std::ofstream                       &file_handle,
    int                                 &nV,
    int                                 &nS,
    std::vector<std::vector<double> >   &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int            err = 0;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// EXPORT MESH DATA                                                           //
// ========================================================================== //

// Vertex coordinate list
file_handle << "VERTEX" << endl;
err = dgf::writeData(file_handle, nV, V);

// Simplex-vertex connectivity
file_handle << "SIMPLEX" << endl;
err = dgf::writeData(file_handle, nS, S);

return(err); }

// -------------------------------------------------------------------------- //
/*!
    Export mesh to dgf file. Overloading of function dgf::writeMesh() for
    container vector<array<double,3> >

    \param[in,out] file_handle output stream to dgf file.
    \param[in] nV number of mesh vertices
    \param[in] nS number of mesh cells
    \param[in] V vertex coordinate list
    \param[in] S cell->vertex connectivity

    \result error flag for I/O errors:
    err = 0: no error(s) encountered
    err = 1: failed to write data to dgf file.
*/
unsigned int dgf::writeMesh(
    std::ofstream                       &file_handle,
    int                                 &nV,
    int                                 &nS,
    std::vector<std::array<double,3> >  &V,
    std::vector<std::vector<int> >      &S
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int            err = 0;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// EXPORT MESH DATA                                                           //
// ========================================================================== //

// Vertex coordinate list
file_handle << "VERTEX" << endl;
err = dgf::writeData(file_handle, nV, V);

// Simplex-vertex connectivity
file_handle << "SIMPLEX" << endl;
err = dgf::writeData(file_handle, nS, S);

return(err); }

}
