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
//                           - STL IO FUNCTIONS -                             //
//                                                                            //
// Input/Output function for STL format.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v3.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "bitpit_common.hpp"
# include "STL.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //
using namespace std;

namespace bitpit{

/*!
    \class STLObj
    \brief Interface to STL I/O function

    This class has been designed to allow an easy interface between end-user
    and STL I/O functions.
*/

// Class STLObj methods ==================================================== //

// Constructors ------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Default constructor for class STLObj.

    Initialize an empty interface to stl file.
*/
STLObj::STLObj(
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
// SET DEFAULT VALUES                                                         //
// ========================================================================== //

// General info
stl_name = "";
stl_type = false;

// Error flags
err = 0;

// stl content
data.n_solids = -1;

return; };

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class STLObj.
    Initialize an interface to stl file with name specified in filename.

    \param[in] filename stl file name
    \param[in] filetype boolean flag for ascii (false) or binary (true) stl file
*/
STLObj::STLObj(
    string                              filename,
    bool                                filetype
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET DEFAULT VALUES                                                         //
// ========================================================================== //

// General info
stl_name = utils::string::trim(filename);
stl_type = filetype;

// Error flags
err = 0;

// stl content
data.n_solids = -1;

return; };

// -------------------------------------------------------------------------- //
/*!
    Constructor #1 for class STLObj.
    Initialize an interface to stl file with name specified in filename.

    \param[in] filename stl file name
*/
STLObj::STLObj(
    string                              filename
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// SET DEFAULT VALUES                                                         //
// ========================================================================== //

// General info
stl_name = utils::string::trim(filename);

// Error flags
err = 0;

// stl content
data.n_solids = -1;

// ========================================================================== //
// DETECT FILE TYPE (ASCII OR BINARY)                                         //
// ========================================================================== //

// Open file
ifile_handle.open(stl_name, ifstream::in);

// Check if the file start with the work "SOLID"
char keyword[] = "solid";
int keywordSize = sizeof(keyword) / sizeof(keyword[0]) - 1;

char c;
int nMatches = 0;
while (ifile_handle.get(c)) {
    if (isblank(c)) {
        continue;
    }

    c = tolower(c);
    if (c == keyword[nMatches]) {
        ++nMatches;
        if (nMatches == keywordSize) {
            break;
        }
    } else {
        break;
    }
}

stl_type = (nMatches != keywordSize);

// Close file
ifile_handle.close();

return; };

// Public methods ----------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Open the file associated to the interface.

    \param[in] mode opening mode ("in": input, "out": output, "app": append mode)
*/
void STLObj::open(
    string                              mode
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

// Open input stream -------------------------------------------------------- //
if (mode.compare("in") == 0) {
    if (!ifile_handle.is_open()) {

        // Close output stream
        close("out");

        // Open input stream
        ifile_handle.open(stl_name, ifstream::in | ifstream::binary);

        // Check stream status
        if (!ifile_handle.good()) { err = 1; return; }
        else                      { err = 0; return; }

    }
}

// Open output stream in "write" mode --------------------------------------- //
else if (mode.compare("out") == 0) {
    if (!ofile_handle.is_open()) {

        // Close input stream
        close("in");

        // Clear infos & flags
        clear();

        // Open output stream
        if (stl_type) { ofile_handle.open(stl_name, ifstream::out | ifstream::binary); }
        else          { ofile_handle.open(stl_name, ifstream::out); }
    }

    // Check stream status
    if (!ofile_handle.good()) { err = 1; return; }
    else                      { err = 0; return; }
}

// Open output stream in "append" mode -------------------------------------- //
else if (mode.compare("app") == 0) {
    if (!ofile_handle.is_open()) {

        // Close input stream
        close("in");

        // Open output stream
        if (stl_type) { ofile_handle.open(stl_name, ifstream::app | ifstream::binary); }
        else          { ofile_handle.open(stl_name, ifstream::app); }
    }

    // Check stream status
    if (!ofile_handle.good()) { err = 1; return; }
    else                      { err = 0; return; }
}

return; };

// -------------------------------------------------------------------------- //
/*!
    Close the current stream to stl file.

    \param[in] mode opening mode used to open stream ("in": input, "out": output,
    "app": append)
*/
void STLObj::close(
    string                              mode
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

// Close input stream ------------------------------------------------------- //
if (mode.compare("in") == 0) {
    ifile_handle.close();
}
else if (mode.compare("out") == 0) {
    ofile_handle.close();
}
else if (mode.compare("app") == 0) {
    ofile_handle.close();
}
else if (mode.compare("") == 0) {
    close("in");
    close("out");
}

return; };

// -------------------------------------------------------------------------- //
/*!
    Clear info and error flags gathered on the stl file associated
    to the interface.
*/
void STLObj::clear(
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
// RESET VALUES TO DEFAULT                                                    //
// ========================================================================== //

// Close input/output stream ------------------------------------------------ //
close();

// Set members values to default -------------------------------------------- //

// Error flags
err = 0;
stl_errors.resize(0);

// stl content
data.n_solids = -1;
data.solid_names.resize(0);
data.solid_facets.resize(0);

return; }

// -------------------------------------------------------------------------- //
/*!
    Display info gathered on the stl file associated to the interface.

    \param[in,out] out output stream
*/
void STLObj::display(
    ostream                             &out
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int                 i;

// ========================================================================== //
// DISPLAY INFO                                                               //
// ========================================================================== //

// General info ------------------------------------------------------------- //
out << "STL object:" << endl;
out << "stl name          : " << stl_name << endl;
if (stl_type) { out << "file type         : binary" << endl; }
else          { out << "file type         : ASCII"  << endl; }
if (ifile_handle.is_open()) { out << "input stream      : open" << endl; }
else                        { out << "input stream      : closed" << endl; }
if (ofile_handle.is_open()) { out << "output stream     : open" << endl; }
else                        { out << "output stream     : closed" << endl; }

// Error flags -------------------------------------------------------------- //
if (err == 1) { out << "**ERROR** stl file is missing!" << endl; return; }

// Stl content -------------------------------------------------------------- //
if (data.n_solids >= 0) {
    out << "# of stl solids   : " << data.n_solids << endl;
    if (data.n_solids > 0) {
        out << "solid's name      : " << data.solid_names << endl;
        out << "# of facet        : " << data.solid_facets    << endl;
    }
}

// Errors ------------------------------------------------------------------- //

if ((data.n_solids >= 0) && (stl_errors.size() > 0)) {
    out << "ERROR report:" << endl;
    for (i = 0; i < data.n_solids; i++) {
        out << "  stl solid: '" << data.solid_names[i] << "'" << endl;
        if (stl_errors[i][0]) {
            out << "    **WARNING** unterminated solid block!!" << endl;
        }
        if (stl_errors[i][1]) {
            out << "    **WARNING** unterminated facet block!!" << endl;
        }
        if (stl_errors[i][2]) {
            out << "    **WARNING** normal data are missing!!" << endl;
        }
        if (stl_errors[i][3]) {
            out << "    **WARNING** wrong number of components for normal data!!" << endl;
        }
        if (stl_errors[i][4]) {
            out << "    **WARNING** wrong number of vertices in facet block!!" << endl;
        }
        if (stl_errors[i][5]) {
            out << "    **WARNING** wrong number of coordinates for vertices!!" << endl;
        }
    } //next i
}

return; };

// -------------------------------------------------------------------------- //
/*!
    Scan and gather info from the stl file associated to the interface
*/
void STLObj::scan(
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
// CLEAR INFOS                                                                //
// ========================================================================== //
clear();

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// SCAN STL FILE                                                              //
// ========================================================================== //
if (stl_type) {
    err = stl::scanBINARY(ifile_handle, data.solid_names, data.solid_facets);
    data.n_solids = data.solid_names.size();
}
else {
    err = stl::scanASCII(ifile_handle, data.solid_names, data.solid_facets);
    data.n_solids = data.solid_names.size();
}

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; }

// -------------------------------------------------------------------------- //
/*!
    Scan and check for format errors in the stl file associated to the interface
*/
void STLObj::check(
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
// SCAN STL FILE                                                              //
// ========================================================================== //
scan();

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");

// ========================================================================== //
// CHECK STL FILE                                                             //
// ========================================================================== //
if (stl_type) { stl::checkBINARY(ifile_handle, stl_errors); }
else          { stl::checkASCII(ifile_handle, stl_errors); }

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Load solid data from the stl file associated to the interface.

    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
*/
void STLObj::load(
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T
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
if (err == 1) { return; }

// ========================================================================== //
// READ STL DATA                                                              //
// ========================================================================== //
if (stl_type)   { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else            { stl::readASCII(ifile_handle, nV, nT, V, N, T); }

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Load solid data from the stl file associated to the interface. Overloading of
    member function STLObj::load() for container vector<array<double, 3>>

    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
*/
void STLObj::load(
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int, 3> >              &T
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
if (err == 1) { return; }

// ========================================================================== //
// READ STL DATA                                                              //
// ========================================================================== //
if (stl_type)   { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else            { stl::readASCII(ifile_handle, nV, nT, V, N, T); }

// ========================================================================== //
// CLOSE INPUT STREAM                                                         //
// ========================================================================== //
close("in");

return; };

// -------------------------------------------------------------------------- //
/*!
    Load single solid data from the stl file associated to the interface.
    This routine assumes that the input stream is already open.

    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
*/
void STLObj::loadSolid(
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// READ STL DATA                                                              //
// ========================================================================== //
if (stl_type)   { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else            { stl::readSolidASCII(ifile_handle, false, nV, nT, V, N, T); }

return; };

// -------------------------------------------------------------------------- //
/*!
    Load single solid solid data from the stl file associated to the interface.
    This routine assumes that the input stream is already open.
    Overloading of member function STLObj::load() for container
    vector<array<double, 3>>

    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
*/
void STLObj::loadSolid(
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int,3> >               &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// READ STL DATA                                                              //
// ========================================================================== //
if (stl_type)   { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else            { stl::readSolidASCII(ifile_handle, false, nV, nT, V, N, T); }

return; };


// Private methods ---------------------------------------------------------- //

// -------------------------------------------------------------------------- //
/*!
    Dummy function for recursive variadic template STLObj::save().
*/
void STLObj::save(
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
close("out");

return; };

// -------------------------------------------------------------------------- //
/*!
    Dummy function for recursive variadic template STLObj::load().
*/
void STLObj::load(
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

return; };

// Scanning routines ======================================================== //

// -------------------------------------------------------------------------- //
/*!
    Scan ascii stl file and retrieve infos about content.

    \param[in,out] file_handle input stream from stl file
    \param[in,out] solid_names list of solid names stored in the stl file
    \param[in] solid_facets number of facets for each stl solid found in the
    stl file (same ordering of solid_names).

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to scan stl file.
*/
unsigned int stl::scanASCII(
    ifstream                            &file_handle,
    vector<string>                      &solid_names,
    vector<int>                         &solid_facets
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int             start_pos;
string               line, word;
stringstream         sline, sword;

// Counters
int                  nF = 0;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
solid_names.resize(0);
solid_facets.resize(0);

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) {
    return(1);
}

// ========================================================================== //
// RESET CURSOR AT FILE BEGIN                                                 //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// SCAN STL FILE.                                                             //
// ========================================================================== //

// Set cursor position at file begin ---------------------------------------- //
file_handle.clear();
file_handle.seekg(0);

// Scan until eof is reached
while (!file_handle.eof()) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Get keyword
    if (sline >> word) {
        if (word.compare("solid") == 0) {

            // Reset counters
            nF = 0;

            // Get solid name
            sword.str("");
            while (sline >> word) {
                sword << word << " ";
            }
            word = sword.str();
            solid_names.push_back(utils::string::trim(word));

            // Get solid info
            stl::scanSolidASCII(file_handle, nF);
            solid_facets.push_back(nF);
            
        }
    }

} //next line

// Restore cursor position -------------------------------------------------- //
file_handle.clear();
file_handle.seekg(start_pos);

return(0); };

// -------------------------------------------------------------------------- //
/*!
    Scan binary stl file and retrieve infos about content.

    \param[in,out] file_handle input stream from stl file
    \param[in,out] solid_names list of solid names stored in the stl file
    \param[in] solid_facets number of facets for each stl solid found in the
    stl file (same ordering of solid_names).

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to scan stl file.
*/
unsigned int stl::scanBINARY(
    ifstream                            &file_handle,
    vector<string>                      &solid_names,
    vector<int>                         &solid_facets
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int           start_pos = file_handle.tellg();
unsigned long int  longint4byte;
unsigned long int *longint4byte_ = &longint4byte;
float              float4byte;
float             *float4byte_   = &float4byte;

// Counters
int                i;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
solid_names.resize(0);
solid_facets.resize(0);

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); };

// ========================================================================== //
// SCAN STL FILE                                                              //
// ========================================================================== //

// Set cursor position at file begin
file_handle.seekg(0);

// Get solid name
for (i = 0; i < 20; i++) {
    file_handle.read(reinterpret_cast<char*>(float4byte_), 4);
} //next i
solid_names.push_back("");

// Read number of facets
file_handle.read(reinterpret_cast<char*>(longint4byte_), 4);
solid_facets.push_back((int)*longint4byte_);

// Reset cursor position
file_handle.clear();
file_handle.seekg(start_pos);

return(0); };

// -------------------------------------------------------------------------- //
/*!
    Scan solid data from ascii stl file and retrieve infos.

    \param[in,out] file_handle input stream from stl file
    \param[in,out] nT number of facets found in the solid data section

    \result returns an error flag for I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to scan stl file.
*/
unsigned int stl::scanSolidASCII(
    ifstream                            &file_handle,
    int                                 &nT
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int        backup_pos;
string          line, word;
stringstream    sline;

// Counters
// none

// ========================================================================== //
// RESET COUNTERS                                                             //
// ========================================================================== //
nT = 0;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); };

// ========================================================================== //
// SCAN STL SOLID                                                             //
// ========================================================================== //

// Get current line
backup_pos = file_handle.tellg();
getline(file_handle, line);
line = utils::string::trim(line);
sline.clear();
sline.str(line);
if (!(sline >> word)) { word = ""; }

// Scan lines until eof or keyword "endsolid" is found
while ((!file_handle.eof())
    && ((word.compare("endsolid") != 0)
    &&  (word.compare("solid") != 0))) {

    // Look for keyword "facet"
    if (word.compare("facet") == 0) {
        nT++;
    }

    // Get next line
    backup_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if (!(sline >> word))  { word = ""; }

} //next line
if (word.compare("endsolid") != 0) {
    file_handle.clear();
    file_handle.seekg(backup_pos);
}


return(0); };

// Check routines =========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Scan ascii stl file and perform check on format error.

    \param[in] file_handle input stream from stl file.
    \param[in,out] err_map for each solid found in the stl file, store a bitmask
    for format errors encountered in the solid data section.
    error code:
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check stl file
*/
unsigned int stl::checkASCII(
    ifstream                            &file_handle,
    vector<vector<bool> >               &err_map
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int            start_pos = file_handle.tellg();
string              line, word;
stringstream        sline;

// Counters
int                 n_solid = 0;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
err_map.resize(0);

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); };

// ========================================================================== //
// CHECK STL FILE                                                             //
// ========================================================================== //

// Set cursor position at file begin
file_handle.clear();
file_handle.seekg(0);

// Scan stl solid
while (getline(file_handle, line)) {

    // Get current line
    line = utils::string::trim(line);
    sline.clear(),
    sline.str(line);

    if ((sline >> word) && (word.compare("solid") == 0)) {
        vector<bool>    _map(6, false);
        n_solid++;
        stl::checkSolidASCII(file_handle, _map);
        err_map.push_back(_map);
    }

} //next line

// Reset cursor position
file_handle.clear();
file_handle.seekg(start_pos);

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Scan solid data section in ascii stl file and perform check on format error.

    \param[in] file_handle input stream from stl file.
    \param[in,out] err_map bitmask for format errors encountered in the solid data section.
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check stl file
*/
unsigned int stl::checkSolidASCII(
    ifstream                            &file_handle,
    vector<bool>                        &err_map
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int            backup_pos;
string              line, word;
stringstream        sline;

// Counters
// none

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
err_map.resize(6, false);

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN SOLID DATA                                                            //
// ========================================================================== //

// Get current line
backup_pos = file_handle.tellg();
getline(file_handle, line);
line = utils::string::trim(line);
sline.clear();
sline.str(line);
if (!(sline >> word)) { word = ""; }

// Scan stl solid data until eof or keyword "endsolid" is found
while ((!file_handle.eof())
    && ((word.compare("endsolid") != 0) 
    &&  (word.compare("solid")    != 0))) {

        // Look for keyword "facet"
        if (word.compare("facet") == 0) {
            file_handle.seekg(backup_pos);
            stl::checkFacetASCII(file_handle, err_map);
        }

        // Get next line
        backup_pos = file_handle.tellg();
        getline(file_handle, line);
        line = utils::string::trim(line);
        sline.clear();
        sline.str(line);
        if (!(sline >> word)) { word = ""; }

} //next line

// Check block temination
if (word.compare("endsolid") != 0) {
    err_map[0] = true;
    file_handle.clear();
    file_handle.seekg(backup_pos);
};

return(0); };

// -------------------------------------------------------------------------- //
/*!
    Scan facet data section in ascii stl file and perform check on format error.

    \param[in] file_handle input stream from stl file.
    \param[in,out] err_map bitmask for format errors encountered in the facet data section.
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check stl file
*/
unsigned int stl::checkFacetASCII(
    ifstream                            &file_handle,
    vector<bool>                        &err_map
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool            check_normal = false;
long int        backup_pos;
string          line, word;
stringstream    sline;

// Counters
int             nxyz = 0, nV = 0;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
err_map.resize(6, false);

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN FACET DATA                                                            //
// ========================================================================== //

// Get current line --------------------------------------------------------- //
backup_pos = file_handle.tellg();
getline(file_handle, line);
line = utils::string::trim(line);
sline.clear(),
sline.str(line);
if ((!(sline >> word)) || (word.compare("facet") == 0)) { word = "begin"; }

// Scan facet data until eof or keyword "endfacet" is found ----------------- //
while ((!file_handle.eof()) &&
    ((word.compare("endfacet")  != 0)
  && (word.compare("facet")     != 0)
  && (word.compare("endsolid")  != 0)
  && (word.compare("solid")     != 0))) {

    // Check facet normal
    if (word.compare("begin") == 0) {
        if ((sline >> word) && (word.compare("normal") == 0)) {
            check_normal = true;
            nxyz = 0;
            while (sline >> word) {
                nxyz++;
            } //next field
            if (nxyz != 3) { err_map[3] = true; }
        }
    }
    else if (word.compare("vertex") == 0) {
        nV++;
        nxyz = 0;
        while (sline >> word) {
            nxyz++;
        } //next field
        if (nxyz != 3) { err_map[5] = true; }
    }

    // Get next line
    backup_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear(),
    sline.str(line);
    if (!(sline >> word)) { word = ""; }

} //next line

// Check error status

// Block termination
if (word.compare("endfacet") != 0) {
    err_map[1] = true;
    file_handle.clear(),
    file_handle.seekg(backup_pos);
}

// Normal data
if (!check_normal) { err_map[2] = true; }

// Number of vertices
if (nV != 3) { err_map[4] = true; }


return(0); }

// -------------------------------------------------------------------------- //
/*!
    Scan binary stl file and perform check on format error.

    \param[in] file_handle input stream from stl file.
    \param[in,out] err_map for each solid found in the stl file, store a bitmask
    for format errors encountered in the solid data section.
    error code:
        0 -> unterminated solid block
        1 -> undeterminated facet block
        2 -> normal data are missing
        3 -> wrong number of components for facet normal
        4 -> wrong number of vertices in facet block
        5 -> wrong number of coordinates for vertex

    \result returns and error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to scan/check stl file
*/
unsigned int stl::checkBINARY(
    ifstream                            &file_handle,
    vector<vector<bool> >               &err_map
) {

// ========================================================================== //
// unsigned int checkBINARY(                                                //
//     ifstream                  &file_handle,                                //
//     bvector2D                 &err_map)                                    //
//                                                                            //
// Check data consistency in binary stl files.                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle : ifstream, input stream to stl file                         //
// - err_map     : bvector2D, list of errors encountered for each solid block //
//                 while checking stl file.                                   //
//                 err_map[i][err_code] for the i-th stl solid is true        //
//                 if an error corresponding to error code "err_code"         //
//                 is detected.                                               //
//                 err_code = 0   --> unterminated solid block                //
//                 err_code = 1   --> unterminated facet block                //
//                 err_code = 2   --> normal data are missing                 //
//                 err_code = 3   --> wrong number of components for normal   //
//                 err_code = 4   --> wrong number of vertices in facet block //
//                 err_code = 5   --> wrong number of coordinates for vertex  //
//                                    data.                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err         : unsigned int, error flag:                                  //
//                 err = 0     --> no errors encountered                      //
//                 err = 1     --> file is missing or is corrupted            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int       int2byte;
unsigned int      *int2byte_ = &int2byte;
unsigned long int  longint4byte;
unsigned long int *longint4byte_ = &longint4byte;
float              float4byte;
float             *float4byte_   = &float4byte;

// Counters
int                i, j, k, nT;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //
err_map.resize(1, vector<bool>(6, false));

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// CHECK STL FILE                                                             //
// ========================================================================== //

// Check file header -------------------------------------------------------- //
for (i = 0; i < 20; i++) {
    file_handle.read(reinterpret_cast<char*>(float4byte_), 4);
} //next i
if (file_handle.eof()) {
    err_map[0][0] = true;
    return(0);
}

// Check number of facets --------------------------------------------------- //
file_handle.read(reinterpret_cast<char*>(longint4byte_), 4);
nT = (int)*longint4byte_;
if (file_handle.eof()) {
    err_map[0][0] = true;
    return(0);
}

// Check facet data --------------------------------------------------------- //
i = 0;
while ((!file_handle.eof()) && (i < nT)) {

    // Read normal
    for (j = 0; j < 3; j++) {
        file_handle.read(reinterpret_cast<char*>(float4byte_), 4);
    } //next j

    // Read vertex coordinates
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            file_handle.read(reinterpret_cast<char*>(float4byte_), 4);
        } //next k
    } //next j

    // Facet closing header
    file_handle.read(reinterpret_cast<char*>(int2byte_), 2);

    // Update facet counter
    i++;

} //next i

// Check error status ------------------------------------------------------- //
if (i < nT) { err_map[0][1] = true; };


return(0); }

// Input routines =========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Read STL facet data from ascii stl file.

    \param[in,out] file_handle input stream from stl file
    \param[in,out] nV on input stores the number of vertices previously acquired from
    the stl file, on output stores the input value incremented by the number
    of vertices read from the facet section.
    \param[in,out] nT on input stores the number of facets previously aquired from
    the stl file, on output stores the input value incremented by 1.
    \param[in,out] V vertex coordinate list. Vertices read from facet section
    are stored in the location V[nV], V[nV+1], etc. nV is the number of vertices
    previously acquired from the stl file and passed as input to this function.
    \param[in,out] N normal list. Normal read from facet section is stored at
    position N[nT], where nT is the number of facet previously acquired from the
    stl file and passed as input to this function.
    \param[in,out] T vertex->facet connectivity. A new connectivity entry for the facet
    begin read is created at position T[nT], where nT is the number of facet previously
    acquired from the stl file and passed as input to this function.

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: one facet has more than 3 vertices
*/  
unsigned int stl::readFacetASCII(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int            cursor_pos;
string              line, word;
stringstream        sline;

// Counters
int                 nv = 0;

// ========================================================================== //
// READ FACET DATA                                                            //
// ========================================================================== //

// Read facet data ---------------------------------------------------------- //

// Get current line
cursor_pos = file_handle.tellg();
getline(file_handle, line);
line = utils::string::trim(line);
sline.clear();
sline.str(line);
if ((!(sline >> word)) || (word.compare("facet") == 0)) { word = "begin"; }

// Loop until end of facet is found
while ((!file_handle.eof())
    && ((word.compare("endfacet")   != 0)
    &&  (word.compare("facet")      != 0)
    &&  (word.compare("solid")      != 0)
    &&  (word.compare("endsolid")   != 0))) {

    // Look for keywords
    if (word.compare("begin") == 0) {
        if ((sline >> word) && (word.compare("normal") == 0)) {
            sline >> N[nT];
        }
    }
    else if (word.compare("vertex") == 0) {
        if (nv >= 3) {
            return(2);
        }

        sline >> V[nV+nv];
        nv++;
    }

    // Get next line
    cursor_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if (!(sline >> word)) { word = ""; }

} //next line

// Restor cursor position --------------------------------------------------- //
if (word.compare("endfacet") != 0) {
    file_handle.clear();
    file_handle.seekg(cursor_pos);
}

// Update triangle-vertex connectivity -------------------------------------- //
for (int i = 0; i < nv; i++) {
    T[nT][i] = nV + i;
} //next i

// Update facet/vertex counters --------------------------------------------- //
nV += 3;
nT++;

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read STL facet data from ascii stl file.
    Overloading of readFacetASCII() for vector<array<double, 3> > container.

    \param[in,out] file_handle input stream from stl file
    \param[in,out] nV on input stores the number of vertices previously acquired from
    the stl file, on output stores the input value incremented by the number
    of vertices read from the facet section.
    \param[in,out] nT on input stores the number of facets previously aquired from
    the stl file, on output stores the input value incremented by 1.
    \param[in,out] V vertex coordinate list. Vertices read from facet section
    are stored in the location V[nV], V[nV+1], etc. nV is the number of vertices
    previously acquired from the stl file and passed as input to this function.
    \param[in,out] N normal list. Normal read from facet section is stored at
    position N[nT], where nT is the number of facet previously acquired from the
    stl file and passed as input to this function.
    \param[in,out] T vertex->facet connectivity. A new connectivity entry for the facet
    begin read is created at position T[nT], where nT is the number of facet previously
    acquired from the stl file and passed as input to this function.

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: one facet has more than 3 vertices
*/  
unsigned int stl::readFacetASCII(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int, 3> >              &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int            cursor_pos;
string              line, word;
stringstream        sline;

// Counters
int                 nv = 0;

// ========================================================================== //
// READ FACET DATA                                                            //
// ========================================================================== //

// Read facet data ---------------------------------------------------------- //

// Get current line
cursor_pos = file_handle.tellg();
getline(file_handle, line);
line = utils::string::trim(line);
sline.clear();
sline.str(line);
if ((!(sline >> word)) || (word.compare("facet") == 0)) { word = "begin"; }

// Loop until end of facet is found
while ((!file_handle.eof())
    && ((word.compare("endfacet")   != 0)
    &&  (word.compare("facet")      != 0)
    &&  (word.compare("solid")      != 0)
    &&  (word.compare("endsolid")   != 0))) {

    // Look for keywords
    if (word.compare("begin") == 0) {
        if ((sline >> word) && (word.compare("normal") == 0)) {
            sline >> N[nT];
        }
    }
    else if (word.compare("vertex") == 0) {
        if (nv >= 3) {
            return(2);
        }

        sline >> V[nV+nv];
        nv++;
    }

    // Get next line
    cursor_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);
    if (!(sline >> word)) { word = ""; }

} //next line

// Restor cursor position --------------------------------------------------- //
if (word.compare("endfacet") != 0) {
    file_handle.clear();
    file_handle.seekg(cursor_pos);
}

// Update triangle-vertex connectivity -------------------------------------- //
for (int i = 0; i < nv; i++) {
    T[nT][i] = nV + i;
} //next i

// Update facet/vertex counters --------------------------------------------- //
nV += 3;
nT++;

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read solid data from ascii stl file.

    \param[in,out] file_handle stream from stl file
    \param[in] wrapAround controls if, when the end of file is reached, the
    search for the solid to read will begin again from the beginning of the
    file.
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the solid.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the solid.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the solid. Vertices acquired from solid are appended at the end
    of V.
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the solid. Normals acquired from solid are appended at the end
    of N.
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T.
    \param[in,out] solid_name name of the stl solid (if available) read from the stl file

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: failed to read a facet
*/
unsigned int stl::readSolidASCII(
    ifstream                            &file_handle,
    bool                                wrapAround,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T,
    string                              solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = false;
long int            start_pos = file_handle.tellg(), current_pos;
string              line, word;
stringstream        sline;

// Counters
int                 nt;

// ========================================================================== //
// CHECK INPUT STREAM                                                         //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN STL FILE UNTIL A SOLID IS FOUND                                       //
// ========================================================================== //

// Parameters --------------------------------------------------------------- //
sline << "solid " << solid_name;
line = sline.str();
solid_name = utils::string::trim(line);



// Scan file until stl solid is found --------------------------------------- //
current_pos = start_pos+1;
while (!check && (start_pos != current_pos)) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Check end of file
    if (file_handle.eof()) {
        if (wrapAround) {
            file_handle.clear();
            file_handle.seekg(0);
        } else {
            check = false;
            break;
        }
    }
    current_pos = file_handle.tellg();

    // Look for keyword "solid"
    if ((sline >> word) && (word.compare("solid") == 0)) {
        if (solid_name.compare("solid") == 0) {
            start_pos = current_pos;
            check = true;
        }
        else {
            if (line.compare(solid_name) == 0) {
                start_pos = current_pos;
                check = true;
            }
        }
    }

} //next line
if (!check) { return(0); }

// Scan stl solid ----------------------------------------------------------- //
file_handle.clear();
file_handle.seekg(start_pos);
stl::scanSolidASCII(file_handle, nt);

// ========================================================================== //
// READ SOLID DATA                                                            //
// ========================================================================== //

// Resize input variables --------------------------------------------------- //
V.resize(nV+3*nt, vector<double>(3, 0.0));
N.resize(nT+nt, vector<double>(3, 0.0));
T.resize(nT+nt, vector<int>(3, -1));

// Read solid data ---------------------------------------------------------- //
file_handle.clear();
file_handle.seekg(start_pos);
word = "begin";
while ((!file_handle.eof())
    && (word.compare("endsolid") != 0)
    && (word.compare("solid") != 0)) {

    // Get current line
    current_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for keyword "facet"
    if ((sline >> word) && (word.compare("facet") == 0)) {
        file_handle.seekg(current_pos);
        int readFacetError = stl::readFacetASCII(file_handle, nV, nT, V, N, T);
        if (readFacetError != 0) {
            return(2);
        }
    }
} //next line
if (word.compare("endsolid") != 0) {
    file_handle.clear();
    file_handle.seekg(current_pos);
}

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read solid data from ascii stl file. Overloading of readSolidASCII()
    for vector<array<double, 3> > container.

    \param[in,out] file_handle stream from stl file
    \param[in] wrapAround controls if, when the end of file is reached, the
    search for the solid to read will begin again from the beginning of the
    file.
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the solid.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the solid.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the solid. Vertices acquired from solid are appended at the end
    of V.
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the solid. Normals acquired from solid are appended at the end
    of N.
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T.
    \param[in,out] solid_name name of the stl solid (if available) read from the stl file

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
        err = 2: failed to read a facet
*/
unsigned int stl::readSolidASCII(
    ifstream                            &file_handle,
    bool                                wrapAround,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int,3> >               &T,
    string                               solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = false;
long int            start_pos = file_handle.tellg(), current_pos;
string              line, word;
stringstream        sline;

// Counters
int                 nt;

// ========================================================================== //
// CHECK INPUT STREAM                                                         //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN STL FILE UNTIL A SOLID IS FOUND                                       //
// ========================================================================== //

// Parameters --------------------------------------------------------------- //
sline << "solid " << solid_name;
line = sline.str();
solid_name = utils::string::trim(line);

// Scan file until stl solid is found --------------------------------------- //
current_pos = start_pos+1;
while (!check && (start_pos != current_pos)) {

    // Get current line
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Check end of file
    if (file_handle.eof()) {
        if (wrapAround) {
            file_handle.clear();
            file_handle.seekg(0);
        } else {
            check = false;
            break;
        }
    }
    current_pos = file_handle.tellg();

    // Look for keyword "solid"
    if ((sline >> word) && (word.compare("solid") == 0)) {
        if (solid_name.compare("solid") == 0) {
            start_pos = current_pos;
            check = true;
        }
        else {
            if (line.compare(solid_name) == 0) {
                start_pos = current_pos;
                check = true;
            }
        }
    }

} //next line
if (!check) { return(0); }

// Scan stl solid ----------------------------------------------------------- //
file_handle.clear();
file_handle.seekg(start_pos);
stl::scanSolidASCII(file_handle, nt);

// ========================================================================== //
// READ SOLID DATA                                                            //
// ========================================================================== //

// Resize input variables --------------------------------------------------- //
array<double,3>         dummyDoubleArray;
dummyDoubleArray.fill(0.) ;

array<int,3>            dummyIntArray;
dummyIntArray.fill(0) ;

V.resize(nV+3*nt, dummyDoubleArray);
N.resize(nT+nt, dummyDoubleArray);
T.resize(nT+nt, dummyIntArray);

// Read solid data ---------------------------------------------------------- //
file_handle.clear();
file_handle.seekg(start_pos);
word = "begin";
while ((!file_handle.eof())
    && (word.compare("endsolid") != 0)
    && (word.compare("solid") != 0)) {

    // Get current line
    current_pos = file_handle.tellg();
    getline(file_handle, line);
    line = utils::string::trim(line);
    sline.clear();
    sline.str(line);

    // Look for keyword "facet"
    if ((sline >> word) && (word.compare("facet") == 0)) {
        file_handle.seekg(current_pos);
        int readFacetError = readFacetASCII(file_handle, nV, nT, V, N, T);
        if (readFacetError != 0) {
            return(2);
        }
    }
} //next line
if (word.compare("endsolid") != 0) {
    file_handle.clear();
    file_handle.seekg(current_pos);
}

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read data from ascii stl file.

    \param[in,out] file_handle stream from stl file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the stl file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the stl file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the stl file. Vertices acquired from the stl file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the stl file. Normals acquired from the stl file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same stl file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int stl::readASCII(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int        start_pos, current_pos;
stringstream    sline;
string          line, word;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(0); }

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// LOAD STL SOLIDS                                                            //
// ========================================================================== //
//while (!file_handle.eof()) {
current_pos = 0;
while (getline(file_handle, line)) {

    // Get current line
    sline.clear();
    sline.str(line);
    if ((sline >> word) && (word.compare("solid") == 0)) {
        file_handle.clear();
        file_handle.seekg(current_pos);
        readSolidASCII(file_handle, true, nV, nT, V, N, T, "");
    }
    current_pos = file_handle.tellg();
}
//} //next solid

// ========================================================================== //
// RESET CURSOR POSITION AT FILE BEGIN                                        //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(0); }


// -------------------------------------------------------------------------- //
/*!
    Read data from ascii stl file. Overloading of readASCII()
    for container vector<array<double,3> >

    \param[in,out] file_handle stream from stl file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the stl file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the stl file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the stl file. Vertices acquired from the stl file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the stl file. Normals acquired from the stl file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same stl file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int stl::readASCII(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int,3> >               &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int        start_pos, current_pos;
stringstream    sline;
string          line, word;

// Counters
// none

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(0); }

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// LOAD STL SOLIDS                                                            //
// ========================================================================== //
//while (!file_handle.eof()) {
current_pos = 0;
while (getline(file_handle, line)) {

    // Get current line
    sline.clear();
    sline.str(line);
    if ((sline >> word) && (word.compare("solid") == 0)) {
        file_handle.clear();
        file_handle.seekg(current_pos);
        stl::readSolidASCII(file_handle, true, nV, nT, V, N, T, "");
    }
    current_pos = file_handle.tellg();
}
//} //next solid

// ========================================================================== //
// RESET CURSOR POSITION AT FILE BEGIN                                        //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read data from binary  stl file.

    \param[in,out] file_handle stream from stl file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the stl file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the stl file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the stl file. Vertices acquired from the stl file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the stl file. Normals acquired from the stl file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same stl file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int stl::readBINARY(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int           start_pos;
unsigned long int  longint4byte;
unsigned int       int2byte;
float              float4byte;
unsigned int      *int2byte_pointer     = &int2byte;
unsigned long int *longint4byte_pointer = &longint4byte;
float             *float4byte_pointer   = &float4byte;

// Counters
int                i, j, k;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) {
    file_handle.close();
    return(1);
}

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// SCAN _BINARY STL                                                            //
// ========================================================================== //

// Title
for (i = 0; i < 20; i++) {
    file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
} //next i

// Read number of elements
file_handle.read(reinterpret_cast<char*>(longint4byte_pointer), 4);
nT = (int)*longint4byte_pointer;
nV = 3*nT;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

// Vertex list
V.resize(nV, vector<double>(3, 0.0));

// Normals
N.resize(nT, vector<double>(3, 0.0));

// Triangle-vector connectivity
T.resize(nT, vector<int>(3, -1));

// ========================================================================== //
// READ DATA                                                                  //
// ========================================================================== //
nV = 0;
for (i = 0; i < nT; i++) {

    // Read normal
    for (j = 0; j < 3; j++) {
        file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
        N[i][j] = (double)*float4byte_pointer;
    } //next j

    // Read vertex coordinates
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
            V[nV][k] = (double)*float4byte_pointer;
        } //next k
        nV++;
    } //next j

    // Triangle-vertex connectivity
    T[i][0] = nV - 3;
    T[i][1] = nV - 2;
    T[i][2] = nV - 1;

    // Facet closing header
    file_handle.read(reinterpret_cast<char*>(int2byte_pointer), 2);

} //next i

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read data from binary stl file. Overloading of readBINARY() for container
    vector<array<double,3> >

    \param[in,out] file_handle stream from stl file
    \param[in,out] nV on input stores the number of vertices currently stored in V.
    On output stores the input value increased by the number of vertices acquired
    from the stl file.
    \param[in,out] nT on input stores the number of facets->vertex connectivity
    entries currently stored in T. On output stores the input value increased
    by the number of facet acquired from the stl file.
    \param[in,out] V vertex coordinate list. On output stores the coordinates of vertices
    acquired from the stl file. Vertices acquired from the stl file are appended at the end
    of V (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] N facet's normals. On output stores the normal unit vector
    to each facet acquired from the stl file. Normals acquired from the stl file are appended at the end
    of N (no distinction is made if multiple solids are stored in the same stl file).
    \param[in,out] T facet->vertex connectivity. On output stores facet->vertex connectivity
    for each facet acquired from the solid. New connectivity entries are appended at the end
    of T (no distinction is made if multiple solids are stored in the same stl file).

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to read from input stream
*/
unsigned int stl::readBINARY(
    ifstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int, 3> >              &T
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
long int           start_pos;
unsigned long int  longint4byte;
unsigned int       int2byte;
float              float4byte;
unsigned int      *int2byte_pointer     = &int2byte;
unsigned long int *longint4byte_pointer = &longint4byte;
float             *float4byte_pointer   = &float4byte;

// Counters
int                i, j, k;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) {
    file_handle.close();
    return(1);
}

// ========================================================================== //
// SET CURSOR POSITION AT FILE BEGIN                                          //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
file_handle.seekg(0);

// ========================================================================== //
// SCAN _BINARY STL                                                            //
// ========================================================================== //

// Title
for (i = 0; i < 20; i++) {
    file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
} //next i

// Read number of elements
file_handle.read(reinterpret_cast<char*>(longint4byte_pointer), 4);
nT = (int)*longint4byte_pointer;
nV = 3*nT;

// ========================================================================== //
// RESIZE INPUT VARIABLES                                                     //
// ========================================================================== //

array<double,3>     dummyDoubleArray;
dummyDoubleArray.fill(0.) ;

array<int,3>        dummyIntArray;
dummyIntArray.fill(0) ;

// Vertex list
V.resize(nV, dummyDoubleArray);

// Normals
N.resize(nT, dummyDoubleArray);

// Triangle-vector connectivity
T.resize(nT, dummyIntArray);

// ========================================================================== //
// READ DATA                                                                  //
// ========================================================================== //
nV = 0;
for (i = 0; i < nT; i++) {

    // Read normal
    for (j = 0; j < 3; j++) {
        file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
        N[i][j] = (double)*float4byte_pointer;
    } //next j

    // Read vertex coordinates
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            file_handle.read(reinterpret_cast<char*>(float4byte_pointer), 4);
            V[nV][k] = (double)*float4byte_pointer;
        } //next k
        nV++;
    } //next j

    // Triangle-vertex connectivity
    T[i][0] = nV - 3;
    T[i][1] = nV - 2;
    T[i][2] = nV - 1;

    // Facet closing header
    file_handle.read(reinterpret_cast<char*>(int2byte_pointer), 2);

} //next i

// ========================================================================== //
// RESET CURSOR POSITION                                                      //
// ========================================================================== //
file_handle.clear();
file_handle.seekg(start_pos);

return(0); }


// Output routines ========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Write solid data to ascii stl file.

    \param[in,out] file_handle stream to stl file
    \param[in,out] nV number of vertices to be written to stl file.
    \param[in,out] nT number of facet to be written to stl file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int stl::writeSolidASCII(
    ofstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T,
    string                              solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
stringstream     sheader;
string           header;
unsigned int     err = 0;

// Counters
int              n, m, i, j, k;

// =========================================================================== //
// CHECK INPUT COHERENCY                                                       //
// =========================================================================== //
if (V.size() < static_cast<unsigned int>(nV))  { return(2); };
if (T.size() < static_cast<unsigned int>(nT))  { return(2); };
if (N.size() < static_cast<unsigned int>(nT))  { return(2); };

// =========================================================================== //
// CHECK STREAM STATUS                                                         //
// =========================================================================== //
if (!file_handle.good()) { return(1); };

// ============================================================================ //
// WRITE STL SOLID                                                              //
// ============================================================================ //

// Solid header --------------------------------------------------------------- //
sheader << "solid " << solid_name;
header = sheader.str();
header = utils::string::trim(header);
file_handle << header << endl;
sheader.str("");

// Loop over triangulation facet ---------------------------------------------- //
for (i = 0; i < nT; i++) {

    // facet normal
    file_handle << "  facet";
    n = N[i].size();
    if (n > 0) {
        file_handle << " normal ";
        for (j = 0; j < n-1; j++) {
            file_handle << scientific << N[i][j] << " ";
        } //next j
        file_handle << scientific << N[i][n-1];;
    }
    file_handle << endl;

    // facet vertices
    file_handle << "    outer loop" << endl;
    n = T[i].size();
    for (j = 0; j < n; j++) {
        m = V[T[i][j]].size();
        file_handle << "      vertex ";
        if (m > 0) {
            for (k = 0; k < m-1; k++) {
                file_handle << scientific << V[T[i][j]][k] << " ";
            } //next k
            file_handle << scientific << V[T[i][j]][m-1];
        }
        file_handle << endl;
    } //next j
    file_handle << "    endloop"    << endl;

    // facet - closing header
    file_handle << "  endfacet" << endl;

} //next i

// Closing header ------------------------------------------------------------- //
sheader << "endsolid " << solid_name;
header = sheader.str();
header = utils::string::trim(header);
file_handle << header << endl;
sheader.str("");

return(err); };

// -------------------------------------------------------------------------- //
/*!
    Write solid data to ascii stl file. Overloading of writeSolidASCII() for
    container vector<array<double, 3> >.

    \param[in,out] file_handle stream to stl file
    \param[in,out] nV number of vertices to be written to stl file.
    \param[in,out] nT number of facet to be written to stl file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int stl::writeSolidASCII(
    ofstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int,3> >               &T,
    string                               solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
stringstream     sheader;
string           header;
unsigned int     err = 0;

// Counters
int              n, m, i, j, k;

// =========================================================================== //
// CHECK INPUT COHERENCY                                                       //
// =========================================================================== //
if (V.size() < static_cast<unsigned int>(nV))  { return(2); };
if (T.size() < static_cast<unsigned int>(nT))  { return(2); };
if (N.size() < static_cast<unsigned int>(nT))  { return(2); };

// =========================================================================== //
// CHECK STREAM STATUS                                                         //
// =========================================================================== //
if (!file_handle.good()) { return(1); };

// ============================================================================ //
// WRITE STL SOLID                                                              //
// ============================================================================ //

// Solid header --------------------------------------------------------------- //
sheader << "solid " << solid_name;
header = sheader.str();
header = utils::string::trim(header);
file_handle << header << endl;
sheader.str("");

// Loop over triangulation facet ---------------------------------------------- //
for (i = 0; i < nT; i++) {

    // facet normal
    file_handle << "  facet";
    n = N[i].size();
    if (n > 0) {
        file_handle << " normal ";
        for (j = 0; j < n-1; j++) {
            file_handle << scientific << N[i][j] << " ";
        } //next j
        file_handle << scientific << N[i][n-1];;
    }
    file_handle << endl;

    // facet vertices
    file_handle << "    outer loop" << endl;
    n = T[i].size();
    for (j = 0; j < n; j++) {
        m = V[T[i][j]].size();
        file_handle << "      vertex ";
        if (m > 0) {
            for (k = 0; k < m-1; k++) {
                file_handle << scientific << V[T[i][j]][k] << " ";
            } //next k
            file_handle << scientific << V[T[i][j]][m-1];
        }
        file_handle << endl;
    } //next j
    file_handle << "    endloop"    << endl;

    // facet - closing header
    file_handle << "  endfacet" << endl;

} //next i

// Closing header ------------------------------------------------------------- //
sheader << "endsolid " << solid_name;
header = sheader.str();
header = utils::string::trim(header);
file_handle << header << endl;
sheader.str("");

return(err); };

// -------------------------------------------------------------------------- //
/*!
    Write solid data to binary stl file.

    \param[in,out] file_handle stream to stl file
    \param[in,out] nV number of vertices to be written to stl file.
    \param[in,out] nT number of facet to be written to stl file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int stl::writeSolidBINARY(
    ofstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<vector<double> >             &V,
    vector<vector<double> >             &N,
    vector<vector<int> >                &T,
    string                               solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int       int2byte;
unsigned long int  longint4byte;
float              float4byte;
unsigned int      *int2byte_pointer = &int2byte;
unsigned long int *longint4byte_pointer = &longint4byte;
float             *float4byte_pointer = &float4byte;

// Counters
int                n, m, i, j, k;

BITPIT_UNUSED(solid_name);

// ============================================================================ //
// CHECK SOLID DEFINITION                                                       //
// ============================================================================ //
if (V.size() < static_cast<unsigned int>(nV))  { return(2); }
if (N.size() < static_cast<unsigned int>(nT))  { return(2); }
if (T.size() < static_cast<unsigned int>(nT))  { return(2); }

// ============================================================================ //
// CHECK STREAM STATUS                                                          //
// ============================================================================ //
if (!file_handle.good()) { return(1); }

// ============================================================================ //
// EXPORT DATA                                                                  //
// ============================================================================ //

// File header ---------------------------------------------------------------- //

// title
float4byte = (float)0.0;
for (i = 0; i < 20; i++) {
    file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
} //next i

// number of elements
longint4byte = (unsigned long int) nT;
file_handle.write(reinterpret_cast<char*>(longint4byte_pointer), 4);

// Loop over triangles -------------------------------------------------------- //
for (i = 0; i < nT; i++) {

    // Normals
    n = N[i].size();
    for (j = 0; j < n; j++) {
        float4byte = (float) N[i][j];
        file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
    } //next j

    // Vertex
    n = T[i].size();
    for (j = 0; j < n; j++) {
        m = V[T[i][j]].size();
        for (k = 0; k < m; k++) {
            float4byte = (float) V[T[i][j]][k];
            file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
        } //next k
    } //next j

    // Block closing header
    int2byte = 0;
    file_handle.write(reinterpret_cast<char*>(int2byte_pointer), 2);
} //next i

return(0); };

// -------------------------------------------------------------------------- //
/*!
    Write solid data to binary stl file. Overloading of writeSolidBINARY()
    for container vector<array<double, 3> >

    \param[in,out] file_handle stream to stl file
    \param[in,out] nV number of vertices to be written to stl file.
    \param[in,out] nT number of facet to be written to stl file.
    \param[in,out] V vertex coordinate list
    \param[in,out] N facet's normals
    \param[in,out] T facet->vertex connectivity
    \param[in] solid_name solid name

    \result returns an error flag for I/O errors:
        err = 0: no error(s) encountered
        err = 1: failed to write data to output stream
*/
unsigned int stl::writeSolidBINARY(
    ofstream                            &file_handle,
    int                                 &nV,
    int                                 &nT,
    vector<array<double,3> >            &V,
    vector<array<double,3> >            &N,
    vector<array<int,3> >               &T,
    string                               solid_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int       int2byte;
unsigned long int  longint4byte;
float              float4byte;
unsigned int      *int2byte_pointer = &int2byte;
unsigned long int *longint4byte_pointer = &longint4byte;
float             *float4byte_pointer = &float4byte;

// Counters
int                n, m, i, j, k;

BITPIT_UNUSED(solid_name);

// ============================================================================ //
// CHECK SOLID DEFINITION                                                       //
// ============================================================================ //
if (V.size() < static_cast<unsigned int>(nV))  { return(2); }
if (N.size() < static_cast<unsigned int>(nT))  { return(2); }
if (T.size() < static_cast<unsigned int>(nT))  { return(2); }

// ============================================================================ //
// CHECK STREAM STATUS                                                          //
// ============================================================================ //
if (!file_handle.good()) { return(1); }

// ============================================================================ //
// EXPORT DATA                                                                  //
// ============================================================================ //

// File header ---------------------------------------------------------------- //

// title
float4byte = (float)0.0;
for (i = 0; i < 20; i++) {
    file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
} //next i

// number of elements
longint4byte = (unsigned long int) nT;
file_handle.write(reinterpret_cast<char*>(longint4byte_pointer), 4);

// Loop over triangles -------------------------------------------------------- //
for (i = 0; i < nT; i++) {

    // Normals
    n = N[i].size();
    for (j = 0; j < n; j++) {
        float4byte = (float) N[i][j];
        file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
    } //next j

    // Vertex
    n = T[i].size();
    for (j = 0; j < n; j++) {
        m = V[T[i][j]].size();
        for (k = 0; k < m; k++) {
            float4byte = (float) V[T[i][j]][k];
            file_handle.write(reinterpret_cast<char*>(float4byte_pointer), 4);
        } //next k
    } //next j

    // Block closing header
    int2byte = 0;
    file_handle.write(reinterpret_cast<char*>(int2byte_pointer), 2);
} //next i

return(0); };

}
