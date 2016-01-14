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
// TEMPLATE IMPLEMENTATIONS                                                   //
// ========================================================================== //

// class DGF_obj templated methods ========================================== //

// -------------------------------------------------------------------------- //
/*!
    Load vertex data from dgf file.

    \param[in] data_name name of vertex data block to be loaded.
    If no name is specified (i.e. data_name = "") the first vertex data block
    found in the dgf file is returned. If no vertex datablock is found with the
    specified name no data is loaded.
    \param[in,out] n number data entries loaded from the data set
    \param[in,out] data on output stores the data loaded from the data block
    \param[in] others other data to be loaded.
*/ 
template< typename T, typename ... T2 >
void DGF_obj::load_vdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
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
// READ DATA SET                                                              //
// ========================================================================== //
err = Read_DGF_VERTEXDATA(ifile_handle, n, data, data_name);

// ========================================================================== //
// ITERATIVELY READ OTHER DATASET                                             //
// ========================================================================== //
load_vdata(others ...);

return; };

// -------------------------------------------------------------------------- //
/*!
    Load cell data from dgf file.

    \param[in] data_name name of cell data block to be loaded.
    If no name is specified (i.e. data_name = "") the first cell data block
    found in the dgf file is returned. If no cell datablock is found with the
    specified name no data is loaded.
    \param[in,out] n number data entries loaded from the data set
    \param[in,out] data on output stores the data loaded from the data block
    \param[in] others other data to be loaded.
*/
template< typename T, typename ... T2 >
void DGF_obj::load_sdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
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
// READ DATA SET                                                              //
// ========================================================================== //
err = Read_DGF_SIMPLEXDATA(ifile_handle, n, data, data_name);

// ========================================================================== //
// ITERATIVELY READ OTHER DATASET                                             //
// ========================================================================== //
load_sdata(others ...);

return; };

// -------------------------------------------------------------------------- //
/*!
    Append vertex data to dgf file.

    \param[in] data_name name of vertex data block to be appended.
    \param[in] n number data entries to be exported to the dgf file.
    \param[in] data container storing the data to be written into the dgf file.
    \param[in] others other data to be exported.
*/ 
template< typename T, typename ... T2 >
void DGF_obj::append_vdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
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
open("app");

// ========================================================================== //
// SAVE DATA SET                                                              //
// ========================================================================== //
err = Write_DGF_VERTEXDATA(ofile_handle, n, data, data_name);

// ========================================================================== //
// RECURSIVELY SAVE OTHERS DATASETS                                           //
// ========================================================================== //
append_vdata(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Append cell data to dgf file.

    \param[in] data_name name of cell data block to be appended.
    \param[in] n number data entries to be exported to the dgf file.
    \param[in] data container storing the data to be written into the dgf file.
    \param[in] others other data to be exported.
*/ 
template< typename T, typename ... T2 >
void DGF_obj::append_sdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
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
open("app");

// ========================================================================== //
// SAVE DATA SET                                                              //
// ========================================================================== //
err = Write_DGF_SIMPLEXDATA(ofile_handle, n, data, data_name);

// ========================================================================== //
// RECURSIVELY SAVE OTHERS DATASETS                                           //
// ========================================================================== //
append_sdata(others ...);

return; }

// Input routines =========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Read data block from dgf file.

    \param[in,out] file_handle input stream from dgf file
    \param[in,out] n number data entries loaded from the data set
    \param[in,out] data on output stores the data loaded from the data block

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to load data from dgf file
*/ 
template< typename T >
unsigned int Read_DGF_data(
    ifstream        &file_handle,
    int             &N,
    vector< T >     &Data
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check;
long int            start_pos;
string              word, line;
stringstream        sline;

// Counters
int                 n = 0;

// ========================================================================== //
// CHECK INPUT STREAM                                                         //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// SCAN DATA BLOCK                                                            //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
Scan_DGF_data(file_handle, n);
file_handle.clear();
file_handle.seekg(start_pos);

// ========================================================================== //
// RESIZE INPUT DATA                                                          //
// ========================================================================== //
Data.resize(n+N);

// ========================================================================== //
// READ DATA                                                                  //
// ========================================================================== //
n = 0;
check = true;
while (!file_handle.eof() && check) {

    // Get current line
    start_pos = file_handle.tellg();
    getline(file_handle, line);
    line = trim(line);
    sline.clear();
    sline.str(line);

    // Read data
    if (sline >> word) {
         check =  ((word.compare("#") != 0)
              &&  (word.compare("VERTEX") != 0)
              &&  (word.compare("SIMPLEX") != 0)
              &&  (word.compare("VERTEXDATA") != 0)
              &&  (word.compare("SIMPLEXDATA") != 0));
        if (check) {
            sline.seekg(0);
            sline >> Data[N+n];
            n++;
        }
    }
} //next line
if (word.compare("#") != 0) {
    file_handle.clear();
    file_handle.seekg(start_pos);
}

// Update counters
N+=n;

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read vertex data from dgf file.

    \param[in,out] file_handle input stream from dgf file
    \param[in,out] n number data entries loaded from the data set
    \param[in,out] data on output stores the data loaded from the data block
    \param[in] data_name name of the vertex data set to be loaded. If no name
    is specified (i.e. data_name = "") the first vertex data set encountered
    in the dgf file is returned. If no vertex data set matching the input
    name is found, no data is loaded from the dgf file.

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to load data from dgf file
*/ 
template <typename T>
unsigned int Read_DGF_VERTEXDATA(
    ifstream        &file_handle,
    int             &n,
    vector< T >     &data,
    string           data_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool            check = false;
long int        current_pos, start_pos;
string          header, line, word;
stringstream    sline;

// Counters

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// PARAMETERS                                                                 //
// ========================================================================== //
sline << "VERTEXDATA " << data_name << endl;
header = sline.str();
header = trim(header);

// ========================================================================== //
// SCAN DGF FILE LOOKING FOR DATASET WITH SPECIFIED NAME                      //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
current_pos = -1;
while (start_pos != current_pos) {

    // Get current line
    getline(file_handle, line);
    line = trim(line);
    sline.clear();
    sline.str(line);

    // Check eof
    if (file_handle.eof()) {
        file_handle.clear();
        file_handle.seekg(0);
    }
    current_pos = file_handle.tellg();

    // Look for keyword
    if ((sline >> word) && (word.compare("VERTEXDATA") == 0)) {
        if (data_name.compare("") == 0) {
            start_pos = current_pos;
            check = true;
        }
        else {
            if (line.compare(header) == 0) {
                start_pos = current_pos;
                check = true;
            }
        }
    }

    
} //next line

// ========================================================================== //
// LOAD DATA                                                                  //
// ========================================================================== //
if (check) {
    Read_DGF_data(file_handle, n, data);
}

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Read cell data from dgf file.

    \param[in,out] file_handle input stream from dgf file
    \param[in,out] n number data entries loaded from the data set
    \param[in,out] data on output stores the data loaded from the data block
    \param[in] data_name name of the cell data set to be loaded. If no name
    is specified (i.e. data_name = "") the first cell data set encountered
    in the dgf file is returned. If no cell data set matching the input
    name is found, no data is loaded from the dgf file.

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to load data from dgf file
*/
template <typename T>
unsigned int Read_DGF_SIMPLEXDATA(
    ifstream        &file_handle,
    int             &n,
    vector< T >     &data,
    string           data_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool            check = false;
long int        current_pos, start_pos;
string          header, line, word;
stringstream    sline;

// Counters

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// PARAMETERS                                                                 //
// ========================================================================== //
sline << "SIMPLEXDATA " << data_name << endl;
header = sline.str();
header = trim(header);

// ========================================================================== //
// SCAN DGF FILE LOOKING FOR DATASET WITH SPECIFIED NAME                      //
// ========================================================================== //
file_handle.clear();
start_pos = file_handle.tellg();
current_pos = -1;
while (start_pos != current_pos) {

    // Get current line
    getline(file_handle, line);
    line = trim(line);
    sline.clear();
    sline.str(line);

    // Check eof
    if (file_handle.eof()) {
        file_handle.clear();
        file_handle.seekg(0);
    }
    current_pos = file_handle.tellg();

    // Look for keyword
    if ((sline >> word) && (word.compare("SIMPLEXDATA") == 0)) {
        if (data_name.compare("") == 0) {
            start_pos = current_pos;
            check = true;
        }
        else {
            if (line.compare(header) == 0) {
                start_pos = current_pos;
                check = true;
            }
        }
    }

    
} //next line

// ========================================================================== //
// LOAD DATA                                                                  //
// ========================================================================== //
if (check) {
    Read_DGF_data(file_handle, n, data);
}

return(0); }

// Output routines ========================================================== //

// -------------------------------------------------------------------------- //
/*!
    Write data block to dgf file.

    \param[in,out] file_handle output stream to dgf file
    \param[in,out] N number data entries to be written to file
    \param[in,out] Data container storing the data

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to write data to dgf file
*/ 
template < typename T >
unsigned int Write_DGF_data(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int             i;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// EXPORT DATA                                                                //
// ========================================================================== //
for (i = 0; i < N; i++) {
    file_handle << Data[i] << endl;
} //next i
file_handle << "#" << endl << endl;

return(0); }

// -------------------------------------------------------------------------- //
/*!
    Write vertex data to dgf file.

    \param[in,out] file_handle output stream to dgf file
    \param[in,out] N number data entries to be written
    \param[in,out] Data container storing the data
    \param[in] Data_name name of the vertex data set

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to load data from dgf file
*/ 
template < typename T >
unsigned int Write_DGF_VERTEXDATA(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data,
    string           Data_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int        err = 0;
stringstream        sheader;
string              header;

// Counters
int             i;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// EXPORT DATA                                                                //
// ========================================================================== //

// Data header -------------------------------------------------------------- //
Data_name = trim(Data_name);
sheader << "VERTEXDATA " << Data_name << endl;
header = sheader.str();
header = trim(header);
file_handle << header << endl;

// Export data -------------------------------------------------------------- //
err = Write_DGF_data(file_handle, N, Data);

return(err); };

// -------------------------------------------------------------------------- //
/*!
    Write cell data to dgf file.

    \param[in,out] file_handle output stream to dgf file
    \param[in,out] N number data entries to be written
    \param[in,out] Data container storing the data
    \param[in] Data_name name of the cell data set

    \result error flag tor I/O errors
        err = 0: no error(s) encountered
        err = 1: failed to load data from dgf file
*/ 
template < typename T >
unsigned int Write_DGF_SIMPLEXDATA(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data,
    string           Data_name
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
unsigned int        err = 0;
stringstream        sheader;
string              header;

// Counters
int             i;

// ========================================================================== //
// CHECK STREAM STATUS                                                        //
// ========================================================================== //
if (!file_handle.good()) { return(1); }

// ========================================================================== //
// EXPORT DATA                                                                //
// ========================================================================== //

// Data header -------------------------------------------------------------- //
Data_name = trim(Data_name);
sheader << "SIMPLEXDATA " << Data_name << endl;
header = sheader.str();
header = trim(header);
file_handle << header << endl;

// Export data -------------------------------------------------------------- //
err = Write_DGF_data(file_handle, N, Data);

return(err); };

