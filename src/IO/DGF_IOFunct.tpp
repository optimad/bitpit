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
template< typename T, typename ... T2 >
void DGF_obj::load_vdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
) {

// ========================================================================== //
// template< typename T, typename ... T2 >                                    //
// void load_vdata(                                                           //
//     string      data_name,                                                 //
//     int        &n,                                                         //
//     vector< T >&data,                                                      //
//     T2     &... others)                                                    //
//                                                                            //
// Load vertex datasets from dgf file whose name is specified in "data_name". //
// If no dataset is found, the input data structure is unchanged.             //
// If no name is specified (i.e. data_name = ""), returns the first vertex    //
// dataset found by circular scanning of the dgf file.                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - data_name    : string, dataset name.                                     //
// - n            : int, number of data in the dataset                        //
// - data         : vector< T >, loaded data.                                 //
// - others       : T2 (optional) other datasets to be loaded                 //
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
template< typename T, typename ... T2 >
void DGF_obj::load_sdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
) {

// ========================================================================== //
// template< typename T, typename ... T2 >                                    //
// void load_sdata(                                                           //
//     string      data_name,                                                 //
//     int        &n,                                                         //
//     vector< T >&data,                                                      //
//     T2     &... others)                                                    //
//                                                                            //
// Load simplex datasets from dgf file whose name is specified in "data_name".//
// If no dataset is found, the input data structure is unchanged.             //
// If no name is specified (i.e. data_name = ""), returns the first simplex   //
// dataset found by circular scanning of the dgf file.                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - data_name    : string, dataset name.                                     //
// - n            : int, number of data in the dataset                        //
// - data         : vector< T >, loaded data.                                 //
// - others       : T2 (optional) other datasets to be loaded                 //
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
template< typename T, typename ... T2 >
void DGF_obj::append_vdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
) {

// ========================================================================== //
// template< typename T, typename ... T2 >                                    //
// void DGF_obj::append_vdata(                                                //
//     string      data_name,                                                 //
//     int        &n,                                                         //
//     vector< T >&data,                                                      //
//     T2     &... others)                                                    //
//                                                                            //
// Append vertex data set to dgf file.                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - data_name      : string, data set name                                   //
// - n              : int, number of data in the dataset                      //
// - data           : vector< T >, dataset to be exported                     //
// - others         : T2 (optional) others data set to be exported.           //
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
template< typename T, typename ... T2 >
void DGF_obj::append_sdata(
    string      data_name,
    int        &n,
    vector< T >&data,
    T2     &... others
) {

// ========================================================================== //
// template< typename T, typename ... T2 >                                    //
// void DGF_obj::append_sdata(                                                //
//     string      data_name,                                                 //
//     int        &n,                                                         //
//     vector< T >&data,                                                      //
//     T2     &... others)                                                    //
//                                                                            //
// Append simplex data set to dgf file.                                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - data_name      : string, data set name                                   //
// - n              : int, number of data in the dataset                      //
// - data           : vector< T >, dataset to be exported                     //
// - others         : T2 (optional) others data set to be exported.           //
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
template< typename T >
unsigned int Read_DGF_data(
    ifstream        &file_handle,
    int             &N,
    vector< T >     &Data
) {

// ========================================================================== //
// template< typename T >                                                     //
// unsigned int Read_DGF_data(                                                //
//     ifstream        &file_handle,                                          //
//     int             &N,                                                    //
//     vector< T >     &Data)                                                 //
//                                                                            //
// Read data in a dgf data block.                                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle    : ifstream, input stream to dgf file.                     //
// - N              : int, number of data currently hosted in Data            //
// - Data           : vector< T >, imported data                              //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err            : unsigned int, error flag:                               //
//                    err = 0    --> no errors encountered                    //
//                    err = 1    --> file is missing or is not accessible     //
// ========================================================================== //

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
template <typename T>
unsigned int Read_DGF_VERTEXDATA(
    ifstream        &file_handle,
    int             &n,
    vector< T >     &data,
    string           data_name
) {

// ========================================================================== //
// template <typename T>                                                      //
// unsigned int Read_DGF_VERTEXDATA(                                          //
//     ifstream        &file_handle,                                          //
//     int             &n,                                                    //
//     vector< T >     &data,                                                 //
//     string           data_name)                                            //
//                                                                            //
// Load vertex data from dgf file, corresponding to name specified in         //
// data_name. If data_name is not specified (i.e. data_name = "") the first   //
// dataset found by circular scanning of the dgf file is loaded. If no        //
// data set is found corresponding to data_name, the input data_structure is  //
// not modified.                                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle  : ifstream, input stream to dgf file.                       //
// - n            : int, number of data stored in data.                       // 
// - data         : vector< T >, loaded data                                  //
// - data_name    : string, dataset name                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err          : unsigned int, error flag:                                 //
//                  err = 0      --> no errors encountered                    //
//                  err = 1      --> file is missing or is not accessible     //
// ========================================================================== //

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
template <typename T>
unsigned int Read_DGF_SIMPLEXDATA(
    ifstream        &file_handle,
    int             &n,
    vector< T >     &data,
    string           data_name
) {

// ========================================================================== //
// template <typename T>                                                      //
// unsigned int Read_DGF_SIMPLEXDATA(                                         //
//     ifstream        &file_handle,                                          //
//     int             &n,                                                    //
//     vector< T >     &data,                                                 //
//     string           data_name)                                            //
//                                                                            //
// Load simplex data from dgf file, corresponding to name specified in        //
// data_name. If data_name is not specified (i.e. data_name = "") the first   //
// dataset found by circular scanning of the dgf file is loaded. If no        //
// data set is found corresponding to data_name, the input data_structure is  //
// not modified.                                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle  : ifstream, input stream to dgf file.                       //
// - n            : int, number of data stored in data.                       // 
// - data         : vector< T >, loaded data                                  //
// - data_name    : string, dataset name                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err          : unsigned int, error flag:                                 //
//                  err = 0      --> no errors encountered                    //
//                  err = 1      --> file is missing or is not accessible     //
// ========================================================================== //

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
template < typename T >
unsigned int Write_DGF_data(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data
) {

// ========================================================================== //
// template < typename T >                                                    //
// unsigned int Write_DGF_data(                                               //
//     ifstream        &file_handle,                                          //
//     int             &N,                                                    //
//     vector< T >     &Data)                                                 //
//                                                                            //
// Write dgf data into dgf file.                                              //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle   : ofstream, output stream to dgf file.                     //
// - N             : int, number of data to be exported                       //
// - Data          : vector< T >, data to be exported                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err           : unsigned int, error flag:                                //
//                   err = 0     --> no errors encountered                    //
//                   err = 1     --> file is missing or is not accessible     //
//                   err = 2     --> input data are badly defined             //
// ========================================================================== //

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
template < typename T >
unsigned int Write_DGF_VERTEXDATA(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data,
    string           Data_name
) {

// ========================================================================== //
// template < typename T >                                                    //
// unsigned int Write_DGF_VERTEXDATA(                                         //
//     ifstream        &file_handle,                                          //
//     int             &N,                                                    //
//     vector< T >     &Data,                                                 //
//     stirng           Data_name)                                            //
//                                                                            //
// Write vertex data into dgf file.                                           //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle   : ofstream, output stream to dgf file.                     //
// - N             : int, number of data to be exported                       //
// - Data          : vector< T >, data to be exported                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err           : unsigned int, error flag:                                //
//                   err = 0     --> no errors encountered                    //
//                   err = 1     --> file is missing or is not accessible     //
//                   err = 2     --> input data are badly defined             //
// ========================================================================== //

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
template < typename T >
unsigned int Write_DGF_SIMPLEXDATA(
    ofstream        &file_handle,
    int             &N,
    vector< T >     &Data,
    string           Data_name
) {

// ========================================================================== //
// template < typename T >                                                    //
// unsigned int Write_DGF_SIMPLEXDATA(                                        //
//     ifstream        &file_handle,                                          //
//     int             &N,                                                    //
//     vector< T >     &Data,                                                 //
//     stirng           Data_name)                                            //
//                                                                            //
// Write simplex data into dgf file.                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - file_handle   : ofstream, output stream to dgf file.                     //
// - N             : int, number of data to be exported                       //
// - Data          : vector< T >, data to be exported                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - err           : unsigned int, error flag:                                //
//                   err = 0     --> no errors encountered                    //
//                   err = 1     --> file is missing or is not accessible     //
//                   err = 2     --> input data are badly defined             //
// ========================================================================== //

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

