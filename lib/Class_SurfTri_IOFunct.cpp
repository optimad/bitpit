// ========================================================================== //
//                         - Class_SurfTri -                                  //
//                                                                            //
// Grid manager for unstructured meshes.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Class_SurfTri.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::Export_dgf(
    string filename
) {

// ========================================================================== //
// void Class_SurfTri::Export_dgf(                                            //
//     string filename)                                                       //
//                                                                            //
// Export mesh data in a dgf file.                                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - filename      : string, .dgf file name.                                  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
DGF_obj                    DGF(trim(filename));

// Counters
// none

// ========================================================================== //
// EXPORT MESH DATA                                                           //
// ========================================================================== //
DGF.save(nVertex, nSimplex, Vertex, Simplex);

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Import_dgf(
    string filename
) {

// ========================================================================== //
// void Class_SurfTri::Import_dgf(                                            //
//     string filename)                                                       //
//                                                                            //
// Import mesh data from a dgf file.                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - filename     : string, .dgf filename.                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATIION                                                     //
// ========================================================================== //

// Local variables
DGF_obj              DGF(trim(filename));

// Counters
// none

// ========================================================================== //
// READ MESH DATA                                                             //
// ========================================================================== //
DGF.load(nVertex, nSimplex, Vertex, Simplex);

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Export_stl(
    string filename,
    bool flag
) {

// ========================================================================== //
// void Class_SurfTri::Export_stl(                                            //
//     string filename,                                                       //
//     bool flag)                                                             //
//                                                                            //
// Export tassellation in a .stl file (triangulation only)                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - filename     : string, .stl file name                                    //
// - flag         : bool, flag for binary (flag = true) or ascii              //
//                  (flag = false) stl                                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
STL_obj                 STL(trim(filename), flag);

// Counters
// none

// ========================================================================== //
// GENERATE NORMALS IF NOT ALREADY BUILT                                      //
// ========================================================================== //
if ((Normal.size() == 0) || (Normal.size() < nSimplex)) {
    GenerateNormals();
}

// ========================================================================== //
// EXPORT MESH IN STL FILE                                                    //
// ========================================================================== //
STL.save("", nVertex, nSimplex, Vertex, Normal, Simplex);

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Import_stl(string filename, bool flag) {

// ========================================================================== //
// void Class_SurfTri::Import_stl(string filename, bool flag)                 //
//                                                                            //
// Import triangulation from a .stl file                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - filename     : string, .stl file name                                    //
// - flag         : bool, flag for binary (true), or ASCII (false) .stl file  //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
STL_obj                   STL(trim(filename), flag);

// Counters
// none

// ========================================================================== //
// READ MESH DATA FROM THE STL FILE                                           //
// ========================================================================== //
STL.load(nVertex, nSimplex, Vertex, Normal, Simplex);

return; };


// ----------------------------------------------------------------------------------- //
void Class_SurfTri::Export_vtu(
    string filename
) {

// =================================================================================== //
// void Class_SurfTri::Export_vtu(                                                     //
//     string filename)                                                                //
//                                                                                     //
// Export tasselation in a .vtu file.                                                  //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - filename      : string, .vtu filename                                             //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
int                  n, dum, connect_size;
vector<short int>    types(nSimplex, 0);
ivector1D            type_lib(6, 0);
ivector1D            offsets(nSimplex, 0);
ivector1D            connectivity;
dvector1D            vertex(3*nVertex, 0.0);
ofstream             file_handle;

// Counters
int                  i, j, k;

        // =================================================================================== //
        // PREPARE OUTPUT                                                                      //
        // =================================================================================== //

        // Element types:
        type_lib[0] = 0;
        type_lib[1] = 1;
        type_lib[2] = 3;
        type_lib[3] = 5;
        type_lib[4] = 9;
        type_lib[5] = 7;

        // Vertex coordinate list ------------------------------------------------------------ //
        k = 0;
        for (i = 0; i < nVertex; i++) {
            for (j = 0; j < 3; j++) {
                vertex[k] = Vertex[i][j];
                k++;
            } //next j
        } //next i

        // Simplex --------------------------------------------------------------------------- //

            // offsets and element types
            connect_size = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                connect_size += n;
                dum = type_lib[min(n,5)];
                types[i] = dum;
                offsets[i] = connect_size;
            } //next i

            // Resize connectivity array
            connectivity.resize(connect_size, 0);

            // Simplex-vertex connectivity
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                for (j = 0; j < n; j++) {
                    connectivity[k] = Simplex[i][j];
                    k++;
                } //next j
            } //next i

        // =================================================================================== //
        // EXPORT DATA IN A .VTU FILE                                                          //
        // =================================================================================== //

        // Open .vtu file
        Open_vtu(file_handle, trim(filename));

        // Export mesh data
        Write_vtuMeshData(file_handle, nVertex, nSimplex, vertex, connectivity, offsets, types);

        // Close .vtu file
        Close_vtu(file_handle);

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::ExportVPData_vtu(
            string       filename,
            string       dataname,
            dvector2D   &Data
        ) {

        // =================================================================================== //
        // void Class_SurfTri::ExportVPData_vtu(                                               //
        //     string       filename,                                                          //
        //     string       dataname,                                                          //
        //     dvector2D   &Data                                                               //
        //                                                                                     //
        // Export vector field evaluated at tasselation vertices in a .vtu file                //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - filename   : string, .vtr file name.                                              //
        // - dataname   : string, data label                                                   //
        // - Data       : [nVertex-by-3] dvector2D with vector field                           //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                  n, dum, connect_size;
        vector<short int>    types(nSimplex, 0);
        ivector1D            type_lib(6, 0);
        ivector1D            offsets(nSimplex, 0);
        ivector1D            connectivity;
        dvector1D            vertex(3*nVertex, 0.0);
        dvector1D            vout(3*nVertex, 0.0);
        ofstream             file_handle;

        // Counters
        int                  i, j, k;

        // =================================================================================== //
        // PREPARE OUTPUT                                                                      //
        // =================================================================================== //

        // Element types:
        type_lib[0] = 0;
        type_lib[1] = 1;
        type_lib[2] = 3;
        type_lib[3] = 5;
        type_lib[4] = 9;
        type_lib[5] = 7;

        // Vertex coordinate list ------------------------------------------------------------ //
        k = 0;
        for (i = 0; i < nVertex; i++) {
            for (j = 0; j < 3; j++) {
                vertex[k] = Vertex[i][j];
                k++;
            } //next j
        } //next i

        // Simplex --------------------------------------------------------------------------- //

            // offsets and element types
            connect_size = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                connect_size += n;
                dum = type_lib[min(n,5)];
                types[i] = dum;
                offsets[i] = connect_size;
            } //next i

            // Resize connectivity array
            connectivity.resize(connect_size, 0);

            // Simplex-vertex connectivity>
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                for (j = 0; j < n; j++) {
                    connectivity[k] = Simplex[i][j];
                    k++;
                } //next j
            } //next i

            // Vector at tasselation vertex
            k = 0;
            for (i = 0; i < nVertex; i++) {
                for (j = 0; j < Data[i].size(); j++) {
                    vout[k] = Data[i][j];
                    k++;
                } //next j
            } //next i

        // =================================================================================== //
        // EXPORT DATA IN A .VTU FILE                                                          //
        // =================================================================================== //

        // Open .vtu file
        Open_vtu(file_handle, trim(filename));

        // Export mesh data
        Write_vtuMeshData(file_handle, nVertex, nSimplex, vertex, connectivity, offsets, types);

        // Export data
        Open_Data(file_handle, "PointData");
        Write_DataArray(file_handle, dataname, Data[0].size(), vout);
        Close_Data(file_handle, "PointData");

        // Close .vtu file
        Close_vtu(file_handle);

        return; };

        // ----------------------------------------------------------------------------------- //
        void Class_SurfTri::ExportVCData_vtu(
            string       filename,
            string       dataname,
            dvector2D   &Data
        ) {

        // =================================================================================== //
        // void Class_SurfTri::ExportVCData_vtu(                                               //
        //     string       filename,                                                          //
        //     string       dataname,                                                          //
        //     dvector2D   &Data                                                               //
        //                                                                                     //
        // Export vector field evaluated at tasselation simplex in a .vtu file                 //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - filename   : string, .vtr file name.                                              //
        // - dataname   : string, data label                                                   //
        // - Data       : [nSimplex-by-3] dvector2D with vector field                          //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                  n, dum, connect_size;
        vector<short int>    types(nSimplex, 0);
        ivector1D            type_lib(6, 0);
        ivector1D            offsets(nSimplex, 0);
        ivector1D            connectivity;
        dvector1D            vertex(3*nVertex, 0.0);
        dvector1D            vout(3*nSimplex, 0.0);
        ofstream             file_handle;

        // Counters
        int                  i, j, k;

        // =================================================================================== //
        // PREPARE OUTPUT                                                                      //
        // =================================================================================== //

        // Element types:
        type_lib[0] = 0;
        type_lib[1] = 1;
        type_lib[2] = 3;
        type_lib[3] = 5;
        type_lib[4] = 9;
        type_lib[5] = 7;

        // Vertex coordinate list ------------------------------------------------------------ //
        k = 0;
        for (i = 0; i < nVertex; i++) {
            for (j = 0; j < 3; j++) {
                vertex[k] = Vertex[i][j];
                k++;
            } //next j
        } //next i

        // Simplex --------------------------------------------------------------------------- //

            // offsets and element types
            connect_size = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                connect_size += n;
                dum = type_lib[min(n,5)];
                types[i] = dum;
                offsets[i] = connect_size;
            } //next i

            // Resize connectivity array
            connectivity.resize(connect_size, 0);

            // Simplex-vertex connectivity>
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                for (j = 0; j < n; j++) {
                    connectivity[k] = Simplex[i][j];
                    k++;
                } //next j
            } //next i

            // Vector at tasselation simplex
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                for (j = 0; j < Data[i].size(); j++) {
                    vout[k] = Data[i][j];
                    k++;
                } //next j
            } //next i

        // =================================================================================== //
        // EXPORT DATA IN A .VTU FILE                                                          //
        // =================================================================================== //

        // Open .vtu file
        Open_vtu(file_handle, trim(filename));

        // Export mesh data
        Write_vtuMeshData(file_handle, nVertex, nSimplex, vertex, connectivity, offsets, types);

        // Export data
        Open_Data(file_handle, "CellData");
        Write_DataArray(file_handle, dataname, Data[0].size(), vout);
        Close_Data(file_handle, "CellData");

        // Close .vtu file
        Close_vtu(file_handle);

        return; };



