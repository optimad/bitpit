// ========================================================================== //
//                 - GRID MANAGER FOR CARTESIAN MESHES -                      //
//                                                                            //
// Grid manager for cartesian meshes                                          //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR Class_UCartMesh2D                             //
// ========================================================================== //

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh2D::CellData2PointData(
    vector< T > &CellData,
    vector< T > &PointData
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh2D::CellData2PointData(                                //
//     vector< T > &CellData,                                                 //
//     vector< T > &PointData)                                                //
//                                                                            //
// Convert cell data into point data                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - CellData    : vector< T >, cell data. T can be of any copy-constructible //
//                 type s.t. summation (difference) and multiplication by a   //
//                 constant are defined.                                      //
// - PointData   : vector< T > with point data.                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //        
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            K, J;
int            ip, jp;

// Counters
int            i, j, k, m;

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
PointData.resize((nx+1)*(ny+1));

// ========================================================================== //
// CONVERT CELL DATA INTO POINT DATA                                          //
// ========================================================================== //
for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {

        // Cell index
        K = AccessCellData(i,j);

        for (k = 0; k < 2; k++) {
            for (m = 0; m < 2; m++) {

                // Point index
                ip = i + k;
                jp = j + m;
                J = AccessPointData(ip,jp);

                if ((ip == 0) || (ip == nx+1)) {
                    if ((jp == 0) || (jp == ny+1)) {
                        // Corner
                        PointData[J] = PointData[J] + CellData[K];
                    }
                    else {
                        // Face
                        PointData[J] = PointData[J] + 0.5 * CellData[K];
                    }
                }
                else {
                    if ((jp == 0) || (jp == ny+1)) {
                        // Face
                        PointData[J] = PointData[J] + 0.5 * CellData[K];
                    }
                    else {
                        // Bulk
                        PointData[J] = PointData[J] + 0.25 * CellData[K];
                    }
                }

            } //next m
        } //next k

    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh2D::CellData2PointData(
    vector< T > &CellData,
    vector< T > &PointData,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh2D::CellData2PointData(                                //
//     vector< T > &CellData,                                                 //
//     vector< T > &PointData,                                                //
//     T2      &... others)                                                   //
//                                                                            //
// Convert an arbitrary set of point cell data into point data.               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - CellData    : vector< T >, cell data. T can be of any copy-constructible //
//                 type s.t. summation (difference) and multiplication by a   //
//                 constant are defined.                                      //
// - PointData   : vector< T > with point data.                               //
// - others      : T2 (optional) other cell data set to be converted.         //
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
// CONVERT FIRST CELL DATA SET INTO POINT DATA                                //
// ========================================================================== //
CellData2PointData(CellData, PointData);

// ========================================================================== //
// CONVERT OTHER CELL DATA SET INTO POINT DATA                                //
// ========================================================================== //
CellData2PointData(others ...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh2D::PointData2CellData(
    vector< T > &PointData,
    vector< T > &CellData
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh2D::PointData2CellData(                                //
//     vector< T > &PointData,                                                //
//     vector< T > &CellData)                                                 //
//                                                                            //
// Convert point data into cell data.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - PointData   : vector< T >, point data. T can be any copy-constructible   //
//                 type s.t. summation (difference), and multiplication by a  //
//                 constant are defined.                                      //
// - CellData    : vector< T >, point data.                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int     K, J;
int     ip, jp;

// Counters
int     i, j, k, m;

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
CellData.resize(nx*ny, 0.0);

// ========================================================================== //
// CONVERT POINT DATA TO CELL DATA                                            //
// ========================================================================== //
for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
        K = AccessCellData(i,j);
        for (k = 0; k < 2; k++) {
            for (m = 0; m < 2; m++) {
                ip = i + k;
                jp = j + m;
                J = AccessPointData(ip,jp);
                CellData[K] = CellData[K] + 0.25 * PointData[J];
            } //next m
        } //next k
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh2D::PointData2CellData(
    vector< T > &PointData,
    vector< T > &CellData,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh2D::PointData2CellData(                                //
//     vector< T > &PointData,                                                //
//     vector< T > &CellData,                                                 //
//     T2      &... others)                                                   //
//                                                                            //
// Convertes and arbitrary set of point data into cell data.                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - PointData   : vector< T >, point data. T can be any copy-constructible   //
//                 type s.t. summation (difference), and multiplication by a  //
//                 constant are defined.                                      //
// - CellData    : vector< T >, point data.                                   //
// - others      : T2 (optional) other data set to be converted.              //
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
// CONVERTES FIRST POINT DATA SET INTO CELL DATA                              //
// ========================================================================== //
PointData2CellData(PointData, CellData);

// ========================================================================== //
// ITERATIVELY CONVERTES OTHER DATA SET                                       //
// ========================================================================== //
PointData2CellData(others ...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh2D::interpolateCellData(
   dvector1D    &P,
   vector< T >  &field,
   T            &value
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh2D::interpolateCellData(                               //
//     dvector1D   &P,                                                        //
//     vector< T > &field,                                                    //
//     T           &value)                                                    //
//                                                                            //
// Interpolate cell data at a given point (1st order interpolation).          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P      : dvector1D, point coordinates                                    //
// - field  : vector< T >, cell data. T can be any copy-constructible type    //
//            s.t. summation and multiplication by a constant are defined.    //
// - value  : T, interpolation result                                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int     i0, j0, ip, jp;
double  wx0, wx1, wy0, wy1;

// Counters
// none

// ========================================================================== //
// INTERPOLATE CELL DATA                                                      //
// ========================================================================== //

// Find cell index
ReturnCellID(P, i0, j0);
if (P[0] > xnode[i0]) { ip = min(i0+1, nx-1);   }
else                  { ip = max(0, i0-1);      }
if (P[1] > ynode[j0]) { jp = min(j0+1, nx-1);   }
else                  { jp = max(0, j0-1);      }

// Interpolation weights
wx1 = max(0.0, min(1.0, abs((P[0] - xnode[i0])/dx)));     wx0 = 1.0 - wx1;
wy1 = max(0.0, min(1.0, abs((P[1] - ynode[j0])/dy)));     wy0 = 1.0 - wy1;

// Interpolation
value = wx0 * wy0 * field[AccessCellData(i0,j0)]
      + wx0 * wy1 * field[AccessCellData(i0,jp)]
      + wx1 * wy0 * field[AccessCellData(ip,j0)]
      + wx1 * wy1 * field[AccessCellData(ip,jp)];

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh2D::interpolateCellData(
   dvector1D    &P,
   vector< T >  &field,
   T            &value,
   T2       &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh2D::interpolateCellData(                               //
//    dvector1D    &P,                                                        //
//    vector< T >  &field,                                                    //
//    T            &value,                                                    //
//    T2       &... others)                                                   //
//                                                                            //
// Interpolate a set of cell data at a specified point.                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P      : dvector1D, point coordinates                                    //
// - field  : vector< T >, cell data. T can be any copy-constructible type    //
//            s.t. summation and multiplication by a constant are defined.    //
// - value  : T, interpolation result                                         //
// - others : T2 (optional) other cell data sets to be used for interpolation //
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
// INTERPOLATION USING THE FIRST CELL DATA SET                                //
// ========================================================================== //
interpolateCellData(P, field, value);

// ========================================================================== //
// INTERPOLATION USING THE THE OTHER DATASETS                                 //
// ========================================================================== //
interpolateCellData(P, others...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh2D::interpolatePointData(
    dvector1D   &P,
    vector<T>   &field,
    T           &value
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh2D::interpolatePointData(                              //
//     dvector1D   &P,                                                        //
//     vector<T>   &field,                                                    //
//     T           &value)                                                    //
//                                                                            //
// Interpolate point data at a given point.                                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : dvector1D, point x, y coordinates                                //
// - field : vector<T>, point data to be interpolated. T can be any copy-     //
//           constructible type s.t. summation (difference) and               //
//           multiplication by a scalar must be defined.                      //
// - value : T, interpolation result                                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         i0, j0, ip, jp;
double      wx0, wx1;
double      wy0, wy1;

// Counters
// none

// ========================================================================== //
// INTERPOLATE POINT DATA                                                     //
// ========================================================================== //

// Closest grid point
i0 = max(0, min(nx, (int) floor((P[0] - xlim[0])/dx)));
j0 = max(0, min(ny, (int) floor((P[1] - ylim[0])/dy)));
if (P[0] < xedge[i0]) { ip = max(0, i0-1); }
else                  { ip = min(nx, i0+1); }
if (P[1] < yedge[j0]) { jp = max(0, j0-1); }
else                  { jp = min(ny, j0+1); }

// Interpolation weights
wx1 = max(0.0, min(1.0, abs((P[0] - xedge[i0])/dx)));     wx0 = 1.0 - wx1;
wy1 = max(0.0, min(1.0, abs((P[1] - yedge[j0])/dy)));     wy0 = 1.0 - wy1;

// Interpolated value
value = wx0 * wy0 * field[AccessPointData(i0,j0)]
      + wx0 * wy1 * field[AccessPointData(i0,jp)]
      + wx1 * wy0 * field[AccessPointData(ip,j0)]
      + wx1 * wy1 * field[AccessPointData(ip,jp)];

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh2D::interpolatePointData(
    dvector1D   &P,
    vector<T>   &field,
    T           &value,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh2D::interpolatePointData(                              //
//     dvector1D   &P,                                                        //
//     vector<T>   &field,                                                    //
//     T           &value,                                                    //
//     T2      &... others)                                                   //
//                                                                            //
// Interpolate point datasets at a given point.                               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P        : dvector1D, point x, y coordinates                             //
// - field    : vector< T >, point dataset to be interpolated. T can be any   //
//              copy-constructible type, s.t. summation (difference) and      //
//              multiplication by a scalare are defined.                      //
// - value    : T, interpolation result.                                      //
// - others   : T2 (optional), interpolation results                          //
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
// INTERPOLATE THE FIRST VERTEX DATASET                                       //
// ========================================================================== //
interpolatePointData(P, field, value);

// ========================================================================== //
// ITERATIVELY INTERPOLATE THE OTHER VERTEX DATASETS                          //
// ========================================================================== //
interpolatePointData(P, others ...);

return; }

        // -------------------------------------------------------------------------------------- //
        template <class T>
        void Class_UCartMesh2D::Export_CellData_vtr(string      filename,
                                                    string      dataname,
                                                    vector< T > &CellData) {

        // ====================================================================================== //
        // template <class T>                                                                     //
        // void Class_UCartMesh2D::Export_CellData_vtr(string       filename,                     //
        //                                             string       dataname,                     //
        //                                             vector< T > &CellData)                     //
        //                                                                                        //
        // Export cell data into a .vtr file.                                                     //
        // ====================================================================================== //
        // INPUT                                                                                  //
        // ====================================================================================== //
        // - filename    : string, .vtr file name                                                 //
        // - dataname    : string, data name                                                      //
        // - CellData    : [nx*ny-by-1] vector< T > with cell data                                //
        // ====================================================================================== //
        // OUTPUT                                                                                 //
        // ====================================================================================== //
        // - none                                                                                 //
        // ====================================================================================== //

        // ====================================================================================== //
        // VARIABLES DECLARATION                                                                  //
        // ====================================================================================== //

        // Local variables
        ofstream        file_handle;
        dvector1D       z(1, 0.0);
        vector< T >     out(nx*ny, (T) 0.0);

        // Counters
        int             i, j, k;

        // ====================================================================================== //
        // PREPARE DATA                                                                           //
        // ====================================================================================== //
        k = 0;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                out[k] = CellData[AccessCellData(i,j)];
                k++;
            } //next i
        } //next j

        // ====================================================================================== //
        // EXPORT CELL DATA                                                                       //
        // ====================================================================================== //

        // Open file
        Open_vtr(file_handle, trim(filename));

        // Write mesh data
        Write_vtrMeshData(file_handle, nx, ny, 0, xedge, yedge, z);

        // Write cell data
        Open_Data(file_handle, "CellData");
        Write_DataArray(file_handle, trim(dataname), 1, out);
        Close_Data(file_handle, "CellData");

        // Close file
        Close_vtr(file_handle);

        return; };

        // -------------------------------------------------------------------------------------- //
        template <class T>
        void Class_UCartMesh2D::Export_PointData_vtr(string      filename,
                                                     string      dataname,
                                                     vector< T > &PointData) {

        // ====================================================================================== //
        // template <class T>                                                                     //
        // void Class_UCartMesh2D::Export_PointData_vtr(string       filename,                    //
        //                                              string       dataname,                    //
        //                                              vector< T > &PointData)                   //
        //                                                                                        //
        // Export point data into a .vtr file.                                                    //
        // ====================================================================================== //
        // INPUT                                                                                  //
        // ====================================================================================== //
        // - filename    : string, .vtr file name                                                 //
        // - dataname    : string, data name                                                      //
        // - PointData   : [(nx+1)*(ny+1)-by-1] vector< T > with point data                       //
        // ====================================================================================== //
        // OUTPUT                                                                                 //
        // ====================================================================================== //
        // - none                                                                                 //
        // ====================================================================================== //

        // ====================================================================================== //
        // VARIABLES DECLARATION                                                                  //
        // ====================================================================================== //

        // Local variables
        ofstream        file_handle;
        dvector1D       z(1, 0.0);
        vector< T >     out((nx+1)*(ny+1), (T) 0.0);

        // Counters
        int             i, j, k;

        // ====================================================================================== //
        // PREPARE DATA                                                                           //
        // ====================================================================================== //
        k = 0;
        for (j = 0; j < ny+1; j++) {
            for (i = 0; i < nx+1; i++) {
                out[k] = PointData[AccessPointData(i,j)];
                k++;
            } //next i
        } //next j

        // ====================================================================================== //
        // EXPORT CELL DATA                                                                       //
        // ====================================================================================== //

        // Open file
        Open_vtr(file_handle, trim(filename));

        // Write mesh data
        Write_vtrMeshData(file_handle, nx, ny, 0, xedge, yedge, z);

        // Write cell data
        Open_Data(file_handle, "PointData");
        Write_DataArray(file_handle, trim(dataname), 1, out);
        Close_Data(file_handle, "PointData");

        // Close file
        Close_vtr(file_handle);

        return; };

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS FOR Class_UCartMesh3D                             //
// ========================================================================== //
 
 // -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh3D::CellData2PointData(
    vector< T > &CellData,
    vector< T > &PointData
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh3D::CellData2PointData(                                //
//     vector< T > &CellData,                                                 //
//     vector< T > &PointData)                                                //
//                                                                            //
// Convert cell data into point data                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - CellData    : vector< T >, cell data. T can be of any copy-constructible //
//                 type s.t. summation (difference) and multiplication by a   //
//                 constant are defined.                                      //
// - PointData   : vector< T > with point data.                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //        
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int            K, J;
int            ip, jp, kp;

// Counters
int            i, j, k, l, m, n;

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
PointData.resize((nx+1)*(ny+1));

// ========================================================================== //
// CONVERT CELL DATA INTO POINT DATA                                          //
// ========================================================================== //
for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {

            // Cell index
            K = AccessCellData(i, j, k);
            for (l = 0; l < 2; l++) {
                for (m = 0; m < 2; m++) {
                    for (n = 0; n < 2; n++) {

                        // Point index
                        ip = i + l;
                        jp = j + m;
                        kp = j + n;
                        J = AccessPointData(ip,jp,kp);

                        if ((ip == 0) || (ip == nx+1)) {
                            if ((jp == 0) || (jp == ny+1)) {
                                if ((kp == 0) || (kp == nz+1)) {
                                    // Corner
                                    PointData[J] = PointData[J] + CellData[K];
                                }
                                else {
                                    // Edge
                                    PointData[J] = PointData[J] + 0.5 * CellData[K];
                                }
                            }
                            else {
                                if ((kp == 0) || (kp == nz+1)) {
                                    //Edge
                                    PointData[J] = PointData[J] + 0.5 * CellData[K];
                                }
                                else {
                                    // Face
                                    PointData[J] = PointData[J] + 0.25 * CellData[K];
                                }
                            }
                        }
                        else {
                            if ((jp == 0) || (jp == ny+1)) {
                                if ((kp == 0) || (kp == nz+1)) {
                                    // Edge
                                    PointData[J] = PointData[J] + 0.5 * CellData[K];
                                }
                                else {
                                    // Face
                                    PointData[J] = PointData[J] + 0.25 * CellData[K];
                                }
                            }
                            else {
                                if ((kp == 0) || (kp == nz+1)) {
                                    // Face
                                    PointData[J] = PointData[J] + 0.25 * CellData[K];
                                }
                                else {
                                    // Bulk
                                    PointData[J] = PointData[J] + 0.125 * CellData[K];
                                }
                            }
                        }
                    } //next n
                } //next m
            } //next l

        } //next k
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh3D::CellData2PointData(
    vector< T > &CellData,
    vector< T > &PointData,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh3D::CellData2PointData(                                //
//     vector< T > &CellData,                                                 //
//     vector< T > &PointData,                                                //
//     T2      &... others)                                                   //
//                                                                            //
// Convert an arbitrary set of point cell data into point data.               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - CellData    : vector< T >, cell data. T can be of any copy-constructible //
//                 type s.t. summation (difference) and multiplication by a   //
//                 constant are defined.                                      //
// - PointData   : vector< T > with point data.                               //
// - others      : T2 (optional) other cell data set to be converted.         //
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
// CONVERT FIRST CELL DATA SET INTO POINT DATA                                //
// ========================================================================== //
CellData2PointData(CellData, PointData);

// ========================================================================== //
// CONVERT OTHER CELL DATA SET INTO POINT DATA                                //
// ========================================================================== //
CellData2PointData(others ...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh3D::PointData2CellData(
    vector< T > &PointData,
    vector< T > &CellData
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh3D::PointData2CellData(                                //
//     vector< T > &PointData,                                                //
//     vector< T > &CellData)                                                 //
//                                                                            //
// Convert point data into cell data.                                         //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - PointData   : vector< T >, point data. T can be any copy-constructible   //
//                 type s.t. summation (difference), and multiplication by a  //
//                 constant are defined.                                      //
// - CellData    : vector< T >, point data.                                   //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int     K, J;
int     ip, jp, kp;

// Counters
int     i, j, k, l, m, n;

// ========================================================================== //
// RESIZE OUTPUT VARIABLES                                                    //
// ========================================================================== //
CellData.resize(nx*ny*nz, 0.0);

// ========================================================================== //
// CONVERT POINT DATA TO CELL DATA                                            //
// ========================================================================== //
for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
            K = AccessCellData(i,j,k);
            for (l = 0; l < 2; l++) {
                for (m = 0; m < 2; m++) {
                    for (n = 0; n < 2; n++) {
                        ip = i + l;
                        jp = j + m;
                        kp = k + n;
                        J = AccessPointData(ip,jp,kp);
                        CellData[K] = CellData[K] + 0.25 * PointData[J];
                    } //next n
                } //next m
            } //next l
        } //next k
    } //next j
} //next i

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh3D::PointData2CellData(
    vector< T > &PointData,
    vector< T > &CellData,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh3D::PointData2CellData(                                //
//     vector< T > &PointData,                                                //
//     vector< T > &CellData,                                                 //
//     T2      &... others)                                                   //
//                                                                            //
// Convertes and arbitrary set of point data into cell data.                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - PointData   : vector< T >, point data. T can be any copy-constructible   //
//                 type s.t. summation (difference), and multiplication by a  //
//                 constant are defined.                                      //
// - CellData    : vector< T >, point data.                                   //
// - others      : T2 (optional) other data set to be converted.              //
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
// CONVERTES FIRST POINT DATA SET INTO CELL DATA                              //
// ========================================================================== //
PointData2CellData(PointData, CellData);

// ========================================================================== //
// ITERATIVELY CONVERTES OTHER DATA SET                                       //
// ========================================================================== //
PointData2CellData(others ...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh3D::interpolateCellData(
   dvector1D    &P,
   vector< T >  &field,
   T            &value
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh3D::interpolateCellData(                               //
//     dvector1D   &P,                                                        //
//     vector< T > &field,                                                    //
//     T           &value)                                                    //
//                                                                            //
// Interpolate cell data at a given point (1st order interpolation).          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P      : dvector1D, point coordinates                                    //
// - field  : vector< T >, cell data. T can be any copy-constructible type    //
//            s.t. summation and multiplication by a constant are defined.    //
// - value  : T, interpolation result                                         //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int     i0, j0, k0, ip, jp, kp;
double  wx0, wx1, wy0, wy1, wz0, wz1;

// Counters
// none

// ========================================================================== //
// INTERPOLATE CELL DATA                                                      //
// ========================================================================== //

// Find cell index
ReturnCellID(P, i0, j0, k0);
if (P[0] > xnode[i0]) { ip = min(i0+1, nx-1);   }
else                  { ip = max(0, i0-1);      }
if (P[1] > ynode[j0]) { jp = min(j0+1, nx-1);   }
else                  { jp = max(0, j0-1);      }
if (P[2] > znode[k0]) { kp = min(k0+1, nz-1);   }
else                  { kp = max(0, k0-1);      }

// Interpolation weights
wx1 = max(0.0, min(1.0, abs((P[0] - xnode[i0])/dx)));     wx0 = 1.0 - wx1;
wy1 = max(0.0, min(1.0, abs((P[1] - ynode[j0])/dy)));     wy0 = 1.0 - wy1;
wz1 = max(0.0, min(1.0, abs((P[2] - znode[k0])/dz)));     wz0 = 1.0 - wz1;

// Interpolation
value = wz0 * wx0 * wy0 * field[AccessCellData(i0,j0,k0)]
      + wz0 * wx0 * wy1 * field[AccessCellData(i0,jp,k0)]
      + wz0 * wx1 * wy0 * field[AccessCellData(ip,j0,k0)]
      + wz0 * wx1 * wy1 * field[AccessCellData(ip,jp,k0)]
      + wz1 * wx0 * wy0 * field[AccessCellData(i0,j0,kp)]
      + wz1 * wx0 * wy1 * field[AccessCellData(i0,jp,kp)]
      + wz1 * wx1 * wy0 * field[AccessCellData(ip,j0,kp)]
      + wz1 * wx1 * wy1 * field[AccessCellData(ip,jp,kp)];

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh3D::interpolateCellData(
   dvector1D    &P,
   vector< T >  &field,
   T            &value,
   T2       &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh3D::interpolateCellData(                               //
//    dvector1D    &P,                                                        //
//    vector< T >  &field,                                                    //
//    T            &value,                                                    //
//    T2       &... others)                                                   //
//                                                                            //
// Interpolate a set of cell data at a specified point.                       //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P      : dvector1D, point coordinates                                    //
// - field  : vector< T >, cell data. T can be any copy-constructible type    //
//            s.t. summation and multiplication by a constant are defined.    //
// - value  : T, interpolation result                                         //
// - others : T2 (optional) other cell data sets to be used for interpolation //
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
// INTERPOLATION USING THE FIRST CELL DATA SET                                //
// ========================================================================== //
interpolateCellData(P, field, value);

// ========================================================================== //
// INTERPOLATION USING THE THE OTHER DATASETS                                 //
// ========================================================================== //
interpolateCellData(P, others...);

return; }

// -------------------------------------------------------------------------- //
template < class T >
void Class_UCartMesh3D::interpolatePointData(
    dvector1D   &P,
    vector<T>   &field,
    T           &value
) {

// ========================================================================== //
// template < class T >                                                       //
// void Class_UCartMesh3D::interpolatePointData(                              //
//     dvector1D   &P,                                                        //
//     vector<T>   &field,                                                    //
//     T           &value)                                                    //
//                                                                            //
// Interpolate point data at a given point.                                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : dvector1D, point x, y, z coordinates                             //
// - field : vector<T>, point data to be interpolated. T can be any copy-     //
//           constructible type s.t. summation (difference) and               //
//           multiplication by a scalar must be defined.                      //
// - value : T, interpolation result                                          //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int         i0, j0, k0, ip, jp, kp;
double      wx0, wx1, wy0, wy1, wz0, wz1;

// Counters
// none

// ========================================================================== //
// INTERPOLATE POINT DATA                                                     //
// ========================================================================== //

// Closest grid point
i0 = max(0, min(nx, (int) floor((P[0] - xlim[0])/dx)));
j0 = max(0, min(ny, (int) floor((P[1] - ylim[0])/dy)));
k0 = max(0, min(nz, (int) floor((P[2] - zlim[0])/dz)));
if (P[0] < xedge[i0]) { ip = max(0, i0-1); }
else                  { ip = min(nx, i0+1); }
if (P[1] < yedge[j0]) { jp = max(0, j0-1); }
else                  { jp = min(ny, j0+1); }
if (P[2] < zedge[k0]) { kp = max(0, k0-1); }
else                  { kp = min(nz, k0+1); }

// Interpolation weights
wx1 = max(0.0, min(1.0, abs((P[0] - xedge[i0])/dx)));     wx0 = 1.0 - wx1;
wy1 = max(0.0, min(1.0, abs((P[1] - yedge[j0])/dy)));     wy0 = 1.0 - wy1;
wz1 = max(0.0, min(1.0, abs((P[2] - zedge[k0])/dz)));     wz0 = 1.0 - wz1;

// Interpolated value
value = wz0 * wx0 * wy0 * field[AccessPointData(i0,j0,k0)]
      + wz0 * wx0 * wy1 * field[AccessPointData(i0,jp,k0)]
      + wz0 * wx1 * wy0 * field[AccessPointData(ip,j0,k0)]
      + wz0 * wx1 * wy1 * field[AccessPointData(ip,jp,k0)]
      + wz1 * wx0 * wy0 * field[AccessPointData(i0,j0,kp)]
      + wz1 * wx0 * wy1 * field[AccessPointData(i0,jp,kp)]
      + wz1 * wx1 * wy0 * field[AccessPointData(ip,j0,kp)]
      + wz1 * wx1 * wy1 * field[AccessPointData(ip,jp,kp)];

return; };

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void Class_UCartMesh3D::interpolatePointData(
    dvector1D   &P,
    vector<T>   &field,
    T           &value,
    T2      &... others
) {

// ========================================================================== //
// template < class T, typename ... T2 >                                      //
// void Class_UCartMesh3D::interpolatePointData(                              //
//     dvector1D   &P,                                                        //
//     vector<T>   &field,                                                    //
//     T           &value,                                                    //
//     T2      &... others)                                                   //
//                                                                            //
// Interpolate point datasets at a given point.                               //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P        : dvector1D, point x, y, z coordinates                          //
// - field    : vector< T >, point dataset to be interpolated. T can be any   //
//              copy-constructible type, s.t. summation (difference) and      //
//              multiplication by a scalare are defined.                      //
// - value    : T, interpolation result.                                      //
// - others   : T2 (optional), interpolation results                          //
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
// INTERPOLATE THE FIRST VERTEX DATASET                                       //
// ========================================================================== //
interpolatePointData(P, field, value);

// ========================================================================== //
// ITERATIVELY INTERPOLATE THE OTHER VERTEX DATASETS                          //
// ========================================================================== //
interpolatePointData(P, others ...);

return; }

        // -------------------------------------------------------------------------------------- //
        template <class T>
        void Class_UCartMesh3D::Export_CellData_vtr(string      filename,
                                                    string      dataname,
                                                    vector< T > &CellData) {

        // ====================================================================================== //
        // template <class T>                                                                     //
        // void Class_UCartMesh3D::Export_CellData_vtr(string       filename,                     //
        //                                             string       dataname,                     //
        //                                             vector< T > &CellData)                     //
        //                                                                                        //
        // Export cell data into a .vtr file.                                                     //
        // ====================================================================================== //
        // INPUT                                                                                  //
        // ====================================================================================== //
        // - filename    : string, .vtr file name                                                 //
        // - dataname    : string, data name                                                      //
        // - CellData    : [nx*ny*nz-by-1] vector< T > with cell data                             //
        // ====================================================================================== //
        // OUTPUT                                                                                 //
        // ====================================================================================== //
        // - none                                                                                 //
        // ====================================================================================== //

        // ====================================================================================== //
        // VARIABLES DECLARATION                                                                  //
        // ====================================================================================== //

        // Local variables
        ofstream        file_handle;
        vector< T >     out(nx*ny*nz, (T) 0.0);

        // Counters
        int             i, j, k, m;

        // ====================================================================================== //
        // PREPARE DATA                                                                           //
        // ====================================================================================== //
        m = 0;
        for (k = 0; k < nz; k++) {
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    out[m] = CellData[AccessCellData(i,j,k)];
                    m++;
                } //next i
            } //next j
        } //next k

        // ====================================================================================== //
        // EXPORT CELL DATA                                                                       //
        // ====================================================================================== //

        // Open file
        Open_vtr(file_handle, trim(filename));

        // Write mesh data
        Write_vtrMeshData(file_handle, nx, ny, nz, xedge, yedge, zedge);

        // Write cell data
        Open_Data(file_handle, "CellData");
        Write_DataArray(file_handle, trim(dataname), 1, out);
        Close_Data(file_handle, "CellData");

        // Close file
        Close_vtr(file_handle);

        return; };

        // -------------------------------------------------------------------------------------- //
        template <class T>
        void Class_UCartMesh3D::Export_PointData_vtr(string      filename,
                                                     string      dataname,
                                                     vector< T > &PointData) {

        // ====================================================================================== //
        // template <class T>                                                                     //
        // void Class_UCartMesh3D::Export_PointData_vtr(string       filename,                    //
        //                                              string       dataname,                    //
        //                                              vector< T > &PointData)                   //
        //                                                                                        //
        // Export point data into a .vtr file.                                                    //
        // ====================================================================================== //
        // INPUT                                                                                  //
        // ====================================================================================== //
        // - filename    : string, .vtr file name                                                 //
        // - dataname    : string, data name                                                      //
        // - PointData   : [(nx+1)*(ny+1)*(nz+1)-by-1] vector< T > with point data                //
        // ====================================================================================== //
        // OUTPUT                                                                                 //
        // ====================================================================================== //
        // - none                                                                                 //
        // ====================================================================================== //

        // ====================================================================================== //
        // VARIABLES DECLARATION                                                                  //
        // ====================================================================================== //

        // Local variables
        ofstream        file_handle;
        vector< T >     out((nx+1)*(ny+1)*(nz+1), (T) 0.0);

        // Counters
        int             i, j, k, m;

        // ====================================================================================== //
        // PREPARE DATA                                                                           //
        // ====================================================================================== //
        m = 0;
        for (k = 0; k < nz+1; k++) {
            for (j = 0; j < ny+1; j++) {
                for (i = 0; i < nx+1; i++) {
                    out[m] = PointData[AccessPointData(i,j,k)];
                    m++;
                } //next i
            } //next j
        } //next k

        // ====================================================================================== //
        // EXPORT CELL DATA                                                                       //
        // ====================================================================================== //

        // Open file
        Open_vtr(file_handle, trim(filename));

        // Write mesh data
        Write_vtrMeshData(file_handle, nx, ny, nz, xedge, yedge, zedge);

        // Write cell data
        Open_Data(file_handle, "PointData");
        Write_DataArray(file_handle, trim(dataname), 1, out);
        Close_Data(file_handle, "PointData");

        // Close file
        Close_vtr(file_handle);

        return; };

