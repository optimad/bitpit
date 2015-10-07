      // =================================================================================== //
      //                         VTK IO FUNCTIONS                                            //
      //                                                                                     //
      // Functions to export data in vtk format.                                             //
      //                                                                                     //
      // LIST OF FUNCTIONS:                                                                  //
      //                                                                                     //
      // -------- .vtr file ---------------------------------------------------------------  //
      // - Open_vtr          : open a .vtr file.                                             //
      // - Close_vtr         : close a .vtr file.                                            //
      // - Write_vtrMeshData : write mesh data in a .vtr file.                               //
      // -------- .vtu file ---------------------------------------------------------------  //
      // - Open_vtu          : open a .vtu file.                                             //
      // - Close_vtu         : close a .vtu file.                                            //
      // - Write_vtuMeshData : write mesh data in a .vtu file.                               //
      // -------- Generic -----------------------------------------------------------------  //
      // - Open_Data         : open a data section in a .vtk file                            //
      // - Close_Data        : close a data section in a .vtk file                           //
      // =================================================================================== //
      // INFO                                                                                //
      // =================================================================================== //
      // Author     : Alessandro Alaia                                                       //
      // Company    : Optimad Engineering srl                                                //
      // Data       : Jul 17, 2013                                                           //
      // Version    : v1.0                                                                   //
      //                                                                                     //
      // All rights reserved.                                                                //
      // =================================================================================== //

      // =================================================================================== //
      // INCLUDES                                                                            //
      // =================================================================================== //
      # include "VTK_IOFunct.hpp"

      // ---------------------------------------------------------------------------------- //
      void Open_vtr(ofstream & file_handle, string const &filename) {

      // =================================================================================== //
      // Open_vtr(ofstream & file_handle, string const &filename)                            //
      //                                                                                     //
      // Open a vtr file and write file header.                                              //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle : ofstream, handle to the .vtr file.                                  //
      // - filename    : string, .vtr filename                                               //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      // none

      // Counters
      // none

      // =================================================================================== //
      // CREATE FILE HEADER                                                                  //
      // =================================================================================== //

      // Open file
      file_handle.open(filename.c_str(), ifstream::out);

      // Prepare file header
      file_handle << "<?xml version=\"1.0\"?>" << endl;
      file_handle << "<VTKFile type=\"RectilinearGrid\" "
                  << "version=\"0.1\" "
                  << "byte_order=\"BigEndian\">" << endl;

      return; };

      // ---------------------------------------------------------------------------------- //
      void Close_vtr(ofstream &file_handle) {

      // =================================================================================== //
      // void Close_vtr(ofstream const &file_handle)                                         //
      //                                                                                     //
      // Close the .vtr file and add closing header.                                         //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle   : ofstream, with opened .vtr file.                                  //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      // none

      // Counter
      // none

      // =================================================================================== //
      // CLOSE VTR FILE                                                                      //
      // =================================================================================== //

      // Add closing header
      file_handle << "    </Piece>"         << endl;
      file_handle << "  </RectilinearGrid>" << endl;
      file_handle << "</VTKFile>"           << endl;

      // Close file
      file_handle.close();

      return; };

      // ---------------------------------------------------------------------------------- //
      void Open_vtu(ofstream & file_handle, string const &filename) {

      // =================================================================================== //
      // void Open_vtu(ofstream & file_handle, string const &filename)                       //
      //                                                                                     //
      // Open a .vtu file.                                                                   //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle   : ofstream, handle to .vtu file.                                    //
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
      // none

      // Conters
      // none

      // =================================================================================== //
      // CREATE FILE HEADER                                                                  //
      // =================================================================================== //

      // Open file
      file_handle.open(filename.c_str(), ifstream::out);

      // Prepare file header
      file_handle << "<?xml version=\"1.0\"?>" << endl;
      file_handle << "<VTKFile type=\"UnstructuredGrid\" "
                  << "version=\"0.1\" "
                  << "byte_order=\"BigEndian\">" << endl;

      return; };

      // ---------------------------------------------------------------------------------- //
      void Close_vtu(ofstream &file_handle) {

      // =================================================================================== //
      // void Close_vtu(ofstream &file_handle)                                               //
      //                                                                                     //
      // Close a .vtu file.                                                                  //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle   : ofstream, handle to the .vtu file                                 //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      // none

      // Counters
      // none

      // =================================================================================== //
      // CLOSING HEADER                                                                      //
      // =================================================================================== //

      // Add closing header
      file_handle << "    </Piece>"         << endl;
      file_handle << "  </UnstructuredGrid>" << endl;
      file_handle << "</VTKFile>"           << endl;

      // Close file
      file_handle.close();

      return; };

      // ---------------------------------------------------------------------------------- //
      void Open_Data(ofstream &file_handle, string datatype) {

      // =================================================================================== //
      // void Open_Data(ofstream const &file_handle, string const &datatype)                 //
      //                                                                                     //
      // Open the data section in the .vtr file.                                             //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle     : ofstream, handle to the .vtr file.                              //
      // - datatype        : string, data type ("CellData" or "PointData")                   //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      string debug;
      bool   debug2;

      // Counters
      // none


      // =================================================================================== //
      // OPEN DATA SECTION                                                                   //
      // =================================================================================== //
      if (trim(datatype).compare("CellData") == 0) {
          file_handle << "      <CellData>" << endl;
      }
      else if (trim(datatype).compare("PointData") == 0) {
          file_handle << "      <PointData>" << endl;
      }

      return; };

      // ---------------------------------------------------------------------------------- //
      void Close_Data(ofstream &file_handle, string datatype) {

      // =================================================================================== //
      // void Close_Data(ofstream const &file_handle, string const &datatype)                //
      //                                                                                     //
      // Close the data section in the .vtr file.                                            //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle     : ofstream, handle to the .vtr file.                              //
      // - datatype        : string, data type ("CellData" or "PointData")                   //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      // none

      // Counters
      // none

      // =================================================================================== //
      // OPEN DATA SECTION                                                                   //
      // =================================================================================== //
      if (trim(datatype).compare("CellData") == 0) {
          file_handle << "      </CellData>" << endl;
      }
      else if (trim(datatype).compare("PointData") == 0) {
          file_handle << "      </PointData>" << endl;
      }

      return; };

      // ---da finire---------------------------------------------------------------------- //
      template <class T>
      void Read_DataArray(ifstream        &file_handle,
                          vector<T>       &data) {

      // =================================================================================== //
      // template <class T>                                                                  //
      // void Read_DataArray(ifstream        &file_handle,                                   //
      //                     vector<T>       &data)                                          //
      //                                                                                     //
      // Read dataarray from VTK file.                                                       //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle   : ifstream, file handle to the vtk file.                            //
      // - data          : vector< T >, vector for data storage.                             //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - data          : vector< T >, vector with stored data                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      stringstream      header;

      // Counters

      return; };

      // ---------------------------------------------------------------------------------- //
      void Write_vtrMeshData(ofstream               &file_handle,
                             int              const &nx,
                             int              const &ny,
                             int              const &nz,
                             vector< double > const &x,
                             vector< double > const &y,
                             vector< double > const &z) {

      // =================================================================================== //
      // void Write_vtrMeshData(ofstream         const &file_handle,                         //
      //                        int              const &nx,                                  //
      //                        int              const &ny,                                  //
      //                        int              const &nz,                                  //
      //                        vector< double > const &x,                                   //
      //                        vector< double > const &y,                                   //
      //                        vector< double > const &z)                                   //
      //                                                                                     //
      // Export mesh data in a .vtr file.                                                    //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle   : ofstream, handle to the .vtr file.                                //
      // - nx            : int, number of mesh cells in the x direction                      //
      // - ny            : int, number of mesh cells in the y direction                      //
      // - nz            : int, number of mesh cells in the z direction                      //
      // - x             : [nx+1] vector< double > with x coordinates of mesh vertex         //
      // - y             : [ny+1] vector< double > with y coordinates of mesh vertex         //
      // - z             : [nz+1] vector< double > with x coordinates of mesh vertex         //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables

      // Counters

      // =================================================================================== //
      // EXPORT MESH DATA                                                                    //
      // =================================================================================== //

      // Mesh header
      file_handle << "  <RectilinearGrid WholeExtent=\""
                  << " " << 1 << " " << nx+1
                  << " " << 1 << " " << ny+1
                  << " " << 1 << " " << nz+1
                  << "\">" << endl;
      file_handle << "    <Piece Extent=\""
                  << " " << 1 << " " << nx+1
                  << " " << 1 << " " << ny+1
                  << " " << 1 << " " << nz+1
                  << "\">" << endl;
      file_handle << "      <Coordinates>" << endl;

      // Export vertex coordinates
      Write_DataArray(file_handle, "X", 1, x);
      Write_DataArray(file_handle, "Y", 1, y);
      Write_DataArray(file_handle, "Z", 1, z);

      // Closing header
      file_handle << "      </Coordinates>" << endl;

      return; }

      // ---------------------------------------------------------------------------------- //
      void Write_vtuMeshData(ofstream                    &file_handle,
                             int                   const &nV,
                             int                   const &nT,
                             vector< double >      const &Vertex,
                             vector< int >         const &Connect,
                             vector< int >         const &offsets,
                             vector< short int >   const &Type) {

      // =================================================================================== //
      // void Write_vtrMeshData(ofstream                  &file_handle,                      //
      //                        int                 const &nV,                               //
      //                        int                 const &nT,                               //
      //                        vector< double >    const &Vertex,                           //
      //                        vector< int >       const &Connect,                          //
      //                        vector< int >       const &offsets,                          //
      //                        vector< short int > const &Type)                             //
      //                                                                                     //
      // Export mesh data in a .vtu file.                                                    //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle    : ofstream, handle to the .vtu file.                               //
      // - nV             : int, number of vertex                                            //
      // - nT             : int, number of simplicies                                        //
      // - Vertex         : [3 nV - by - 1] vector< double >, vertex coordinate list         //
      //                    Vertex[i] stores the x, y, z coordinates of the i-th vertex.     //
      // - Connect        : [3 nT - by - 1] vector< int >, simplex-vertex coordinate list.   //
      //                    Connect[i] stores the indexes of vertexes of the i-th simplex    //
      // - offsets        : [nT - by - 1] vector< int >, number of vertex in each simplex.   //
      //                    offsets[i] stores the number of vertexes in the i-th simplex.    //
      //                    Connect[i] stores the indexes of vertexes of the i-th simplex    //
      // - Type           : [nT - by - 1] vector< short int >, simplex type. Type[i] stores  //
      //                    a flag for each simplex indicating the simplex type.             //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //


      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      // none

      // Counters
      // none

      // =================================================================================== //
      // WRITE MESH DATA                                                                     //
      // =================================================================================== //

      // Data header
      file_handle << "  <UnstructuredGrid>" << endl;
      file_handle << "    <Piece NumberOfPoints=\" "
                  << nV
                  << "\" NumberOfCells=\" "
                  << nT
                  << "\">" << endl;

      // Vertex coordinates
      file_handle << "      <Points>" << endl;
      Write_DataArray(file_handle, "Point", 3, Vertex);
      file_handle << "      </Points>" << endl;

      // Simplicies
      file_handle << "      <Cells>" << endl;
      Write_DataArray(file_handle, "connectivity", 1, Connect);
      Write_DataArray(file_handle, "offsets", 1, offsets);
      Write_DataArray(file_handle, "types", 1, Type);
      file_handle << "      </Cells>" << endl;

      return; };