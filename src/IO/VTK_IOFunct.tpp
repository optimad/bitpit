      // =================================================================================== //
      //                    VTK IO FUNCTIONS - TEMPLATES                                     //
      //                                                                                     //
      // Templeta functions to export data in vtk format.                                    //
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
      // TEMPLATES                                                                           //
      // =================================================================================== //

      // ---------------------------------------------------------------------------------- //
      template <class T>
      void Write_DataArray(ofstream        &file_handle,
                           string           dataname,
                           int              ncomponents,
                           vector<T> const &data) {

      // =================================================================================== //
      // template <class T>                                                                  //
      // void Write_DataArray(ofstream  const &file_handle,                                  //
      //                      string    const &dataname,                                     //
      //                      int       const &ncomponents,                                  //
      //                      vector<T> const &data)                                         //
      //                                                                                     //
      // Export data in a vtk file.                                                          //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - file_handle    : ofstream, handle to the .vtk file.                               //
      // - dataname       : string, data fields name                                         //
      // - ncomponents    : int, number of components                                        //
      // - data           : vector<T>, number of data                                        //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - none                                                                              //
      // =================================================================================== //

      // =================================================================================== //
      // VARIABLES DECLARATION                                                               //
      // =================================================================================== //

      // Local variables
      string data_id;
      T      dummy;

      // Counters
      int    i, j;

      // =================================================================================== //
      // DETERMINE DATA ID                                                                   //
      // =================================================================================== //
      switch (*typeid(dummy).name()) {
          case 'a': data_id = "Int8";    break;
          case 'c': data_id = "Int8";    break;
          case 'd': data_id = "Float64"; break;
          case 'e': data_id = "Float64"; break;
          case 'f': data_id = "Float32"; break;
          case 'i': data_id = "Int32";   break;
          case 'j': data_id = "";        break;
          case 'l': data_id = "Int32";   break;
          case 'm': data_id = "";        break;
          case 's': data_id = "Int8";    break;
      }


      // =================================================================================== //
      // EXPORT DATA                                                                         //
      // =================================================================================== //

      // Data header
      file_handle << "<DataArray ";
      file_handle << "type=\"" << data_id << "\" ";
      file_handle << "Name=\"" << trim(dataname) << "\" ";
      file_handle << "NumberOfComponents=\"" << ncomponents << "\" ";
      file_handle << "format=\"ascii\">" << endl;

      // Write data
      i = 0;
      while (i < data.size()) {
          for (j = 0; j < ncomponents; j++) {
              file_handle << setprecision(8) << data[i] << " ";
              i++;
          } //next j
          file_handle << endl;
      } //next i

      // Closing header
      file_handle << "</DataArray>" << endl;

      return; }
