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
      // - Write_DataArray   : write data in a .vtk file.                                    //
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
      # ifndef __VTK_IOFUNCT_HH__
      # define __VTK_IOFUNCT_HH__

      // =================================================================================== //
      // INCLUDES                                                                            //
      // =================================================================================== //

      // Standard Template Library
      # include <vector>
      # include <string>
      # include <cstdlib>
      # include <fstream>
      # include <typeinfo>
      # include <iomanip>

      // Optimad C++ library
      # include "Operators.hpp"

      // =================================================================================== //
      // NAME SPACES                                                                         //
      // =================================================================================== //
      using namespace std;

      // =================================================================================== //
      // PROTOTYPES                                                                          //
      // =================================================================================== //

      // ----------------------------------------------------------------------------------- //
      void Open_vtr(ofstream &, string const &);

      // ----------------------------------------------------------------------------------- //
      void Close_vtr(ofstream &);

      // ----------------------------------------------------------------------------------- //
      void Open_vtu(ofstream &, string const &);

      // ----------------------------------------------------------------------------------- //
      void Close_vtu(ofstream &);

      // ----------------------------------------------------------------------------------- //
      void Open_Data(ofstream &, string);

      // ----------------------------------------------------------------------------------- //
      void Close_Data(ofstream &, string);

      // ----------------------------------------------------------------------------------- //
      template <class T>
      void Read_DataArray(ifstream         &,
                           vector<T> const &);

      // ----------------------------------------------------------------------------------- //
      template <class T>
      void Write_DataArray(ofstream        &,
                           string           ,
                           int              ,
                           vector<T> const &);

      // ----------------------------------------------------------------------------------- //
      void Write_vtrMeshData(ofstream               &,
                             int              const &,
                             int              const &,
                             int              const &,
                             vector< double > const &,
                             vector< double > const &,
                             vector< double > const &);

      // ----------------------------------------------------------------------------------- //
      void Write_vtuMeshData(ofstream                    &,
                             int                   const &,
                             int                   const &,
                             vector< double >      const &,
                             vector< int >         const &,
                             vector< int >         const &,
                             vector< short int >   const &);

      // =================================================================================== //
      // TEMPLATES                                                                           //
      // =================================================================================== //
      # include "VTK_IOFunct.tpp"

      # endif