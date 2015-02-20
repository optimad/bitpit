
#ifndef __CLASS_VTK__HH__
#define __CLASS_VTK__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

#include "Generic_IO.hpp"
#include "Class_FH.hpp"


using namespace std;


// =================================================================================== //
// TYPEDEFs

// =================================================================================== //
// VTK BASE CLASS DEFINITION                                                           //
class VTK{

    protected:
    // =================================================================================== //
    // FIELD Class DEFINITIONS                                                        //
    class Field_C{
    
      //members
      protected:
        string                   name;          // name of the field
        string                   type;          // type of data, options [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
        string                   location;      // cell or point data, [ Cell, Point ]
        string                   codification ; // Type of codification [ascii, appended]
        int                      components;    // nr of components of field, options[ 1, 3 ]
        int                      nr_elements;   // nr of cells or points
        int                      size;          // total number of data element = components *nr_elements
        int                      offset;        // offset in the appended section
        int                      nbytes;        // size in bytes of field
        fstream::pos_type        position;      // position in file
    
      //methods
      public:
        Field_C();
        Field_C( string name_, int comp_, string type_, string loc_ );
        Field_C( string name_, int comp_, string type_, string loc_, string cod_, int nr_elements_);
       ~Field_C();
    
        string                   Get_Name();
        string                   Get_Type();
        string                   Get_Location();
        string                   Get_Codification();
        int                      Get_Components();
        int                      Get_Elements();
        int                      Get_Size();
        int                      Get_Offset();
        int                      Get_Nbytes();
        fstream::pos_type        Get_Position(); 

        void                     Set_Name( string name_ ) ;
        void                     Set_Type( string type_ ) ;
        void                     Set_Location( string loc_ ) ;
        void                     Set_Codification( string cod_ ) ;
        void                     Set_Components( int comp_ ) ;
        void                     Set_Elements( int elem_ ) ;
        void                     Set_Offset( int offs_ ) ;
        void                     Set_Position( fstream::pos_type pos_ ) ;
    
        int                      size_of_type( string type ) ;
    };
     
    // members ---------------------------------------------------------------------- //
    protected:
      FileHandler_C        fh ;                     // File_Handler for Input&Output
      int                  nr_points ;              // Number of vertices
      int                  nr_cells  ;              // Number of Cells
      int                  nr_procs  ;              // Number of parallel processes 
      int                  my_proc   ;              // My process id

      vector< Field_C >    geometry ;               // Geometry fields

      int                  nr_data ;                // Nr of data fields
      vector< Field_C >    data ;                   // Data fields

    // methods ----------------------------------------------------------------------- //
    public:
      VTK( );
      VTK( string dir_, string name_ );
      virtual ~VTK( );

      void    Set_Names( string dir_ , string name_  ) ;
      void    Set_Counter( int c_ ) ;
      void    Set_Parallel( int nr, int my ) ;

      void    Add_Data( string name_, int comp_, string type_, string loc_, string cod_ ) ;
      void    Remove_Data( string name_ ) ;

      virtual void    Write() = 0 ;
      virtual void    Read() = 0 ;

    protected:
      //General Purpose
      VTK::Field_C*  get_field_by_name( vector<Field_C> &fields_, const string &name_ ) ;

      //For Writing
      void    Write_Data_Header( fstream &str, bool parallel ) ;

      void    Write_DataArray( fstream &str, Field_C &data_ ) ;
      void    Write_PDataArray( fstream &str, Field_C &data_ ) ;

      void    Write_All_Appended( fstream &str ) ;

      void    Calc_Appended_Offsets() ;

      virtual
      void    Flush( fstream &str, string codex_, string name  ) ; //CRTP


      //For Reading
      void    Read_Data_Header( fstream &str ) ;
      bool    Read_DataArray( string &str, Field_C &data_ ) ;
      void    Read_FieldValues( fstream &str ) ;

      bool    Seek_and_Read( fstream &str, const string &name_, Field_C &field_  );

      virtual
      void    Absorb( fstream &str, string codex_, string name  ) ; //CRTP
       
};

//----------------------------------------------------------------------
// Derived Classes -----------------------------------------------------
//----------------------------------------------------------------------
template <class Derived >
class VTK_UnstructuredGrid : public VTK{

  public:
    typedef Derived b_t;

  protected:
    int              ncells, npoints, nconnectivity ;

  protected:
    VTK_UnstructuredGrid();
    VTK_UnstructuredGrid( string dir_, string name_, string cod_, int ncell_, int npoints_, int nconn_  ) ;
   ~VTK_UnstructuredGrid();

    void      Write_pvtu() ;

  public:
    // For Writing
    void      Write() ;
    void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP


    // For Reading
    void      Read() ;
    void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP

    // General Purpose
    void      Set_Dimensions( int ncells_, int npoints_, int nconn_ ) ;
    int       numberofelements( int t) ;

};

//------------------------------------------------------------------
template <class Derived>
class VTK_RectilinearGrid : public VTK{

  protected:
    int                    n1, n2, m1, m2, l1, l2 ;
    array<int,6>           global_index ;
    vector<array<int,6>>   proc_index ;

  protected:
    VTK_RectilinearGrid();
    VTK_RectilinearGrid( string dir_, string name_ );
    VTK_RectilinearGrid( string dir_, string name_, string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ );
   ~VTK_RectilinearGrid();

    void      Write_pvtr() ;

  public:
    // For Writing
    void      Write() ;
    void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP

    void      Set_Parallel_Index( array<int,6> glo_, vector<array<int,6>> loc_ ) ;

    // For Reading
    void      Read() ;
    void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP

    // General Purpose
    void      Set_Dimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ) ;

};

#include"Class_VTK_Unstructured.tpp"
#include"Class_VTK_Rectilinear.tpp"


#endif
