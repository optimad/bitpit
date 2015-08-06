
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
    
        string                   GetName();
        string                   GetType();
        string                   GetLocation();
        string                   GetCodification();
        int                      GetComponents();
        int                      GetElements();
        int                      GetSize();
        int                      GetOffset();
        int                      GetNbytes();
        fstream::pos_type        GetPosition(); 

        void                     SetName( string name_ ) ;
        void                     SetType( string type_ ) ;
        void                     SetLocation( string loc_ ) ;
        void                     SetCodification( string cod_ ) ;
        void                     SetComponents( int comp_ ) ;
        void                     SetElements( int elem_ ) ;
        void                     SetOffset( int offs_ ) ;
        void                     SetPosition( fstream::pos_type pos_ ) ;
    
        int                      SizeOfType( string type ) ;
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

      void    SetNames( string dir_ , string name_  ) ;
      void    SetCounter( int c_ ) ;
      void    SetParallel( int nr, int my ) ;

      void    SetGeomCodex( string cod_ );
      void    SetDataCodex( string cod_ );

      void    AddData( string name_, int comp_, string type_, string loc_, string cod_ ) ;
      void    RemoveData( string name_ ) ;

      void    Read() ;
      virtual 
      void    ReadMetaData() = 0 ;
      void    ReadData() ;


      void    Write()  ;
      virtual 
      void    WriteMetaData() = 0 ;
      void    WriteData() ;

      virtual
      void    WriteCollection() = 0 ;  

    protected:
      //General Purpose
      bool    GetFieldByName( vector<Field_C> &fields_, const string &name_, VTK::Field_C*& the_field ) ;
      void    CalcAppendedOffsets() ;

      bool    StringToDataArray( string &str, Field_C &data_ ) ;
      void    DataArrayToString( string &str, Field_C &data_ ) ;
      void    PDataArrayToString( string &str, Field_C &data_ ) ;


      //For Writing
      void    WriteDataHeader( fstream &str, bool parallel ) ;
      void    WriteDataArray( fstream &str, Field_C &data_ ) ;
      void    WritePDataArray( fstream &str, Field_C &data_ ) ;

      virtual
      void    Flush( fstream &str, string codex_, string name  ) ; //CRTP


      //For Reading
      void    ReadDataHeader( fstream &str ) ;
      bool    ReadDataArray( fstream &str, Field_C &field_  );

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

    void      WriteCollection() ;  

    void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP
    void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP

  public:
    void      ReadMetaData() ;
    void      WriteMetaData() ;

    void      WritePMetaData() ;

    void      SetDimensions( int ncells_, int npoints_, int nconn_ ) ;
    int       NumberOfElements( int t) ;

};

//------------------------------------------------------------------
template <class Derived>
class VTK_RectilinearGrid : public VTK{

  protected:
    int                    n1, n2, m1, m2, l1, l2 ;
    array<int,6>           global_index ;
    vector<array<int,6> >  proc_index ;

  protected:
    VTK_RectilinearGrid();
    VTK_RectilinearGrid( string dir_, string name_ );
    VTK_RectilinearGrid( string dir_, string name_, string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ );
   ~VTK_RectilinearGrid();

    void      WriteCollection() ;  

    void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP
    void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP

  public:
    void      ReadMetaData() ;
    void      WriteMetaData() ;

    void      SetParallelIndex( array<int,6> glo_, vector<array<int,6>> loc_ ) ;
    void      SetDimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ) ;

};

#include"Class_VTK_Unstructured.tpp"
#include"Class_VTK_Rectilinear.tpp"


#endif
