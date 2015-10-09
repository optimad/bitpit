
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


uint8_t SizeOfType( string type ) ;


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
        uint8_t                  components;    // nr of components of field, options[ 1, 3 ]
        string                   type;          // type of data, options [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
        string                   location;      // cell or point data, [ Cell, Point ]
        string                   codification ; // Type of codification [ascii, appended]
        uint64_t                 nr_elements;   // nr of cells or points
        uint64_t                 offset;        // offset in the appended section
        fstream::pos_type        position;      // position in file
    
      //methods
      public:
        Field_C();
        Field_C( string name_, uint8_t comp_, string loc_ );
        Field_C( string name_, uint8_t comp_, string type_, string loc_ );
        Field_C( string name_, uint8_t comp_, string type_, string loc_, string cod_, uint64_t nr_elements_);
       ~Field_C();
    
        string                   GetName();
        string                   GetType();
        string                   GetLocation();
        string                   GetCodification();
        uint8_t                  GetComponents();
        uint64_t                 GetElements();
        uint64_t                 GetSize();
        uint64_t                 GetOffset();
        uint64_t                 GetNbytes();
        fstream::pos_type        GetPosition(); 

        void                     SetName( string name_ ) ;
        void                     SetType( string type_ ) ;
        void                     SetLocation( string loc_ ) ;
        void                     SetCodification( string cod_ ) ;
        void                     SetComponents( uint8_t comp_ ) ;
        void                     SetElements( uint64_t elem_ ) ;
        void                     SetOffset( uint64_t offs_ ) ;
        void                     SetPosition( fstream::pos_type pos_ ) ;
    
    };
     
    // members ---------------------------------------------------------------------- //
    protected:
      FileHandler_C        fh ;                     // File_Handler for Input&Output
      uint64_t             nr_points ;              // Number of vertices
      uint64_t             nr_cells  ;              // Number of Cells
      int                  nr_procs  ;              // Number of parallel processes 
      int                  my_proc   ;              // My process id

      string               HeaderType ;            // UInt32 or UInt64_t

      vector< Field_C >    geometry ;               // Geometry fields
      string               GeomCodex ;

      unsigned             nr_data ;                // Nr of data fields
      vector< Field_C >    data ;                   // Data fields
      string               DataCodex ;

    // methods ----------------------------------------------------------------------- //
    public:
      VTK( );
      VTK( string dir_, string name_ );
      virtual ~VTK( );

      void    SetHeaderType( string sg_ );
      string  GetHeaderType(  );

      void    SetNames( string dir_ , string name_  ) ;
      void    SetCounter( int c_ ) ;
      void    SetParallel( int nr, int my ) ;

      void    SetCodex( string cod_ );
      void    SetGeomCodex( string cod_ );
      void    SetDataCodex( string cod_ );

      Field_C*  AddData( string name_, int comp_, string type_, string loc_ ) ;
      Field_C*  AddData( string name_, int comp_, string type_, string loc_, string cod_ ) ;
      void      RemoveData( string name_ ) ;

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
      bool    GetFieldByName( const string &name_, VTK::Field_C*& the_field ) ;
      void    CalcAppendedOffsets() ;

      bool    StringToDataArray( string &str, Field_C &data_ ) ;
      void    DataArrayToString( string &str, Field_C &data_ ) ;
      void    PDataArrayToString( string &str, Field_C &data_ ) ;

      template<class T>
      string  WhichType( T dum_) ;

      template<class T>
      string  WhichType( vector<T> dum_) ;

      template<class T, size_t d>
      string  WhichType( array<T,d> dum_) ;

      //For Writing
      void    WriteDataHeader( fstream &str, bool parallel ) ;
      void    WriteDataArray( fstream &str, Field_C &data_ ) ;
      void    WritePDataArray( fstream &str, Field_C &data_ ) ;

      virtual
      void    Flush( fstream &str, string codex_, string name  ) =0 ; 

      //For Reading
      void    ReadDataHeader( fstream &str ) ;
      bool    ReadDataArray( fstream &str, Field_C &field_  );

      virtual
      void    Absorb( fstream &str, string codex_, string name  ) =0 ; 
       
};

//----------------------------------------------------------------------
// Derived Classes -----------------------------------------------------
//----------------------------------------------------------------------
template <class Derived >
class VTK_UnstructuredGrid : public VTK{

  protected:
    uint64_t  nconnectivity ;

  protected:
    VTK_UnstructuredGrid();
    VTK_UnstructuredGrid( string dir_, string name_  ) ;
   ~VTK_UnstructuredGrid();

    void      WriteCollection() ;  

    void      Flush(  fstream &str, string codex_, string name  ) ; //CRTP
    void      Absorb( fstream &str, string codex_, string name  ) ; //CRTP
    uint64_t  CalcSizeConnectivity( ) ;

  public:
    void      ReadMetaData() ;
    void      WriteMetaData() ;

    void      WritePMetaData() ;

    void      SetDimensions( uint64_t ncells_, uint64_t npoints_, uint64_t nconn_ ) ;
    void      SetGeomTypes( string Ptype, string Otype, string Ttype, string Ctype ) ;

    uint8_t   NumberOfElements( uint8_t t ) ;
    uint64_t  GetNConnectivity( ) ; 

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

#include"Class_VTK.tpp"
#include"Class_VTK_Unstructured.tpp"
#include"Class_VTK_Rectilinear.tpp"


#endif
