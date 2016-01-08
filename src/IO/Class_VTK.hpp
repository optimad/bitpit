
#ifndef __CLASS_VTK__HH__
#define __CLASS_VTK__HH__

#include <vector>
#include <array>

#include "Generic_IO.hpp"
#include "Class_FH.hpp"



uint8_t SizeOfType( std::string type ) ;


// =================================================================================== //
// VTK BASE CLASS DEFINITION                                                           //
class VTK{

    protected:
    // =================================================================================== //
    // FIELD Class DEFINITIONS                                                        //
    class Field_C{
    
      //members
      protected:
        std::string                   name;          // name of the field
        uint8_t                  components;    // nr of components of field, options[ 1, 3 ]
        std::string                   type;          // type of data, options [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
        std::string                   location;      // cell or point data, [ Cell, Point ]
        std::string                   codification ; // Type of codification [ascii, appended]
        uint64_t                 nr_elements;   // nr of cells or points
        uint64_t                 offset;        // offset in the appended section
        fstream::pos_type        position;      // position in file
    
      //methods
      public:
        Field_C();
        Field_C( std::string name_, uint8_t comp_, std::string loc_ );
        Field_C( std::string name_, uint8_t comp_, std::string type_, std::string loc_ );
        Field_C( std::string name_, uint8_t comp_, std::string type_, std::string loc_, std::string cod_, uint64_t nr_elements_);
       ~Field_C();
    
        std::string                   GetName();
        std::string                   GetType();
        std::string                   GetLocation();
        std::string                   GetCodification();
        uint8_t                  GetComponents();
        uint64_t                 GetElements();
        uint64_t                 GetSize();
        uint64_t                 GetOffset();
        uint64_t                 GetNbytes();
        fstream::pos_type        GetPosition(); 

        void                     SetName( std::string name_ ) ;
        void                     SetType( std::string type_ ) ;
        void                     SetLocation( std::string loc_ ) ;
        void                     SetCodification( std::string cod_ ) ;
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

      std::string               HeaderType ;            // UInt32 or UInt64_t

      std::vector< Field_C >    geometry ;               // Geometry fields
      std::string               GeomCodex ;

      unsigned             nr_data ;                // Nr of data fields
      std::vector< Field_C >    data ;                   // Data fields
      std::string               DataCodex ;

    // methods ----------------------------------------------------------------------- //
    public:
      VTK( );
      VTK( std::string dir_, std::string name_ );
      virtual ~VTK( );

      void    SetHeaderType( std::string sg_ );
      std::string  GetHeaderType(  );

      void    SetNames( std::string dir_ , std::string name_  ) ;
      void    SetCounter( int c_ ) ;
      void    SetParallel( int nr, int my ) ;

      void    SetCodex( std::string cod_ );
      void    SetGeomCodex( std::string cod_ );
      void    SetDataCodex( std::string cod_ );

      Field_C*  AddData( std::string name_, int comp_, std::string type_, std::string loc_ ) ;
      Field_C*  AddData( std::string name_, int comp_, std::string type_, std::string loc_, std::string cod_ ) ;
      void      RemoveData( std::string name_ ) ;

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
      bool    GetFieldByName( const std::string &name_, VTK::Field_C*& the_field ) ;
      void    CalcAppendedOffsets() ;

      bool    StringToDataArray( std::string &str, Field_C &data_ ) ;
      void    DataArrayToString( std::string &str, Field_C &data_ ) ;
      void    PDataArrayToString( std::string &str, Field_C &data_ ) ;

      template<class T>
      std::string  WhichType( T dum_) ;

      template<class T>
      std::string  WhichType( std::vector<T> dum_) ;

      template<class T, size_t d>
      std::string  WhichType( std::array<T,d> dum_) ;

      //For Writing
      void    WriteDataHeader( fstream &str, bool parallel ) ;
      void    WriteDataArray( fstream &str, Field_C &data_ ) ;
      void    WritePDataArray( fstream &str, Field_C &data_ ) ;

      virtual
      void    Flush( fstream &str, std::string codex_, std::string name  ) =0 ; 

      //For Reading
      void    ReadDataHeader( fstream &str ) ;
      bool    ReadDataArray( fstream &str, Field_C &field_  );

      virtual
      void    Absorb( fstream &str, std::string codex_, std::string name  ) =0 ; 
       
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
    VTK_UnstructuredGrid( std::string dir_, std::string name_  ) ;
   ~VTK_UnstructuredGrid();

    void      WriteCollection() ;  

    void      Flush(  fstream &str, std::string codex_, std::string name  ) ; //CRTP
    void      Absorb( fstream &str, std::string codex_, std::string name  ) ; //CRTP
    uint64_t  CalcSizeConnectivity( ) ;

  public:
    void      ReadMetaData() ;
    void      WriteMetaData() ;

    void      WritePMetaData() ;

    void      SetDimensions( uint64_t ncells_, uint64_t npoints_, uint64_t nconn_ ) ;
    void      SetGeomTypes( std::string Ptype, std::string Otype, std::string Ttype, std::string Ctype ) ;

    uint8_t   NumberOfElements( uint8_t t ) ;
    uint64_t  GetNConnectivity( ) ; 

};

//------------------------------------------------------------------
template <class Derived>
class VTK_RectilinearGrid : public VTK{

  protected:
    int                    n1, n2, m1, m2, l1, l2 ;
    std::array<int,6>           global_index ;
    std::vector<std::array<int,6> >  proc_index ;

  protected:
    VTK_RectilinearGrid();
    VTK_RectilinearGrid( std::string dir_, std::string name_ );
    VTK_RectilinearGrid( std::string dir_, std::string name_, std::string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ );
   ~VTK_RectilinearGrid();

    void      WriteCollection() ;  

    void      Flush(  fstream &str, std::string codex_, std::string name  ) ; //CRTP
    void      Absorb( fstream &str, std::string codex_, std::string name  ) ; //CRTP

  public:
    void      ReadMetaData() ;
    void      WriteMetaData() ;

    void      SetParallelIndex( std::array<int,6> glo_, std::vector<std::array<int,6>> loc_ ) ;
    void      SetDimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ) ;

};

#include"Class_VTK.tpp"
#include"Class_VTK_Unstructured.tpp"
#include"Class_VTK_Rectilinear.tpp"


#endif
