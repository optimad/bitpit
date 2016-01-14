

#ifndef __CLASS_VTK__HH__
#define __CLASS_VTK__HH__

#include <vector>
#include <array>

#include "Generic_IO.hpp"
#include "Class_FH.hpp"



uint8_t SizeOfType( std::string type ) ;


class VTK{

    protected:

        class Field_C{

            //members
            protected:
                std::string              name;                      /**< name of the field */
                uint8_t                  components;                /**< nr of components of field, options[ 1, 3 ] */
                std::string              type;                      /**< type of data, options [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ] */
                std::string              location;                  /**< cell or point data, [ Cell, Point ] */
                std::string              codification ;             /**< Type of codification [ascii, appended] */
                uint64_t                 nr_elements;               /**< nr of cells or points */
                uint64_t                 offset;                    /**< offset in the appended section */
                std::fstream::pos_type   position;                  /**< position in file */

                //methods
            public:
                Field_C();
                Field_C( std::string , uint8_t , std::string );
                Field_C( std::string , uint8_t , std::string , std::string );
                Field_C( std::string , uint8_t , std::string , std::string , std::string , uint64_t );
                ~Field_C();

                std::string              GetName();
                std::string              GetType();
                std::string              GetLocation();
                std::string              GetCodification();
                uint8_t                  GetComponents();
                uint64_t                 GetElements();
                uint64_t                 GetSize();
                uint64_t                 GetOffset();
                uint64_t                 GetNbytes();
                std::fstream::pos_type   GetPosition(); 

                void                     SetName( std::string ) ;
                void                     SetType( std::string ) ;
                void                     SetLocation( std::string ) ;
                void                     SetCodification( std::string ) ;
                void                     SetComponents( uint8_t ) ;
                void                     SetElements( uint64_t ) ;
                void                     SetOffset( uint64_t ) ;
                void                     SetPosition( std::fstream::pos_type ) ;

        };

        // members ---------------------------------------------------------------------- //
    protected:
        FileHandler_C                   fh ;                        /**< File_Handler for Input and Output */
        uint64_t                        nr_points ;                 /**< Number of vertices */
        uint64_t                        nr_cells  ;                 /**< Number of Cells */
        int                             nr_procs  ;                 /**< Number of parallel processes  */
        int                             my_proc   ;                 /**< My process id */

        std::string                     HeaderType ;                /**< UInt32 or UInt64_t */

        std::vector< Field_C >          geometry ;                  /**< Geometry fields */
        std::string                     GeomCodex ;                 /**< Geometry codex */

        unsigned                        nr_data ;                   /**< Nr of data fields */
        std::vector< Field_C >          data ;                      /**< Data fields */
        std::string                     DataCodex ;                 /**< Data codex */

        // methods ----------------------------------------------------------------------- //
    public:
        VTK( );
        VTK( std::string dir_, std::string name_ );
        virtual ~VTK( );

        void                            SetHeaderType( std::string );
        std::string                     GetHeaderType(  );

        void                            SetNames( std::string , std::string ) ;
        void                            SetCounter( int c_=0 ) ;
        void                            SetParallel( int , int ) ;

        void                            SetCodex( std::string );
        void                            SetGeomCodex( std::string );
        void                            SetDataCodex( std::string );

        Field_C*                        AddData( std::string , int , std::string , std::string ) ;
        Field_C*                        AddData( std::string , int , std::string , std::string , std::string ) ;
        void                            RemoveData( std::string ) ;

        void                            Read() ;

        virtual void                    ReadMetaData() = 0 ; 
        void                            ReadData() ;

        void                            Write()  ;
        virtual void                    WriteMetaData() = 0 ;
        void                            WriteData() ;

        virtual void                    WriteCollection() = 0 ;

    protected:
        //For Writing
        void                            WriteDataHeader( std::fstream &, bool parallel=false ) ;
        void                            WriteDataArray( std::fstream &, Field_C &) ;
        void                            WritePDataArray( std::fstream &, Field_C &) ;

        virtual void                    Flush( std::fstream &, std::string , std::string ) =0 ;

        //For Reading
        void                            ReadDataHeader( std::fstream &) ;
        bool                            ReadDataArray( std::fstream &, Field_C &);

        virtual  void                   Absorb( std::fstream &, std::string , std::string ) =0 ;

        //General Purpose
        bool                            GetFieldByName( const std::string &, VTK::Field_C*& ) ;
        void                            CalcAppendedOffsets() ;

        bool                            StringToDataArray( std::string &, Field_C &) ;
        void                            DataArrayToString( std::string &, Field_C &) ;
        void                            PDataArrayToString( std::string &, Field_C &) ;

        template<class T> 
        std::string                     WhichType( T ) ;

        template<class T> 
        std::string                     WhichType( std::vector<T> ) ;

        template<class T, size_t d>
        std::string                     WhichType( std::array<T,d> ) ;

};

template <class Derived >
class VTK_UnstructuredGrid : public VTK{

    protected:
        uint64_t                        nconnectivity ;             /**< size of the connectivity information */

    protected:
        VTK_UnstructuredGrid();
        VTK_UnstructuredGrid( std::string , std::string ) ;
        ~VTK_UnstructuredGrid();

        void                            WriteCollection() ;  

        void                            Flush(  std::fstream &, std::string , std::string ) ; //CRTP
        void                            Absorb( std::fstream &, std::string , std::string ) ; //CRTP
        uint64_t                        CalcSizeConnectivity( ) ;

    public:
        void                            ReadMetaData() ;
        void                            WriteMetaData() ;

        void                            SetDimensions( uint64_t , uint64_t , uint64_t ) ;
        void                            SetGeomTypes( std::string , std::string , std::string , std::string ) ;

        uint8_t                         NumberOfElements( uint8_t t ) ;
        uint64_t                        GetNConnectivity( ) ; 

};

template <class Derived>
class VTK_RectilinearGrid : public VTK{

    typedef std::array<std::array<int,2>,2> extension2D_t ;         /**< typedef to describe min and max indices in 2D of restilinear grid */
    typedef std::array<std::array<int,2>,3> extension3D_t ;         /**< typedef to describe min and max indices in 3D of restilinear grid */

    protected:
        int                             dimensions ;                /**< dimensions of the grid [2/3] */
        extension3D_t                   local_index ;               /**< min and max indices of local grid */
        extension3D_t                   global_index ;              /**< min and max indices of global grid */
        std::vector<extension3D_t>      proc_index ;                /**< global indices of each processors */

    protected:
        VTK_RectilinearGrid();
        VTK_RectilinearGrid( std::string , std::string  );
        VTK_RectilinearGrid( std::string , std::string , std::string, int, int, int, int, int, int );
        VTK_RectilinearGrid( std::string , std::string , std::string, int, int, int );
        VTK_RectilinearGrid( std::string , std::string , std::string, int, int, int, int );
        VTK_RectilinearGrid( std::string , std::string , std::string, int, int );
        ~VTK_RectilinearGrid();

        void                            WriteCollection() ;  

        void                            Flush(  std::fstream &, std::string , std::string ) ; //CRTP
        void                            Absorb( std::fstream &, std::string , std::string ) ; //CRTP

    public:
        void                            ReadMetaData() ;
        void                            WriteMetaData() ;

        void                            SetDimensions( int, int, int, int, int, int ) ;
        void                            SetDimensions( int, int, int ) ;
        void                            SetDimensions( int, int, int, int ) ;
        void                            SetDimensions( int, int ) ;

        void                            SetGlobalDimensions( int, int, int ) ;
        void                            SetGlobalDimensions( int, int ) ;

        void                            SetGeomTypes( std::string ) ;

        void                            SetGlobalIndex( std::vector<extension3D_t> ) ;
        void                            SetGlobalIndex( std::vector<extension2D_t> ) ;

};

#include"Class_VTK.tpp"
#include"Class_VTK_Unstructured.tpp"
#include"Class_VTK_Rectilinear.tpp"


#endif
