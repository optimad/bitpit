

#ifndef __VTK__HH__
#define __VTK__HH__

#include <vector>
#include <array>

//#include "BitPit_common.hpp"
#include "Generic_IO.hpp"
#include "FileHandler.hpp"


class VTKField{

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
        VTKField();
        VTKField( std::string , uint8_t , std::string );
        VTKField( std::string , uint8_t , std::string , std::string );
        VTKField( std::string , uint8_t , std::string , std::string , std::string , uint64_t );
        ~VTKField();

        std::string              getName();
        std::string              getType();
        std::string              getLocation();
        std::string              getCodification();
        uint8_t                  getComponents();
        uint64_t                 getElements();
        uint64_t                 getSize();
        uint64_t                 getOffset();
        uint64_t                 getNbytes();
        std::fstream::pos_type   getPosition(); 

        void                     setName( std::string ) ;
        void                     setType( std::string ) ;
        void                     setLocation( std::string ) ;
        void                     setCodification( std::string ) ;
        void                     setComponents( uint8_t ) ;
        void                     setElements( uint64_t ) ;
        void                     setOffset( uint64_t ) ;
        void                     setPosition( std::fstream::pos_type ) ;

};

class VTK{

    protected:


        // members ---------------------------------------------------------------------- //
    protected:
        FileHandler                     fh ;                        /**< File_Handler for Input and Output */
        uint64_t                        nr_points ;                 /**< Number of vertices */
        uint64_t                        nr_cells  ;                 /**< Number of Cells */
        int                             nr_procs  ;                 /**< Number of parallel processes  */
        int                             my_proc   ;                 /**< My process id */

        std::string                     HeaderType ;                /**< UInt32 or UInt64_t */

        std::vector<VTKField>           geometry ;                  /**< Geometry fields */
        std::string                     GeomCodex ;                 /**< Geometry codex */

        std::vector<VTKField>           data ;                      /**< Data fields */
        std::string                     DataCodex ;                 /**< Data codex */

        // methods ----------------------------------------------------------------------- //
    public:
        VTK( );
        VTK( std::string dir_, std::string name_ );
        virtual ~VTK( );

        void                            setHeaderType( std::string );
        std::string                     getHeaderType(  );

        void                            setNames( std::string , std::string ) ;
        void                            setCounter( int c_=0 ) ;
        void                            setParallel( int , int ) ;

        void                            setCodex( std::string );
        void                            setGeomCodex( std::string );
        void                            setDataCodex( std::string );

        VTKField*                       addData( std::string , int , std::string , std::string ) ;
        VTKField*                       addData( std::string , int , std::string , std::string , std::string ) ;
        void                            removeData( std::string ) ;

        void                            read() ;

        virtual void                    readMetaData() = 0 ; 
        void                            readData() ;

        void                            write()  ;
        virtual void                    writeMetaData() = 0 ;
        void                            writeData() ;

        virtual void                    writeCollection() = 0 ;

    protected:
        //For Writing
        void                            writeDataHeader( std::fstream &, bool parallel=false ) ;
        void                            writeDataArray( std::fstream &, VTKField &) ;
        void                            writePDataArray( std::fstream &, VTKField &) ;

        virtual void                    flush( std::fstream &, std::string , std::string ) =0 ;

        //For Reading
        void                            readDataHeader( std::fstream &) ;
        bool                            readDataArray( std::fstream &, VTKField &);

        virtual  void                   absorb( std::fstream &, std::string , std::string ) =0 ;

        //General Purpose
        bool                            getFieldByName( const std::string &, VTKField*& ) ;
        void                            calcAppendedOffsets() ;

};

template <class Derived >
class VTKUnstructuredGrid : public VTK{

    protected:
        uint64_t                        nconnectivity ;             /**< size of the connectivity information */

    protected:
        VTKUnstructuredGrid();
        VTKUnstructuredGrid( std::string , std::string ) ;
        ~VTKUnstructuredGrid();

        void                            writeCollection() ;  

        void                            flush(  std::fstream &, std::string , std::string ) ; //CRTP
        void                            absorb( std::fstream &, std::string , std::string ) ; //CRTP
        uint64_t                        calcSizeConnectivity( ) ;

    public:
        void                            readMetaData() ;
        void                            writeMetaData() ;

        void                            setDimensions( uint64_t , uint64_t , uint64_t ) ;
        void                            setGeomTypes( std::string , std::string , std::string , std::string ) ;

        uint64_t                        getNConnectivity( ) ; 

};

template <class Derived>
class VTKRectilinearGrid : public VTK{

    typedef std::array<std::array<int,2>,2> extension2D_t ;         /**< typedef to describe min and max indices in 2D of restilinear grid */
    typedef std::array<std::array<int,2>,3> extension3D_t ;         /**< typedef to describe min and max indices in 3D of restilinear grid */

    protected:
    int                             dimensions ;                /**< dimensions of the grid [2/3] */
    extension3D_t                   local_index ;               /**< min and max indices of local grid */
    extension3D_t                   global_index ;              /**< min and max indices of global grid */
    std::vector<extension3D_t>      proc_index ;                /**< global indices of each processors */

    protected:
    VTKRectilinearGrid();
    VTKRectilinearGrid( std::string , std::string  );
    VTKRectilinearGrid( std::string , std::string , std::string, int, int, int, int, int, int );
    VTKRectilinearGrid( std::string , std::string , std::string, int, int, int );
    VTKRectilinearGrid( std::string , std::string , std::string, int, int, int, int );
    VTKRectilinearGrid( std::string , std::string , std::string, int, int );
    ~VTKRectilinearGrid();

    void                            writeCollection() ;  

    void                            flush(  std::fstream &, std::string , std::string ) ; //CRTP
    void                            absorb( std::fstream &, std::string , std::string ) ; //CRTP

    public:
    void                            readMetaData() ;
    void                            writeMetaData() ;

    void                            setDimensions( int, int, int, int, int, int ) ;
    void                            setDimensions( int, int, int ) ;
    void                            setDimensions( int, int, int, int ) ;
    void                            setDimensions( int, int ) ;

    void                            setGlobalDimensions( int, int, int ) ;
    void                            setGlobalDimensions( int, int ) ;

    void                            setGeomTypes( std::string ) ;

    void                            setGlobalIndex( std::vector<extension3D_t> ) ;
    void                            setGlobalIndex( std::vector<extension2D_t> ) ;

};

namespace VTKUtils{
    uint8_t                         sizeOfType( std::string type ) ;
    uint8_t                         getNNodeInElement( uint8_t t ) ;
    bool                            convertStringToDataArray( std::string &, VTKField &) ;
    void                            convertDataArrayToString( std::string &, VTKField &) ;
    void                            convertPDataArrayToString( std::string &, VTKField &) ;
    
    template<class T> 
    std::string                     whichType( T ) ;
    
    template<class T> 
    std::string                     whichType( std::vector<T> ) ;
    
    template<class T, size_t d>
    std::string                     whichType( std::array<T,d> ) ;
}


#include"VTKUtils.tpp"
#include"VTKUnstructured.tpp"
#include"VTKRectilinear.tpp"


#endif
