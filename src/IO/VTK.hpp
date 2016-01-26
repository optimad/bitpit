

#ifndef __VTK__HH__
#define __VTK__HH__

#include <typeinfo>
#include <vector>
#include <array>

#include "BitPit_common.hpp"
#include "Generic_IO.hpp"
#include "FileHandler.hpp"

enum class VTKFieldType {
    UNDEFINED = -1,
    SCALAR = 1,
    VECTOR = 3
enum class VTKDataType {
    UNDEFINED,
    Int8     ,
    Int16    ,
    Int32    ,
    Int64    ,
    UInt8    ,
    UInt16   ,
    UInt32   ,
    UInt64   ,
    Float32  ,
    Float64      
};

enum class VTKFormat {
    UNDEFINED,
    ASCII,
    APPENDED
};

enum class VTKLocation {
    UNDEFINED,
    CELL,
    POINT
};

enum class VTKElementType {
    UNDEFINED  = -1,
    VERTEX     = 1,
    LINE       = 3,
    TRIANGLE   = 5,
    POLYGON    = 7,
    PIXEL      = 8,
    QUAD       = 9,
    TETRA      = 10,
    VOXEL      = 11,
    HEXAHEDRON = 12,
    WEDGE      = 13,
    PYRAMID    = 14,
    POLYHEDRON = 42
};

class VTKFieldMetaData{

    private:
    uint64_t                m_size ;
    const std::type_info&   m_type ;

    public:
    VTKFieldMetaData( uint64_t, const std::type_info &);
    uint64_t                getSize() const;
    const std::type_info&   getType() const;

};

class VTKField{

    //members
    protected:
        std::string              name;                      /**< name of the field */
        VTKFieldType             components;                /**< nr of components of field, options[ 1, 3 ] */
        VTKDataType              type;                      /**< type of data, options [  VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ]] */
        VTKLocation              location;                  /**< cell or point data, [ VTKLocation::CELL/VTKLocation::Point] */
        VTKFormat                codification ;             /**< Type of codification [VTKFormat::ASCII, VTKFormat::APPENDED] */
        uint64_t                 nr_elements;               /**< nr of cells or points */
        uint64_t                 offset;                    /**< offset in the appended section */
        std::fstream::pos_type   position;                  /**< position in file */

        //methods
    public:
        VTKField();
        VTKField( std::string, VTKFieldType, VTKLocation );
        VTKField( std::string, VTKFieldType, VTKLocation, VTKDataType );
        VTKField( std::string, VTKFieldType, VTKLocation, VTKDataType , VTKFormat , uint64_t );
        ~VTKField();

        std::string              getName() const;
        VTKDataType              getType() const;
        VTKLocation              getLocation() const;
        VTKFormat                getCodification() const;
        VTKFieldType             getComponents() const;
        uint64_t                 getElements() const;
        uint64_t                 getSize() const;
        uint64_t                 getOffset() const;
        uint64_t                 getNbytes() const;
        std::fstream::pos_type   getPosition() const; 
        bool                     hasAllMetaData() const ;

        void                     setName( std::string ) ;
        void                     setType( VTKDataType ) ;
        void                     setLocation( VTKLocation ) ;
        void                     setCodification( VTKFormat ) ;
        void                     setComponents( VTKFieldType ) ;
        void                     setElements( uint64_t ) ;
        void                     setOffset( uint64_t ) ;
        void                     setPosition( std::fstream::pos_type ) ;

        void                     importMetaData( const VTKFieldMetaData & ) ;
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
        VTKFormat                       GeomCodex ;                 /**< Geometry codex */

        std::vector<VTKField>           data ;                      /**< Data fields */
        VTKFormat                       DataCodex ;                 /**< Data codex */

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

        void                            setCodex( VTKFormat );
        void                            setGeomCodex( VTKFormat );
        void                            setDataCodex( VTKFormat );

        VTKField*                       addData( std::string, VTKFieldType, VTKLocation ) ;
        VTKField*                       addData( std::string, VTKFieldType, VTKLocation, VTKDataType ) ;
        VTKField*                       addData( std::string, VTKFieldType, VTKLocation, VTKDataType, VTKFormat ) ;

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

        //For Reading
        void                            readDataHeader( std::fstream &) ;
        bool                            readDataArray( std::fstream &, VTKField &);

        //General Purpose
        bool                            getFieldByName( const std::string &, VTKField*& ) ;
        void                            calcAppendedOffsets() ;
        void                            getMissingMetaData() ;

        //Interface methods
        virtual void                    flushData( std::fstream &, VTKFormat , std::string )  ;
        virtual  void                   absorbData( std::fstream &, VTKFormat , std::string )  ;
        virtual const VTKFieldMetaData  getFieldMetaData( std::string )  ;

};

class VTKUnstructuredGrid : public VTK{

    protected:
        uint64_t                        nconnectivity ;             /**< size of the connectivity information */

    protected:
        VTKUnstructuredGrid();
        VTKUnstructuredGrid( std::string , std::string ) ;
        ~VTKUnstructuredGrid();

        void                            writeCollection() ;  
        uint64_t                        calcSizeConnectivity( ) ;

    public:
        void                            readMetaData() ;
        void                            writeMetaData() ;

        void                            setDimensions( uint64_t , uint64_t , uint64_t ) ;
        void                            setGeomTypes( VTKDataType , VTKDataType , VTKDataType , VTKDataType ) ;

        uint64_t                        getNConnectivity( ) ; 

};

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
    VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int, int, int, int );
    VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int );
    VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int, int );
    VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int );
    ~VTKRectilinearGrid();

    void                            writeCollection() ;  


    public:
    void                            readMetaData() ;
    void                            writeMetaData() ;

    void                            setDimensions( int, int, int, int, int, int ) ;
    void                            setDimensions( int, int, int ) ;
    void                            setDimensions( int, int, int, int ) ;
    void                            setDimensions( int, int ) ;

    void                            setGlobalDimensions( int, int, int ) ;
    void                            setGlobalDimensions( int, int ) ;

    void                            setGeomTypes( VTKDataType ) ;

    void                            setGlobalIndex( std::vector<extension3D_t> ) ;
    void                            setGlobalIndex( std::vector<extension2D_t> ) ;

};

namespace VTKUtils{
    uint8_t                         sizeOfType( const VTKDataType & ) ;
    uint8_t                         getNNodeInElement( const VTKElementType & ) ;

    std::string                     convertDataArrayToString( const VTKField & ) ;
    std::string                     convertPDataArrayToString( const VTKField & ) ;

    bool                            convertStringToDataArray( const std::string &, VTKField &) ;

    std::string                     convertEnumToString( const VTKLocation & ) ;
    std::string                     convertEnumToString( const VTKFormat & ) ;
    std::string                     convertEnumToString( const VTKDataType & ) ;

    bool                            convertStringToEnum( const std::string &, VTKLocation & ) ;
    bool                            convertStringToEnum( const std::string &, VTKFormat & ) ;
    bool                            convertStringToEnum( const std::string &, VTKDataType &) ;

    VTKDataType                     whichType( const std::type_info & ) ;

    template<class T> 
    VTKDataType                     whichType( T ) ;

    template<class T> 
    VTKDataType                     whichType( std::vector<T> ) ;

    template<class T, size_t d>
    VTKDataType                     whichType( std::array<T,d> ) ;
}


#include"VTKUtils.tpp"


#endif
