/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/


#ifndef __BITPIT_VTK_HPP__
#define __BITPIT_VTK_HPP__

#include <typeinfo>
#include <type_traits>
#include <vector>
#include <array>
#include <typeindex>
#include <unordered_map>

#include "bitpit_common.hpp"
#include "GenericIO.hpp"
#include "FileHandler.hpp"

namespace bitpit{

class VTK;
class VTKField;

/*!
 * @ingroup VTKEnums
 * Enum class defining different modes for writing VTK files
 */
enum class VTKWriteMode {
    DEFAULT = 0,
    NO_INCREMENT =1,
    NO_SERIES =2
};

/*!
 * @ingroup VTKEnums
 * Enum class defining types of fields whic may be written through class VTK
 */
enum class VTKFieldType {
    UNDEFINED = -1,
    SCALAR = 1,
    VECTOR = 3,
    KNOWN_BY_CLASS = 4,
};

/*!
 * @ingroup VTKEnums
 * Enum class defining basic data types of fields which may be written through VTK
 */
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

/*!
 * @ingroup VTKEnums
 * Enum class defining the VTK format to be used for writing fields 
 */
enum class VTKFormat {
    UNDEFINED,
    ASCII,
    APPENDED
};

/*!
 * @ingroup VTKEnums
 * Enum class defining wheather data is stored at cells or nodes
 */
enum class VTKLocation {
    UNDEFINED,
    CELL,
    POINT
};

/*!
 * @ingroup VTKEnums
 * Enum class listing different element types supported by VTKUnstructuredGrid
 */
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

/*!
 * @ingroup VTKEnums
 * Enum class listing different geometry fields used by VTKUnstructuredGrid
 */
enum class VTKUnstructuredField {
    POINTS      = 0,
    OFFSETS     = 1,
    TYPES       = 2,
    CONNECTIVITY =3 
};

/*!
 * @ingroup VTKEnums
 * Enum class listing different geometry fields used by VTKRectilinearGrid
 */
enum class VTKRectilinearField {
    X_COORDS    = 0,
    Y_COORDS    = 1,
    Z_COORDS    = 2 
};

class VTKTypes{

    private:
        static std::unordered_map<std::type_index, VTKDataType> m_types;  /**< map conatining registered data types */

    public:
        static uint8_t          sizeOfType( const VTKDataType & type );

        template<typename T>
        static VTKDataType      registerType();

        template<typename T>
        static VTKDataType      registerType(VTKDataType VTKType);

        static VTKDataType      whichType( const std::type_info & ) ;

        template<class T>
        static VTKDataType      whichType( const T &) ;

        template<class T>
        static VTKDataType      whichType( const std::vector<T> & ) ;

        template<class T, size_t d>
        static VTKDataType      whichType( const std::array<T,d> & ) ;
};

class VTKBaseContainer{
    private:

    public:
        VTKBaseContainer( ) ;
        VTKBaseContainer( const VTKBaseContainer &) = default;
        virtual ~VTKBaseContainer( ) ;

        virtual VTKBaseContainer *  clone() const = 0 ;

        virtual void                flushData( std::fstream &, VTKFormat) =0 ;
        virtual void                absorbData( std::fstream &, VTKFormat, uint64_t, uint8_t) =0 ;
};

template<class T>
class VTKVectorContainer : public VTKBaseContainer{
    private:
        std::vector<T>*         m_ptr ;                     /**< pointer to data */

    public:
        VTKVectorContainer( std::vector<T> &) ;
        VTKVectorContainer( const VTKVectorContainer &);
        ~VTKVectorContainer( ) ;

        VTKVectorContainer*     clone() const ;

        void                    flushData( std::fstream &, VTKFormat) ;
        void                    absorbData( std::fstream &, VTKFormat, uint64_t, uint8_t) ;
        void                    resize( std::true_type, uint64_t , uint8_t) ;
        void                    resize( std::false_type, uint64_t , uint8_t) ;
};

class VTKBaseStreamer{ 

    private:

    public:
        virtual void            flushData( std::fstream &, std::string, VTKFormat)  ;
        virtual void            absorbData( std::fstream &, std::string, VTKFormat, uint64_t, uint8_t )  ;
};

class VTKNativeStreamer : public VTKBaseStreamer {

    private:
        std::unordered_map<std::string,std::unique_ptr<VTKBaseContainer> >         m_field ; /**< association between name of field and conatiner */

    public:
        VTKNativeStreamer();
        VTKNativeStreamer( const VTKNativeStreamer & );
        ~VTKNativeStreamer();

        template< class T>
        void                    addData( std::string, std::vector<T> & ) ;
        void                    removeData( std::string ) ;
        void                    flushData( std::fstream &, std::string, VTKFormat) ;
        void                    absorbData( std::fstream &, std::string, VTKFormat, uint64_t, uint8_t) ;
};

class VTKField{

    //members
    protected:
        std::string             m_name;                     /**< name of the field */
        VTKFieldType            m_fieldType;                /**< type of field [ VTKFieldType::SCALAR/VECTOR/KNOWN_BY_CLASS ] */
        VTKDataType             m_dataType;                 /**< type of data [  VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ]] */
        VTKLocation             m_location;                 /**< cell or point data [ VTKLocation::CELL/VTKLocation::POINT] */
        VTKFormat               m_codification ;            /**< Type of codification [VTKFormat::ASCII, VTKFormat::APPENDED] */
        uint64_t                m_offset;                   /**< offset in the appended section */
        std::fstream::pos_type  m_position;                 /**< position in file */
        VTKBaseStreamer*        m_streamer;                 /**< pointer to streamer */
        bool                    m_enabled;                  /**< if field is enabled for read/write */
                                

        //methods
    public:
        virtual ~VTKField();

        VTKField();
        VTKField( const VTKField &);
        VTKField( std::string );

        VTKField& operator=( const VTKField & );
                                
        std::string             getName() const;
        VTKDataType             getDataType() const;
        VTKFieldType            getFieldType() const;
        VTKLocation             getLocation() const;
        VTKFormat               getCodification() const;
        uint64_t                getOffset() const;
        std::fstream::pos_type  getPosition() const; 
        bool                    isEnabled() const ;
        bool                    hasAllMetaData() const ;

        void                    setName( std::string ) ;
        void                    setDataType( VTKDataType ) ;
        void                    setFieldType( VTKFieldType ) ;
        void                    setLocation( VTKLocation ) ;
        void                    setCodification( VTKFormat ) ;
        void                    setOffset( uint64_t ) ;
        void                    setPosition( std::fstream::pos_type ) ;
        void                    setStreamer( VTKBaseStreamer& ) ;
        void                    enable() ;
        void                    disable() ;

        void                    read( std::fstream&, uint64_t, uint8_t ) const ;
        void                    write( std::fstream& ) const ;
};

class VTK{

    protected:
        FileHandler             m_fh ;                      /**< File_Handler for Input and Output */
        uint64_t                m_points ;                  /**< Number of vertices */
        uint64_t                m_cells  ;                  /**< Number of Cells */
        uint16_t                m_procs  ;                  /**< Number of parallel processes  */
        uint16_t                m_rank   ;                  /**< My process id */

        std::string             m_headerType ;              /**< UInt32 or UInt64_t */

        std::vector<VTKField>   m_geometry ;                /**< Geometry fields */
        VTKFormat               m_geomCodex ;               /**< Geometry codex */

        std::vector<VTKField>   m_data ;                    /**< Data fields */
        VTKFormat               m_dataCodex ;               /**< Data codex */

        VTKNativeStreamer       m_nativeStreamer;           /**< native streamer for streaming data stored in std::vector<> */

    public:
        VTK( );
        VTK( std::string, std::string );
        virtual ~VTK( );

        void                    setHeaderType( std::string );
        std::string             getHeaderType(  ) const;

        std::string             getName(  ) const;
        std::string             getDirectory(  ) const;
        int                     getCounter(  ) const;

        void                    setNames( std::string , std::string ) ;
        void                    setName( std::string ) ;
        void                    setDirectory( std::string ) ;
        void                    setCounter( int c_=0 ) ;
        int                     unsetCounter( ) ;
        void                    setParallel( uint16_t , uint16_t ) ;

        void                    setCodex( VTKFormat );
        void                    setGeomCodex( VTKFormat );
        void                    setDataCodex( VTKFormat );

        VTKField&               addData( std::string, VTKBaseStreamer* = NULL ) ;
        VTKField&               addData( std::string, VTKFieldType, VTKLocation, VTKDataType, VTKBaseStreamer* =NULL ) ;

        template<class T>
        VTKField&               addData( std::string, std::vector<T> & ) ;
        template<class T>
        VTKField&               addData( std::string, VTKFieldType, VTKLocation, std::vector<T> & ) ;

        void                    removeData( std::string ) ;
        void                    enableData( std::string ) ;
        void                    disableData( std::string ) ;


        void                    read() ;
        virtual void            readMetaInformation() = 0 ; 
        void                    readData() ;

        void                    write( VTKWriteMode writeMode=VTKWriteMode::DEFAULT )  ;
        void                    write( std::string, VTKWriteMode writeMode=VTKWriteMode::NO_INCREMENT )  ;

        virtual void            writeMetaInformation() = 0 ;
        void                    writeData() ;
        virtual void            writeCollection() = 0 ;

    protected:
        //For Writing
        void                    writeDataHeader( std::fstream &, bool parallel=false ) ;
        void                    writeDataArray( std::fstream &, VTKField &) ;
        void                    writePDataArray( std::fstream &, VTKField &) ;

        //For Reading
        void                    readDataHeader( std::fstream &) ;
        bool                    readDataArray( std::fstream &, VTKField &);

        //General Purpose
        bool                    getDataByName( const std::string &, VTKField*& ) ;
        bool                    getGeomByName( const std::string &, VTKField*& ) ;
        void                    calcAppendedOffsets() ;
        virtual uint64_t        calcFieldSize( const VTKField &) =0;
        virtual uint64_t        calcFieldEntries( const VTKField &) =0;
        virtual uint8_t         calcFieldComponents( const VTKField &) =0;
        void                    checkAllFields() ;

};

class VTKUnstructuredGridStreamer :public VTKBaseStreamer{

    private:
        VTKElementType          m_homogeneousType ;         /**< the type of cells */
        long                    m_cells ;                   /**< numer of cells */
    
        void                    flushData( std::fstream &, std::string, VTKFormat) ;
    
    public:
        void                    setGrid( VTKElementType, long) ;
};

class VTKUnstructuredGrid : public VTK {

    protected:
        uint64_t                m_connectivity ;            /**< size of the connectivity information */
        VTKElementType          m_homogeneousType ;         /**< type of element mesh is made of */
        VTKUnstructuredGridStreamer     m_unstructuredStreamer;     /**< streamer if unstructured grid is of homogenous type */

    public:
    ~VTKUnstructuredGrid();

    VTKUnstructuredGrid();
    VTKUnstructuredGrid( std::string , std::string ) ;
    VTKUnstructuredGrid( std::string , std::string, VTKElementType ) ;

    protected:
        void                    writeCollection() ;  
        uint64_t                calcSizeConnectivity( ) ;

    public:
        void                    readMetaInformation() ;
        void                    writeMetaInformation() ;

        void                    setElementType( VTKElementType ) ;
        void                    setDimensions( uint64_t , uint64_t , uint64_t nconn_=0 ) ;
        void                    setDimensions( uint64_t , uint64_t , VTKElementType ) ;

        template<class T>
        void                    setGeomData( VTKUnstructuredField, std::vector<T> & ) ;
        void                    setGeomData( VTKUnstructuredField, VTKBaseStreamer* = NULL ) ;
        void                    setGeomData( VTKUnstructuredField, VTKDataType, VTKBaseStreamer* =NULL ) ;

        uint64_t                getNConnectivity( ) ; 
        uint64_t                calcFieldSize( const VTKField &) ;
        uint64_t                calcFieldEntries( const VTKField &) ;
        uint8_t                 calcFieldComponents( const VTKField &) ;
};

class VTKRectilinearGrid : public VTK{

    typedef std::array<std::array<int,2>,2> extension2D_t ; /**< typedef to describe min and max indices in 2D of restilinear grid */
    typedef std::array<std::array<int,2>,3> extension3D_t ; /**< typedef to describe min and max indices in 3D of restilinear grid */

    protected:
        int                     m_dimensions ;              /**< dimensions of the grid [2/3] */
        extension3D_t           m_localIndex ;              /**< min and max indices of local grid */
        extension3D_t           m_globalIndex ;             /**< min and max indices of global grid */
        std::vector<extension3D_t>      m_procIndex ;       /**< global indices of each processors */

    protected:
        VTKRectilinearGrid();
        VTKRectilinearGrid( std::string , std::string  );
        VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int, int, int, int );
        VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int );
        VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int, int, int );
        VTKRectilinearGrid( std::string , std::string , VTKFormat, int, int );
        ~VTKRectilinearGrid();

        void                    writeCollection() ;  

    public:
        void                    readMetaInformation() ;
        void                    writeMetaInformation() ;

        void                    setDimensions( int, int, int, int, int, int ) ;
        void                    setDimensions( int, int, int ) ;
        void                    setDimensions( int, int, int, int ) ;
        void                    setDimensions( int, int ) ;

        void                    setGlobalDimensions( int, int, int ) ;
        void                    setGlobalDimensions( int, int ) ;

        template<class T>
        void                    setGeomData( VTKRectilinearField, std::vector<T> & ) ;
        void                    setGeomData( VTKRectilinearField, VTKBaseStreamer* = NULL ) ;
        void                    setGeomData( VTKRectilinearField, VTKDataType, VTKBaseStreamer* =NULL ) ;

        void                    setGlobalIndex( std::vector<extension3D_t> ) ;
        void                    setGlobalIndex( std::vector<extension2D_t> ) ;

        uint64_t                calcFieldSize( const VTKField &) ;
        uint64_t                calcFieldEntries( const VTKField &) ;
        uint8_t                 calcFieldComponents( const VTKField &) ;
};

/*!
 * @ingroup  VisualizationToolKit
 * @brief Utility fuctions for VTK
 */
namespace vtk{
    uint8_t                     getElementNodeCount( const VTKElementType & ) ;

    std::string                 convertDataArrayToString( const VTKField & ) ;
    std::string                 convertPDataArrayToString( const VTKField & ) ;

    bool                        convertStringToDataArray( const std::string &, VTKField &) ;

    std::string                 convertEnumToString( const VTKLocation & ) ;
    std::string                 convertEnumToString( const VTKFormat & ) ;
    std::string                 convertEnumToString( const VTKDataType & ) ;

    bool                        convertStringToEnum( const std::string &, VTKLocation & ) ;
    bool                        convertStringToEnum( const std::string &, VTKFormat & ) ;
    bool                        convertStringToEnum( const std::string &, VTKDataType &) ;

    template<class T>
    void                        allocate( std::vector<T> &, int) ;

    template<class T>
    void                        allocate( T &, int) ;
}

}

#include"VTK.tpp"
#include"VTKTypes.tpp"
#include"VTKStreamer.tpp"
#include"VTKUtils.tpp"


#endif
