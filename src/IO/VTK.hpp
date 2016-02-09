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

/*!
 * @ingroup VTKEnums
 * Enum class defining types of fields whic may be written through class VTK
 */
enum class VTKFieldType {
    UNDEFINED = -1,
    SCALAR = 1,
    VECTOR = 3,
    CONSTANT = 4,
    VARIABLE = 5
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

class VTKTypes{

    private:
        static std::unordered_map<std::type_index, VTKDataType> m_types;  /**< map conatining registered data types */

    public:
        static uint8_t                  sizeOfType( const VTKDataType & type );

        template<typename T>
        static VTKDataType              registerType();

        template<typename T>
        static VTKDataType              registerType(VTKDataType VTKType);

        static VTKDataType              whichType( const std::type_info & ) ;

        template<class T>
        static VTKDataType              whichType( T ) ;

        template<class T>
        static VTKDataType              whichType( std::vector<T> ) ;

        template<class T, size_t d>
        static VTKDataType              whichType( std::array<T,d> ) ;


};

class VTKFieldMetaData{

    private:
        uint64_t                m_size ;                        /**< size of the field */
        const std::type_info&   m_type ;                        /**< tye of the field */

    public:
        VTKFieldMetaData( uint64_t, const std::type_info &);
        uint64_t                getSize() const;
        const std::type_info&   getType() const;

};

class VTKField{

    //members
    protected:
        std::string              name;                      /**< name of the field */
        VTKFieldType             fieldType;                      /**< type of field [ VTKFieldType::SCALAR/VECTOR/CONSTANT/VARIABLE ] */
        uint8_t                  components;                /**< type of field [ VTKFieldType::SCALAR/VECTOR/CONSTANT/VARIABLE ] */
        VTKDataType              dataType;                      /**< type of data [  VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ]] */
        VTKLocation              location;                  /**< cell or point data [ VTKLocation::CELL/VTKLocation::POINT] */
        VTKFormat                codification ;             /**< Type of codification [VTKFormat::ASCII, VTKFormat::APPENDED] */
        uint64_t                 nr_elements;               /**< nr of cells or points */
        uint64_t                 offset;                    /**< offset in the appended section */
        std::fstream::pos_type   position;                  /**< position in file */

        bool                     derived;                   /**< true if derived class (storing a pointer to data vector) is used, false if base class (using interface) */
        bool                     implicitKnown;             /**< true if class storing the Field is awre about the data to be written */

        //methods
    public:
        virtual ~VTKField();

        VTKField();
        VTKField( const VTKField &);
        VTKField( std::string );

        VTKField& operator=( const VTKField & );

        std::string              getName() const;
        VTKDataType              getDataType() const;
        VTKFieldType             getFieldType() const;
        VTKLocation              getLocation() const;
        VTKFormat                getCodification() const;
        uint8_t                  getComponents() const;
        uint64_t                 getElements() const;
        uint64_t                 getSize() const;
        uint64_t                 getOffset() const;
        uint64_t                 getNbytes() const;
        std::fstream::pos_type   getPosition() const; 
        bool                     hasAllMetaData() const ;

        bool                     usesInterface() const ;
        bool                     autoWrite() const ;

        void                     setName( std::string ) ;
        void                     setDataType( VTKDataType ) ;
        void                     setFieldType( VTKFieldType ) ;
        void                     setLocation( VTKLocation ) ;
        void                     setCodification( VTKFormat ) ;
        void                     setComponents( uint8_t ) ;
        void                     setElements( uint64_t ) ;
        void                     setOffset( uint64_t ) ;
        void                     setPosition( std::fstream::pos_type ) ;
        void                     setImplicit( bool ) ;

        void                     importMetaData( const VTKFieldMetaData & ) ;
        virtual void             flushData( std::fstream &) const ;
        virtual void             absorbData( std::fstream &) const ;
};

template<class T>
class VTKFieldWithVector : public VTKField{
    private:
    std::vector<T>*     m_ptr;            /**< pointer to data */

    public:
    ~VTKFieldWithVector( ) ;

    VTKFieldWithVector( );
    VTKFieldWithVector( const VTKField &, std::vector<T> & );
    VTKFieldWithVector( std::string, std::vector<T> & );

    void    setData( std::vector<T> &) ;
    void    flushData( std::fstream &) const;
    void    absorbData( std::fstream &) const;

};

class VTK{

    protected:


        // members ---------------------------------------------------------------------- //
    protected:
        FileHandler                     fh ;                        /**< File_Handler for Input and Output */
        uint64_t                        nr_points ;                 /**< Number of vertices */
        uint64_t                        nr_cells  ;                 /**< Number of Cells */
        uint16_t                        nr_procs  ;                 /**< Number of parallel processes  */
        uint16_t                        my_proc   ;                 /**< My process id */

        std::string                     HeaderType ;                /**< UInt32 or UInt64_t */

        std::vector<VTKField*>          geometry ;                  /**< Geometry fields */
        VTKFormat                       GeomCodex ;                 /**< Geometry codex */

        std::vector<VTKField*>          data ;                      /**< Data fields */
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
        void                            setParallel( uint16_t , uint16_t ) ;

        void                            setCodex( VTKFormat );
        void                            setGeomCodex( VTKFormat );
        void                            setDataCodex( VTKFormat );

        VTKField**                      addData( std::string ) ;
        VTKField**                      addData( std::string, VTKFieldType, VTKLocation ) ;
        VTKField**                      addData( std::string, VTKFieldType, VTKLocation, VTKDataType ) ;

        template<class T>
        VTKField**                      addData( std::string, std::vector<T> & ) ;
        template<class T>
        VTKField**                      addData( std::string, VTKFieldType, VTKLocation, std::vector<T> & ) ;

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
        virtual void                    writeFieldData( std::fstream &, VTKField &) ; 

        //For Reading
        void                            readDataHeader( std::fstream &) ;
        bool                            readDataArray( std::fstream &, VTKField &);
        virtual void                    readFieldData( std::fstream &, VTKField &) ; 

        //General Purpose
        bool                            getFieldByName( const std::string &, VTKField**& ) ;
        void                            calcAppendedOffsets() ;
        void                            getMissingMetaData() ;
        virtual void                    setMissingGlobalData() ;

        //Interface methods
        virtual void                    flushData( std::fstream &, VTKFormat , std::string )  ;
        virtual  void                   absorbData( std::fstream &, VTKFormat , std::string )  ;
        virtual const VTKFieldMetaData  getMetaData( std::string )  ;

};

class VTKUnstructuredGrid : public VTK{

    protected:
    uint64_t                        nconnectivity ;             /**< size of the connectivity information */
    VTKElementType                  homogeneousType ;           /**< type of element mesh is made of */

    public:
    ~VTKUnstructuredGrid();

    VTKUnstructuredGrid();
    VTKUnstructuredGrid( std::string , std::string ) ;
    VTKUnstructuredGrid( std::string , std::string, VTKElementType ) ;

    template<class T0, class T1>
    VTKUnstructuredGrid( std::string , std::string, VTKElementType, std::vector<T0> &, std::vector<T1> & ) ;

    protected:
    void                            writeFieldData( std::fstream &, VTKField &) ; 
    void                            readFieldData( std::fstream &, VTKField &) ; 

    void                            writeCollection() ;  
    uint64_t                        calcSizeConnectivity( ) ;
    void                            setMissingGlobalData() ;

    public:
    void                            readMetaData() ;
    void                            writeMetaData() ;

    void                            setElementType( VTKElementType ) ;
    void                            setDimensions( uint64_t , uint64_t , uint64_t nconn_=0 ) ;
    void                            setDimensions( uint64_t , uint64_t , VTKElementType ) ;

    void                            setGeomTypes( VTKDataType , VTKDataType , VTKDataType , VTKDataType ) ;

    template<class T0>
    void                            setGeomData( std::string, std::vector<T0> & ) ;
    template<class T0, class T1>
    void                            setGeomData( std::vector<T0> &, std::vector<T1> & ) ;

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

    void                            setMissingGlobalData() ;

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

/*!
 * @ingroup  VisualizationToolKit
 * @brief Utility fuctions for VTK
 */
namespace vtk{
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

    template<class T>
    void                            allocate( std::vector<T> &, int) ;

    template<class T>
    void                            allocate( T &, int) ;
}

}

#include"VTK.tpp"
#include"VTKTypes.tpp"
#include"VTKField.tpp"
#include"VTKUtils.tpp"


#endif
