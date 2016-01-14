
#ifndef __CLASS_VTK_WRAP_HH__
#define __CLASS_VTK_WRAP_HH__

#include <boost/mpl/vector.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/variant.hpp>

#include"Class_VTK.hpp"


/*!
 * @ingroup    VisualizationToolKit
 * @{
 */

/*!
 * @typedef     VtkVecOfVector
 * All supported variants of vector of POD and vector of vector of POD
 */
typedef boost::mpl::vector<
    std::vector<int8_t>*  ,
    std::vector<int16_t>* ,
    std::vector<int32_t>* ,
    std::vector<int64_t>* ,
    std::vector<uint8_t>* ,
    std::vector<uint16_t>*, 
    std::vector<uint32_t>*, 
    std::vector<uint64_t>*,

    std::vector<float>*   ,
    std::vector<double>*  ,

    std::vector< std::vector<int8_t> >*  , 
    std::vector< std::vector<int16_t> >*  , 
    std::vector< std::vector<int32_t> >*  , 
    std::vector< std::vector<int64_t> >*  , 
    std::vector< std::vector<uint8_t> >*  , 
    std::vector< std::vector<uint16_t> >* ,
    std::vector< std::vector<uint32_t> >* ,
    std::vector< std::vector<uint64_t> >* ,

    std::vector< std::vector<float> >*  , 
    std::vector< std::vector<double> >*
>::type VtkVecOfVector;

/*!
 * @typedef     VtkVecOfArray
 * All supported variants of vector of array of POD
 */
typedef boost::mpl::vector<
    std::vector< std::array<int16_t,3> >*  , 
    std::vector< std::array<int32_t,3> >*  , 
    std::vector< std::array<int64_t,3> >*  , 
    std::vector< std::array<uint16_t,3> >* ,
    std::vector< std::array<uint32_t,3> >* ,
    std::vector< std::array<uint64_t,3> >* ,

    std::vector< std::array<int16_t,4> >*  ,
    std::vector< std::array<int32_t,4> >*  ,
    std::vector< std::array<int64_t,4> >*  ,
    std::vector< std::array<uint16_t,4> >* ,
    std::vector< std::array<uint32_t,4> >* ,
    std::vector< std::array<uint64_t,4> >* ,

    std::vector< std::array<int16_t,8> >*  ,
    std::vector< std::array<int32_t,8> >*  ,
    std::vector< std::array<int64_t,8> >*  ,
    std::vector< std::array<uint16_t,8> >* ,
    std::vector< std::array<uint32_t,8> >* ,
    std::vector< std::array<uint64_t,8> >* ,

    std::vector< std::array<float,3> >*    ,
    std::vector< std::array<double,3> >*   
>::type VtkVecOfArray;

/*!
 * @typedef     VtkVecContainers
 * All supported containers by VtkUnstrVec
 */
typedef boost::mpl::copy< VtkVecOfVector::type, boost::mpl::back_inserter<VtkVecOfArray> > ::type VtkVecContainers ;

/*!
 * @typedef     VtkVecVariants
 * Boost variant over VtkVecContainers
 */
typedef boost::make_variant_over< VtkVecContainers >::type  VtkVecVariants ;

/*!
 * @}
 */

class VtkUnstrVec : public VTK_UnstructuredGrid<VtkUnstrVec>{


    friend VTK_UnstructuredGrid<VtkUnstrVec> ;                      /**< provides friendship to base class */

    private:
    struct ufield{
        std::string             name;                               /**< name of the field */
        VtkVecVariants          DPtr ;                              /**< pointer to the field */

        ufield() ;

        template<class T>
        ufield( std::string name_, std::vector<T>& data_) ;

    };

    struct stream_visitor : public boost::static_visitor<>{

        private:
            std::fstream*       str ;                               /**< stream for reading/writing */
            uint8_t             components ;                        /**< number of components of field */
            uint64_t            size ;                              /**< number of elements */
            std::string         codex ;                             /**< codex of VTK file ["appended"/"ascii"] */
            std::string         name ;                              /**< name of the field */
            std::string         task ;                              /**< task to be performed ["read"/"write"] */

        public:

            template <typename T>
                void                operator()( T* ) const;

            template <typename T>
                void                operator()( std::vector< std::vector<T> > * ) const;


            void                SetStream( std::fstream& );
            void                SetCodex( std::string );
            void                SetTask( std::string );
            void                SetName( std::string );
            void                SetSize( uint64_t );
            void                SetComponents( uint8_t );


    };

    uint8_t                 type ;                                  /**< Type of element used in unstructured grid */
    std::vector<ufield>     adata;                                  /**< All field data */

    void Flush( std::fstream &str, std::string codex_, std::string name ) ;
    void Absorb( std::fstream &str, std::string codex_, std::string name ) ;

    public:
    VtkUnstrVec( ) ;
    VtkUnstrVec( std::string , std::string , std::string , uint8_t ) ;

    template< class T0, class T1>
    VtkUnstrVec( std::string , std::string , std::string , uint8_t , std::vector<T0> &, std::vector<T1> &) ;

    ~VtkUnstrVec( ) ;

    template< class T>
    void    AddData( std::vector<T> &, std::string , std::string ) ; 

    template< class T>
    void    AddData( std::vector< std::array<T,3> > &, std::string , std::string ) ; 

    template< class T>
    void    AddData( std::vector< std::vector<T> > &, std::string , std::string ) ; 

    template< class T>
    void    LinkData( std::vector<T> &, std::string ) ; 

    void    Write() ;

    private:
    bool    GetFieldByName( const std::string &, ufield*& ) ;


};


#include"Class_VTK_Wrappers.tpp"

#endif
