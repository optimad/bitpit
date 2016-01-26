
#ifndef __VTKWRAP_HH__
#define __VTKWRAP_HH__

#include <boost/mpl/vector.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/variant.hpp>

#include "VTK.hpp"
#include "VTKContainers.hpp"

/*!
 * @ingroup    VisualizationToolKit
 * @{
 */


class VTKUnstructuredVec : public VTKUnstructuredGrid{


    friend VTKUnstructuredGrid ;                                        /**< provides friendship to base class */

    private:
    struct ufield{
        std::string                 name;                               /**< name of the field */
        VTKContainers::Variants     DPtr ;                              /**< pointer to the field */

        ufield() ;

        template<class T>
        ufield( std::string name_, std::vector<T>& data_) ;

    };

    struct stream_visitor : public boost::static_visitor<>{

        private:
            std::fstream*       str ;                               /**< stream for reading/writing */
            int                 components ;                        /**< number of components of field */
            uint64_t            size ;                              /**< number of elements */
            VTKFormat           codex ;                             /**< codex of VTK file [VTKFormat::ASCII/VTKFormat::APPENDED */
            std::string         name ;                              /**< name of the field */
            std::string         task ;                              /**< task to be performed ["read"/"write"] */

        public:

            template <typename T>
                void                operator()( T* ) const;

            template <typename T>
                void                operator()( std::vector< std::vector<T> > * ) const;


            void                setStream( std::fstream& );
            void                setCodex( VTKFormat );
            void                setTask( std::string );
            void                setName( std::string );
            void                setSize( uint64_t );
            void                setComponents( int );


    };

    VTKElementType          type ;                                  /**< Type of element used in unstructured grid */
    std::vector<ufield>     adata;                                  /**< All field data */

    void flushData( std::fstream &str, VTKFormat codex_, std::string name ) ;
    void absorbData( std::fstream &str, VTKFormat codex_, std::string name ) ;

    public:
    VTKUnstructuredVec( ) ;
    VTKUnstructuredVec( std::string, std::string, VTKFormat, VTKElementType ) ;

    template< class T0, class T1>
    VTKUnstructuredVec( std::string, std::string, VTKFormat, VTKElementType, std::vector<T0> &, std::vector<T1> &) ;

    ~VTKUnstructuredVec( ) ;

    template< class T>
    void    addData( std::vector<T> &, std::string, VTKLocation ) ; 

    template< class T>
    void    addData( std::vector< std::array<T,3> > &, std::string, VTKLocation ) ; 

    template< class T>
    void    addData( std::vector< std::vector<T> > &, std::string, VTKLocation ) ; 

    template< class T>
    void    linkData( std::vector<T> &, std::string ) ; 

    void    write() ;

    private:
    bool    getFieldByName( const std::string &, ufield*& ) ;


};

/*!
 * @}
 */

#include"VTKWrappers.tpp"

#endif
