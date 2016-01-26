#include "VTK.hpp"

/*!
 * @ingroup     VisualizationToolKit
 * @{
 *
 * @class       VTKFieldMetaData
 * @brief       A class used in the interface VTK::getFieldMetaData( std::string )
 *
 */

/*! 
 * Constructor with all members set
 * @param[in] size the entire size of field date, e.g. in case of a vector field on nodes size = NNodes x 3 
 * @param[in] type the type of the basic data used
 */
VTKFieldMetaData::VTKFieldMetaData( uint64_t size, const std::type_info &type): m_size(size), m_type(type){
};

/*! 
 * Get the size of field
 * @return size of field
 */
uint64_t VTKFieldMetaData::getSize() const{
    return m_size;
}

/*! 
 * Get the type of data used in field
 * @return data type of field
 */
const std::type_info& VTKFieldMetaData::getType() const{
    return m_type;
}

/*! 
 * @} 
 */
