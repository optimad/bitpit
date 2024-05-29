/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

# ifndef __BITPIT_LEVELSET_COMMON_HPP__
# define __BITPIT_LEVELSET_COMMON_HPP__

// Standard Template Library
# include <array>
# include <limits>
# include <unordered_map>
# include <vector>

namespace bitpit{

/*!
 * Setting the size of the narrow band to LEVELSET_NARROW_BAND_UNLIMITED means that the whole
 * domain belongs to the narrow band.
 */
const double LEVELSET_NARROW_BAND_UNLIMITED = std::numeric_limits<double>::max();

/*!
 * @ingroup levelset
 * @brief namespace containing default values
 */
namespace levelSetDefaults{
    const double                            VALUE = 1.e18 ;             /**< Default value for levelset function */
    const std::array<double,3>              GRADIENT = {{0.,0.,0.}};    /**< Default value for levelset gradient */
    const std::array<double,3>              POINT = {{0.,0.,0.}};       /**< Default value for levelset projection point */
    const short                             SIGN = 1;                   /**< Default value for the sign */
    const double                            SIZE = 1.e18 ;              /**< Default value for surface feature */
    const int                               OBJECT = -1 ;               /**< Default value for closest object  */
    const int                               PART  = -1 ;                /**< Default value for closest patch  */
    const long                              SUPPORT = -1 ;              /**< Default value for closest support element */
    const std::array<double,3>              NORMAL = {{0.,0.,0.}};      /**< Default value for closest surface normal */
    const double                            NARROW_BAND_SIZE = -1 ;     /**< Default value for the narrow band size */
};

struct LevelSetInfo{
    double                                  value ;                     /**< Levelset value */
    std::array<double,3>                    gradient ;                  /**< Levelset gradient */

    LevelSetInfo() ;
    LevelSetInfo( double , const std::array<double,3> &) ;
};

/*!
 * @ingroup levelsetEnums
 * Enum class defining different type of zones.
 */
enum class LevelSetZone {
    NARROW_BAND,        //!< Narrow band zone
    BULK,               //!< Bulk zone
};

/*!
 * @ingroup levelsetEnums
 * Enum class defining different type of cell locations.
 */
enum class LevelSetCellLocation {
    UNKNOWN = -1,                //!< Unknown location
    NARROW_BAND_DISTANCE,        //!< Narrow band zone, the distance of the cell from the surface
                                 //!< is less than the narrow band size
    NARROW_BAND_INTERSECTED,     //!< Narrow band zone, the cell intersects the surface
    NARROW_BAND_NEIGHBOUR,       //!< Narrow band zone, on of the cell face neighbours intersect
                                 //!< the surface
    NARROW_BAND_UNDEFINED,       //!< Narrow band zone, the reason why the cell is inside the
                                 //!< narrow band is not defined
    BULK,                        //!< Bulk zone
};

/*!
 * @ingroup levelsetEnums
 * Enum class defining different boolean operations
 */
enum class LevelSetBooleanOperation{
    UNION =0,                                                           /**< Union between two objects */
    INTERSECTION =1,                                                    /**< Intersection between two objects */
    SUBTRACTION =2                                                      /**< Substract object two from object one */
};

/*!
 * @ingroup levelsetEnums
 * Enum class defining the intersection between a cell and the zero levelset.
 * For details see LevelSetObject::intersectSurface()
 */
enum class LevelSetIntersectionStatus{
    FALSE=0,                                                            /**< Cell does not intersect zero levelset */
    TRUE=1,                                                             /**< Cell does intersect zero levelset */
    CLOSE=2                                                             /**< Zero levelset lies within tangent and bounding radius of a cell */
};

/*!
 * @ingroup levelsetEnums
 * Enum class describing the how the intersection between cell and zero levelset 
 * should be computed
 */
enum class LevelSetIntersectionMode{
    FAST_FUZZY=0,                   /**< Compares levelset value to tangent and bounding radius of a cell */
    FAST_GUARANTEE_TRUE=1,          /**< All LevelSetIntersectionStatus::TRUE are accurate but LevelSetIntersectionStatus::FALSE may be wrong */
    FAST_GUARANTEE_FALSE=2,         /**< All LevelSetIntersectionStatus::FALSE are accurate but LevelSetIntersectionStatus::TRUE may be wrong */
    ACCURATE=3                      /**< Accurate but more costly checks */
};

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible fill-in modes for the levelset.
 */
enum class LevelSetFillIn{
    SPARSE,                         /**< Sparse fill-in, to be used when the levelset will be evaluated on a small portion of the domain */
    DENSE,                          /**< Dense fill-in, to be used when the levelset will be evaluated on almost all the domain */
};

typedef LevelSetFillIn LevelSetCacheType;

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible modes for the caches.
 */
enum class LevelSetCacheMode{
    BEGIN = 0,
    NONE = 0,              //!< No caching will be performed
    ON_DEMAND,             //!< Data are cached only where explicitly evaluated
    NARROW_BAND,           //!< Data are cached only inside the narrow band
    FULL,                  //!< Data are cached in the whole domain
    END,
    COUNT = END - BEGIN,
};

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible evaluation modes for cell data in the bulk.
 */
enum class LevelSetBulkEvaluationMode{
    BEGIN = 0,
    NONE = 0,              //!< No data is evaluated
    SIGN_PROPAGATION,      //!< Sign is propagated from the narrow band, no other data will be evaluated
    EXACT,                 //!< Exact data is evaluated
    END,
    COUNT = END - BEGIN,
};

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible level set fields
 */
enum class LevelSetField{
    BEGIN = 0,
    VALUE = BEGIN,                  /**< level set value */
    SIGN,                           /**< level set sign */
    GRADIENT,                       /**< level set gradient */
    SUPPORT,                        /**< facet that contains the projection point */
    PART,                           /**< part identifier at projection point */
    NORMAL,                         /**< body normal at projection point */
    END,
    COUNT = END - BEGIN,
    UNDEFINED
};

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible order of surface smooting
 */
enum class LevelSetSurfaceSmoothing{
    LOW_ORDER,                      /**< low order surface smoothing */
};

/*!
 * Hasher for the LevelSetField enum.
 */
struct LevelSetFieldHasher
{
    std::size_t operator()(const LevelSetField &field) const
    {
        return static_cast<int>(field);
    }
};

/*!
 * Set of field
 */
typedef std::vector<LevelSetField> LevelSetFieldset;

/*!
 * Map of write fields.
 */
template<typename value_t>
using LevelSetFieldMap = std::unordered_map<LevelSetField, value_t, LevelSetFieldHasher>;

/*!
 * @ingroup levelsetEnums
 * Enum class containing the possible fields to be added to the VTK file
 */
enum class LevelSetWriteField{
    VALUE    = static_cast<int>(LevelSetField::VALUE),     /**< adds level set value to VTK*/
    SIGN     = static_cast<int>(LevelSetField::SIGN),      /**< adds level set sign to VTK*/
    GRADIENT = static_cast<int>(LevelSetField::GRADIENT),  /**< adds level set gradient to VTK*/
    SUPPORT  = static_cast<int>(LevelSetField::SUPPORT),   /**< adds facet that contains the projection point to VTK*/
    PART     = static_cast<int>(LevelSetField::PART),      /**< adds part identifier at projection point to VTK*/
    NORMAL   = static_cast<int>(LevelSetField::NORMAL),    /**< adds body normal at projection point to VTK*/
    ALL,                                                   /**< adds level set value, gradient, normal and projection point to VTK*/
    DEFAULT                                                /**< adds levelset value and gradient to VTK*/
};

}

#endif
