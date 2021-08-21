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

# ifndef __BITPIT_LEVELSET_OBJECT_FACTORY_TPP__
# define __BITPIT_LEVELSET_OBJECT_FACTORY_TPP__

namespace bitpit {

/*!
	@class      LevelSetObjectFactory
	@ingroup    levelset
	@brief      Allows to create levelset ojbects.
*/

/*!
 * Create a new boolean object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createBooleanObject(const LevelSetKernel *kernel, Args&&... args)
{
    BITPIT_UNUSED(kernel);

    return std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject(std::forward<Args>(args)...));
}

/*!
 * Create a new immutable object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createImmutableObject(const LevelSetKernel *kernel, Args&&... args)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createImmutableObject(cartesianKernel, std::forward<Args>(args)...);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createImmutableObject(octreeKernel, std::forward<Args>(args)...);
    }

    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new immutable object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param object is the source for the immutable object
 */
template<typename kernel_t>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createImmutableObject(const kernel_t *kernel, LevelSetObject *object)
{
    BITPIT_UNUSED(kernel);

    if (LevelSetBooleanObject *booleanObject = dynamic_cast<LevelSetBooleanObject *>(object)) {
        return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject(booleanObject));
    } else if (LevelSetCachedObject *cachedObject = dynamic_cast<LevelSetCachedObject *>(object)) {
        return std::unique_ptr<LevelSetObject>(new LevelSetImmutableObject(cachedObject));
    }

    BITPIT_UNREACHABLE("Unable to create an immutable object from the specified object.");
}

/*!
 * Create a new mask object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createMaskObject(const LevelSetKernel *kernel, Args&&... args)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createMaskObject(cartesianKernel, std::forward<Args>(args)...);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createMaskObject(octreeKernel, std::forward<Args>(args)...);
    }

    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new mask object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename kernel_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createMaskObject(const kernel_t *kernel, Args&&... args)
{
    BITPIT_UNUSED(kernel);

    return std::unique_ptr<LevelSetObject>(new LevelSetMaskObject(std::forward<Args>(args)...));
}

/*!
 * Create a new segmentation object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::createSegmentationObject(const LevelSetKernel *kernel, Args&&... args)
{
    if( const LevelSetCartesianKernel *cartesianKernel = dynamic_cast<const LevelSetCartesianKernel *>(kernel) ){
        return _createSegmentationObject(cartesianKernel, std::forward<Args>(args)...);
    } else if ( const LevelSetOctreeKernel *octreeKernel = dynamic_cast<const LevelSetOctreeKernel *>(kernel) ){
        return _createSegmentationObject(octreeKernel, std::forward<Args>(args)...);
    }

    BITPIT_UNREACHABLE("Kernel type not supported");
}

/*!
 * Create a new segmentation object for the specified kernel.
 *
 * \param kernel is the kernel
 * \param args are the arguments that will be forwarded to the constructor
 */
template<typename kernel_t, typename... Args>
std::unique_ptr<LevelSetObject> LevelSetObjectFactory::_createSegmentationObject(const kernel_t *kernel, Args&&... args)
{
    BITPIT_UNUSED(kernel);

    return std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(std::forward<Args>(args)...));
}

}

# endif
