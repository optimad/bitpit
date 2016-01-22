/*!
 * @defgroup IO Input/Output
 * @{
 * @defgroup Generic Generic
 * @defgroup DuneGridFormat DuneGridFormat (DGF)
 * @defgroup STereoLithography STereoLithography (STL)
 * @defgroup VisualizationToolKit VisualizationToolKit (VTK)
 * @}
 *
 */

#ifndef __BITPIT_IO__
#define __BITPIT_IO__

#include "Class_FH.hpp"
#include "Class_VTK.hpp"
#ifdef IO_ENABLE_VTK_WRAPPERS
#include "Class_VTK_Wrappers.hpp"
#endif
#include "DGF_IOFunct.hpp"
#include "Generic_IO.hpp"
#include "STL_IOFunct.hpp"
#include "VTK_IOFunct.hpp"

#endif
