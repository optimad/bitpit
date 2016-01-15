# Change Log
All notable changes to this project will be documented in this file.
This project *tries* to adhere to [Semantic Versioning](http://semver.org/).

# v1.0.0 - 2016-01-13

## COMMON

## PABLO

### Added
- Unique Octant 2D/3D class, no more template class.
- Unique LocalTree 2D/3D class, no more template class.
- Unique ParaTree 2D/3D classes, no more template class.
- New basic get methods for parallel and local tree. 

### Changed
- Major changes in ParaTree interface. Main change: use of array instead of vector where possible.
- Class Global rearranged and new methods to manipulate it.
- Coding style and syntax modified according to the Coding Style Guide of this project.
- Pre-compilation variable NOMPI renamed ENABLE_MPI.
- CMakeLists files changed with new variables ENABLE_MPI.

### Fixed
- Bug-fixing in adapt tracking the changes with the mapper. In this version, when the mapping is active, even if some markers are set <-1 or >1, the octants will be adapted one times in refinement and coarsening; the resulting markers will be decreased or increased by 1 respctively.

### Removed
- Old deprecated use of mapping now removed.
- Private overloaded adapt methods removed.

## PATCHMAN

## SURFTRI

## UCARTMESH

## VOLTRI

# v0.4.0 - 2015-11-10

## COMMON

## PABLO

### Added
- This CHANGELOG file.
- CMAKE folders, files and new custom variables for ccmake precompiling.
- make check for testing serial and parallel version of PABLO.
- Active option for building documentation in ccmake setup.

### Changed
- Installation and building instructions.

### Fixed
- Include instruction in Class_Log.

### Removed
- Temporary tests for bug fixing.

## PATCHMAN

## SURFTRI

## UCARTMESH

## VOLTRI


