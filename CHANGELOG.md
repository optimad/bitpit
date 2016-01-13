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

### Fixed
- Possibility to adapt the octree, tracking the changes with the mapper, even with marker <-1 and >1.

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


