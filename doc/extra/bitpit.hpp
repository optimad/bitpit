/*! @mainpage bitpit

<B>bitpit</B> is a C++ library for scientific High Performance Computing.
Within bitpit different modules factorize the typical effort which is needed to derived a real-life application code.

Efforts is dedicated to handle different types of computational meshes, their runtime adaptation and data transfer for parallel applications.

# Basic modules

## common common
The <B>bitpit::utils</B> namespace provides miscellaneous functions for the bitpit framework

## operators
- <B> + - * / </B> operators for the std::vector and std::array classes
- basic mathematic functions like <B> dotProduct crossProduct norm abs ...</B>  for the std::vector and std::array classes
- <B> << >> </B> operators and <B> display </B> functions for the std::vector and std::array classes

## containers
- <B>bitpit::PiercedVector</B> is a container which allows the cancellation and insertion on the fly of elements.
- <B>bitpit::CollapsedArray2D & bitpit::CollapsedVector2D</B> are 2D linearized vectors which avoid the overhead of the vector infrastructures.
- <B>bitpit::IBinaryStream & bitpit::OBinaryStream</B> allow to copy chunks of memory

## Input Output (IO)
- <B>bitpit::genericIO</B> is a namespace which contains basic routines to read/write in ASCII or BINARY
- <B>bitpit::DGFObj</B> class and <B>bitpit::dgf</B> namespace for reading and writing Dune Grid Format files
- <B>bitpit::STLObj</B> class and <B>bitpit::stl</B> namespace for reading and writing Stereo Litography files
- <B>bitpit::VTK</B> classes and <B>bitpit::vtk</B> namespace for reading and writing Visualization ToolKit files
- <B>bitpit::Log</B> classes and <B>bitpit::log</B> namespace for unified and coordinated output throughout all modules

## Linear Algebra (LA)
<B>LA</B> providesi methods for small dense linear systems stored as std::vector<std::vector> or std::array<std::array>
- creation of basic matrices, like identity and diagonal matrices
- basic matrix operations like transpose and matrix multiplication
- solution 

## Sort algorithms (SA)
- <B>bitpit::LIFOStack</B> manage insertion and extraction uding the Last-In-First-Out pragma
- <B>bitpit::KdTree</B> sorts vertices in a d-dimensional Euclidean space in a Kd-tree structure
- <B>bitpit::MinPQueue & bitpit::MaxPQueue</B> manage a binary tree which ensures that the root element has the smalles/largest value in the tree

## Computational Geometry (CG)
<B>CG</B> module provides methods for computational geometry.

# Mesh modules

## PArallel Balanced Linear Octree (PABLO)
<B>PABLO</B> is a stand-alone module which provides a parallel linear octree/quadtree. 
<B>bitpit::ParaTree</B> provides connectivity/adjecency information only,
whereas <B>bitpit::PabloUniform</B> provides aditionally all geometrical information for an Octree within and rectangular domain.
Message passing paradigm is transparent to the user since PABLO has embedded MPI calls. By this way, the
user can easily perform data communications and dynamic load-balance by calling
straightforward high level methods.

PABLO allows adaptive mesh refinement by generating non-conforming grid with
hanging nodes.
Additional features available in PABLO are: 2:1 balancing between octants
and a easy way to generate and store intersections between octants.

## Patch kernel 
<B>bitpit::PatchKernel</B> is the base mesh container of bitpit.
Basic elements like <B>bitpit::Vertex</B>, <B>bitpit::Interface</B> and <B>bitpit::Cell</B> are defined here, together with <B>bitpit::Adaption</B> which is used for dynamic mesh adaptation. 
It provides a homogenous interface class to all types of meshes and two specialized derived classes, <B>bitpit::SurfaceKernel</B> and <B>bitpit::VolumeKernel</B>, for surface and volume meshes.

## Unstructured surface patch 
<B>bitpit::SurfUnstructured</B> is the principal container for surface segmentations and has methods in order to read/write surface triangulations.

## Cartesian/Octree/Unstructured volume patch 
<B>bitpit::VolCartesian, bitpit::VolOctree & bitpit::VolUnstructured</B> are the derived voume meshes ifor 2D and 3D in bitpit.
They share the common interface through <B>bitpit::VolumeKernel</B> but each grid provides specific optimized methods.
*/
