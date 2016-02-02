/*!

\mainpage bitpit

<B>bitpit</B> is a C++ library for scientific high-performance computing.
<B>bitpit</B> is composed of several modules:

Basic modules
========================

1. common
------------------------
<B>common</B> provides generic functions for the bitpit framework

2. operators
------------------------
- <B> + - * / </B> operators for the std::vector and std::array classes
- basic mathematic functions like <B> dotProduct crossProduct norm abs ...</B>  for the std::vector and std::array classes
- <B> << >> </B> operators and <B> display </B> functions for the std::vector and std::array classes
- basic functions to manipulate std::strings

3. containers
------------------------
- <B>PiercedVector</B> is a container which allows the cancellation and insertion on the fly of elements.
- <B>CollapsedArray2D & CollapsedVector2D</B> are 2D linearized vectors which avoid the overhead of the vector infrastructures.
- <B>I/O BinaryStream</B> allow to copy chunks of memory

4. Input Output (IO)
------------------------
- <B>GenericIO</B> basic routines to read/write in ASCII or BINARY
- <B>DGF</B> classes and routines for reading and writing Dune Grid Format files
- <B>STL</B> classes and routines for reading and writing Stereo Litography files
- <B>VTK</B> classes and routines for reading and writing Visualization ToolKit files

5. Linear Algebra (LA)
------------------------
<B>LA</B> providesi methods for small dense linear systems stored as std::vector<std::vector> or std::array<std::array>
- creation of basic matrices, like identity and diagonal matrices
- basic matrix operations like transpose and matrix multiplication
- solution 


6. Sort algorithms (SA)
------------------------
- <B>LIFOStack</B> manage insertion and extraction uding the Last-In-First-Out pragma
- <B>KdTree</B> sorts vertices in a d-dimensional Euclidean space in a Kd-tree structure
- <B>MinPQueue & MaxPQueue</B> manage a binary tree which ensures that the root element has the smalles/largest value in the tree

Mesh modules
========================

1. PArallel Balanced Linear Octree (PABLO)
------------------------
<B>PABLO</B> provides a parallel linear octree/quadtree. Message passing paradigm
is transparent to the user since PABLO has embedded MPI calls. By this way, the
user can easily perform data communications and dynamic load-balance by calling
straightforward high level methods.

PABLO allows adaptive mesh refinement by generating non-conforming grid with
hanging nodes.

Additional features available in PABLO are: 2:1 balancing between octants
and a easy way to generate and store intersections between octants.


*/
