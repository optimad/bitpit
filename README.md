PABLO
=====

PArallel Balanced Linear Octree

=====

PABLO is C++/MPI library for parallel linear octree/quadtree developed by Optimad Engineering srl under the GNU Lesser General Public License. The aim of the project is to provide users
with a ready-to-use tool for parallel adaptive grid of quadrilaterals/hexahedra.
Transparent MPI communications are implemented in PABLO for dynamic load-balancing of both grid and user data. Dynamic load-balancing is possible and easy to perform with simple calls to methods provided with PABLO.
PABLO allows adaptive mesh refinement by generating non-conforming grid with hanging nodes. One of the main feature of PABLO is the low memory consumption in the basic configuration (approx. 30B per octant in 3D).
Additional features available in PABLO are:  2:1 balancing between octants and a easy way to generate and store intersections between octants.