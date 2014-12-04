PABLO
=====

PArallel Balanced Linear Octree

=====

PABLO is C++/MPI library for parallel linear octree/quadtree developed by Optimad Engineering srl under the GNU Lesser General Public License. 

The aim of the project is to provide users with a ready-to-use tool for parallel adaptive grid of quadrilaterals/hexahedra. 

Message passing paradigm is transparent to the user since PABLO has embedded MPI calls. By this way, the user can easily perform data communications and dynamic load-balance by calling straightforward high level methods. 

Moreover, the user can feel free to customize his data in whatever way he likes. 

PABLO allows adaptive mesh refinement by generating non-conforming grid with hanging nodes. 

One of the main feature of PABLO is the low memory consumption in the basic configuration (approx. 30B per octant in 3D). 

Additional features available in PABLO are: 2:1 balancing between octants and a easy way to generate and store intersections between octants.

##Installation
Please, see INSTALL.md for installation instructions.

