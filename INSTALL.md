# BITPIT installation

bitpit runs on Linux and Mac OSX platforms.

## Dependencies
bitpit depends on
* c++ compiler supporting `-std=c++11`. It has been tested with g++ >= 4.7.3
* cmake >= 2.8
* (optionally) MPI implementation. It has been tested with OpenMPI >= 1.6.5.

Some additional dependencies are required for building specific modules
* libxml2 and its development headers are needed when compiling the 'IO'
  module (please note that the 'IO' module is a dependecies for many other
  bitpit modules, the only modules that do not depend on 'IO' are the low
  level modules like 'operators', 'containers', 'LA', and 'SA');
* blas, lapack, and lapacke are needed when compiling 'CG', 'RBF', and 'POD'
  modules.

When compiling bitpit, both shared library files and header library files are
required. Therefore, in addition to the library packages, also the corresponding
'*-devel' packages ('*-dev' in Debian-based distributions) are required.

## Confguring BITPIT
bitpit uses ccmake as building tool.
In the bitpit's root folder make a building folder, e.g. build
```bash
    bitpit$ mkdir build
```
Enter the `build` folder
```bash
    bitpit$ cd build
```
 In order to configure it, run:
```bash
    bitpit/build$ ccmake ../
```
 By this way, bitpit can be configured for production and installation.
Setting some variable in ccmake interface you can customize a bit your configuration.

The `CMAKE_BUILD_TYPE` variable has to be used to set the type of build. The possible options are : `None`, the environment compiler flags are used; `Release`, using compiler optimization flag `-O2`; `Debug`, related to compiler flags `-O0 -fmessage-length=0`, `RelWithDebInfo`, that uses compilation flags `-O2 -g` and `MinSizeRel` to have the smallest binary size.

In addition the `ENABLE_PROFILING` variable can be set to `ON` in order to add profiling flag `-pg` during the compilation.

The `ENABLE_MPI` variable can be used to compile the parallel implementation of the bitpit packages and to allow the dependency on MPI libraries.

The `BUILD_EXAMPLES` can be used to compile examples sources in `bitpit/examples`. Note that the tests sources in `bitpit/test`are necessarily compiled and successively available at `bitpit/build/test/` as well as the compiled examples are available at `bitpit/build/examples/`.

The module variables (available in the advanced mode) can be used to compile each module singularly by setting the related varible `ON/OFF` (BITPIT_MODULE_CONTAINERS, BITPIT_MODULE_IO, BITPIT_MODULE_LA, BITPIT_MODULE_SA...). Possible dependencies between bitpit modules are automatically resolved. 

Finally, you can choose the installation folder setting the cmake variable `CMAKE_INSTALL_PREFIX`. The default installation folder is `/usr/local/`.

Remember that if you choose the default installation path or another path without write permission you will need administration privileges to install bitpit in.

## Building and Installing
Once cmake has configured bitpit.s building just do
```bash
    bitpit/build$ make   
```
to build and
```bash
    bitpit/build$ make install   
```
to install.

If you have just built bitpit, its headers will be available at `bitpit/include/` folder and a static library `libbitpit.a` (or `libbitpit_MPI.a` in case of parallel compilation) will be available at `bitpit/build/lib/` folder.

If you have also installed bitpit, its headers will be available at `/my/installation/folder/bitpit/include/` folder and a static library `libbitpit.a` will be available at `/my/installation/folder/lib/` folder.

## Building Documentation
In order to build properly the documentation Doxygen (>=1.8.6) and Graphviz (>=2.20.2) are needed.

In the ccmake interface the variable `BUILD_DOCUMENTATION` can be set to `ON` in order to build the documentation during the library compilation. 
If turned on the new variable `DOC_EXTRACT_PRIVATE` can be used to include all the private class members in the documentation.
  
After the `make` or `make install` the doxygen documentation will be built. You can chose to compile only the documentation with command 
```bash
    bitpit/build$ make doc   
```
You can now browse the html documentation with your favorite browser by opening 'html/index.html'.

## Help
For any problem, please contact <a href="http://www.optimad.it">Optimad engineering srl</a> at info@optimad.it. 
