# PABLO installation

PABLO runs on Linux and Mac OSX platforms.

## Dependencies
PABLO depends on
* c++ compiler supporting `-std=c++11`. It has been tested with g++ >= 4.7.3, clang++ >=3.5 and icpc >= 15.0
* cmake >= 2.8
* (optionally) MPI implementation. It has been tested with OpenMPI >= 1.6.5. 

## Confguring PABLO
PABLO uses cmake as building tool. The use of curses interface ccmake is reccommended.
In the PABLO's root folder make a building folder, e.g. build
```bash
	PABLO$ mkdir build
```
Enter the `build` folder
```bash
	PABLO$ cd build
```
 In order to configure it with default options, run:
```bash
	PABLO/build$ cmake ../
```
 By this way, PABLO is configured for production (using compiler optimization flags for release cmake building type, `-O3`), the test sources in `PABLO/test/` will be compiled and successively available at `PABLO/build/test/` while the examples sources in `PABLO/examples/` will not. 

The default installation folder is `/opt/bitpit/PABLO`. The folder `/opt/bitpit/` is the default installation folder for bitpit libraries developed by Optimad engineering and available on GitHub.

Setting some variable in ccmake you can customize a bit your configuration.

You can change the `BITPIT_DIR` to choose alternatives location for bitpit libraries. Note that if you choose the default bitpit installation path or another path without write permission you will need administration privileges to install bitpit libraries in.

You can even choose different installation folder for PABLO setting the cmake variable `CMAKE_INSTALL_PREFIX`.
Remember that if you choose the default installation path or another path without write permission you will need administration privileges to install PABLO in.

The `CMAKE_BUILD_TYPE` variable can be used to set the compiler flags for production/debugging, then you can set it to Release/Debug to obtain a desired version of PABLO. `CMAKE_BUILD_TYPE` has not a default value.

The `WITHOUT_MPI` variable can be used to compile the serial implementation of PABLO and to avoid the dependency on MPI libraries, then you can set to `ON` to obtain a serial version of PABLO. `WITHOUT_MPI` default value is `OFF`.

The `BUILD_EXAMPLES` variable can be use to compile examples, then set to `ON` and the building procedure will compile the examples sources. `BUILD_EXAMPLES` default value is `OFF`.

You can change the `COMPILER` variable and use `gcc` or `intel` option to compile PABLO with system primary compiler or forcing intel compiler.

The `ENABLE_PROFILING` variable can be use to activate `-Wall -Wextra` compilation flags.


## Building and Installing
Once cmake has configured PABLO's building just do
```bash
	PABLO/build$ make	
```
to build and
```bash
	PABLO/build$ make install	
```
to install.

If you have just built PABLO, its headers will be available at `PABLO/include/` folder and a static library `libPABLO.a` will be available at `PABLO/build/lib/` folder.

If you have also installed PABLO, its headers will be available at `/my/installation/folder/PABLO/include/` folder and a static library `libPABLO.a` will be available at `/my/installation/folder/lib/` folder.

To build and run the compilation tests do
```bash
	PABLO/build$ make check
```
and the tests will start immediately after the compilation of the library.
 

## Building Documentation
In order to build properly the documentation Doxygen (>=1.8.6) and Graphviz (>=2.20.2) are needed.
In `ccmake` setup set `BUILD_DOCUMENTATION` to `ON` to build Doxygen documentation.

You can now browse the html documentation with your favorite browser by opening `build/doc/html/index.html`.


## Help
For any problem, please join the <a href="https://groups.google.com/forum/#!forum/pablo-users" target="pablousers">PABLO Users Google Group</a> and post your requests.

