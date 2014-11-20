# PABLO installation

PABLO runs on Linux and Mac OSX platforms.

## Dependencies
PABLO depends on
* c++ compiler supporting `-std=c++11`. It has been tested with g++ >= 4.7.3 and clang++ >=3.5
* cmake >= 2.8
* MPI implementation. It has been tested with OpenMPI >= 1.6.5. 

## Confguring PABLO
PABLO uses cmake as building tool.
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
 By this way, PABLO is configured for production (using compiler optimization flags, `-O3`) and the test sources in `PABLO/test/` will be compiled and successively available at `PABLO/build/test/`; moreover, the default installation folder is `/usr/local/`.

Passing some variable to cmake you can customize a bit your configuration.

The `DEBUG` variable can be used to set the compiler flags `-ggdb -O0 -fmessage-length=0`, then you can call
```bash
	PABLO/build$ cmake -DDEBUG=1 ../	
```
to obtain a debug version of PABLO. `DEBUG` default value is 0.

The `COMPILE_TESTS` variable can be use to avoid tests compilation, then
```bash
	PABLO/build$ cmake -DCOMPILE_TESTS=0 ../	
```
and the building procedure will not compile the test sources. `COMPILE_TESTS` default value is 1.

Finally, you can choose the installation folder setting the cmake variable `CMAKE_INSTALL_PREFIX`
```bash
	PABLO/build$ cmake -DCMAKE_INSTALL_PREFIX=/my/installation/folder/ ../	
```
Remember that if you choose the default installation path or another path without write permission you will need administration privileges to install PABLO in.

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

## Help
For any problem, please join the [PABLO Users Google Group](https://groups.google.com/forum/#!forum/pablo-users) and post your requests.