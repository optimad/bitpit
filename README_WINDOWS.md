#  bitpit on Microsoft Windows

This is a short guide to set a bitpit compliant 64bit Windows environment using MSYS2/MingGW64.

## Requirements:
- <B>Windows 8/10</B> (8.1 pro/10 tested)

- <B>MSYS2</B> (64bit) https://www.msys2.org/

- <B>MSYS2 packages</B>:
  - base-devel
  - binutils
  - python2
  - git (optional)


- <B>MinGW64 packages</B>:

  - mingw-w64-x86_64-toolchain
  - mingw-w64-x86_64-lapack
  - mingw-w64-x86_64-msmpi
  - mingw-w64-x86_64-libxml2
  - mingw-w64-x86_64-cmake
  - mingw-w64-x86_64-qt5 (optional, but needed to run cmake-gui. Consider 1.22 GB to install it.)
  - mingw-w64-x86_64-doxygen (to generate Doxygen documentation)
  - mingw-w64-x86_64-graphviz (to generate Doxygen documentation)


- <B>Microsoft MPI</B>: MSMpiSetup.exe or msmpisdk.msi (dowloadable for free from https://www.microsoft.com/ searching for "MicrosoftMPI v x.x.x". Here the exact version x.x.x must be compliant with the version of MinGW64 package mingw-w64-x86_64-msmpi. See MSMPI section in          procedure chapter for details)

- <B>PETSc library</B> https://www.mcs.anl.gov/petsc/download/index.html. Version 3.12.5 tested (older version 3.8.4 tested too).

## Procedure
#### Install MSYS2 (64 bit)/ MinGW64

1. Download MSYS2 from https://www.msys2.org/ and install it following the instructions on the site to update and get the environment ready to be used. A more detailed install guide can be found also here: https://www.msys2.org/wiki/MSYS2-installation/. This will install 2 subsystems in particular:
    - *msys2*: general POSIX compliant emulating environment
    - *mingw64*: environment, based on *msys2* and specifically tuned to provide better interoperability with native 64bit Windows software, using MingGW.


2. If you did not run MSYS2 just after installation, look for msys2 or mingw64 shells among your installed programs or launch them from cmd (Windows Prompt):
    ```bash
    >C:\msys64\msys2_shell.cmd -msys2 (for msys2 shell)
    >C:\msys64\msys2_shell.cmd -mingw64 (for mingw64 shell)
    ```

3. Open a *msys2* shell and install base *msys2* packages base-devel, binutils, python2, git. On usage of *pacman* please refer to https://www.msys2.org/wiki/Using-packages/ :
    ```bash
    user@machine MSYS2 ~
    > pacman -S base-devel binutils python2 git
    ```
    As in a Linux system, libs and binaries installed can be found navigating C:/msys64/usr/lib or C:/msys64/usr/bin

4. Open a *mingw64* shell and install the base development toolchain for *mingw64*. This will include C,C++, Fortran compilers and libraries :
    ```bash
    user@machine MINGW64 ~
    > pacman -S mingw-w64-x86_64-toolchain
    ```
    As in a Linux system, libs and binaries installed can be found navigating C:/msys64/mingw64/lib or C:/msys64/mingw64/bin. Verify the installation with:
    ```bash
    user@machine MINGW64 ~
    > gcc -v
    ```
    It should print something like this: (the gcc version will depend on the update status of your MSYS2 package repository. Please refer to MSYS2 documentation for instruction on how update it)
    ```
    Using built-in specs.
    COLLECT_GCC=C:\msys64\mingw64\bin\gcc.exe
    COLLECT_LTO_WRAPPER=C:/msys64/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/9.3.0/lto-wrapper.exe
    Target: x86_64-w64-mingw32
    Configured with: ../gcc-9.3.0/configure --prefix=/mingw64 --with-local-prefix=/mingw64/local --build=x86_64-w64-mingw32 --  host=x86_64-w64-mingw32 --target=x86_64-w64-mingw32 --with-native-system-header-dir=/mingw64/x86_64-w64-mingw32/include --libexecdir=/mingw64/lib --enable-bootstrap --with-arch=x86-64 --with-tune=generic --enable-languages=c,lto,c++,fortran,ada,objc,obj-c++ --enable-shared --enable-static --enable-libatomic --enable-threads=posix --enable-graphite --enable-fully-dynamic-string --enable- libstdcxx-filesystem-ts=yes --enable-libstdcxx-time=yes --disable-libstdcxx-pch --disable-libstdcxx-debug --disable-isl-version-check --enable-lto --enable-libgomp --disable-multilib --enable-checking=release --disable-rpath --disable-win32-registry --disable-nls --  disable-werror --disable-symvers --enable-plugin --with-libiconv --with-system-zlib --with-gmp=/mingw64 --with-mpfr=/mingw64 --with- mpc=/mingw64 --with-isl=/mingw64 --with-pkgversion='Rev2, Built by MSYS2 project' --with-bugurl=https://sourceforge.net/projects/msys2 - -with-gnu-as --with-gnu-ld
    Thread model: posix
    gcc version 9.3.0 (Rev2, Built by MSYS2 project)
    ```

5. In order to make your Windows OS *see* the *mingw64* shell binaries, modify the Windows  *PATH* environment variable by adding `C:\msys64\mingw64\bin`. The same output of 4) should appear launching from a Windows PowerShell or a Windows prompt the command:
    ```bash
    > gcc.exe -v
    ```
    If not, check if `C:\msys64\mingw64\bin` is present in windows PATH typing:
    ```bash
    (prompt)     > echo %PATH%
    (PowerShell) > $Env:Path
    ```

#### MS-MPI:

This section will install MPI on your machine, both libraries and compiler wrappers (*mpic, mpicxx*) through mingw-w64-x86_64-msmpi package and *mpiexec* command through native Windows package MicrosoftMPI

1. Open a *mingw64* shell and type:
    ```bash
    user@machine MINGW64 ~
    > pacman -Ss mingw-w64-x86_64-msmpi
    ```
    It should return something like:
    ```bash
    mingw64/mingw-w64-x86_64-msmpi 10.1.1-1
    Microsoft MPI SDK (mingw-w64)
    ```
    Please take note of the MSMPIVersion number, in the specific case 10.1.1, but it can vary according to the update status of your MSYS2 package repository. Install the msmpi from *mingw64*:
    ```bash
    user@machine MINGW64 ~
    > pacman -S mingw-w64-x86_64-msmpi
    ```
    Close the *mingw64* shell.

2. Download Microsoft MPI from https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi. If on site version does not match the MSMPIVersion noted earlier, a simple browser search of "MSMPI" and target version number will redirect you to the correct version download page on the official Microsoft site. Then install **both** *MSMpiSetup.exe* and *msmpisdk.msi*.

3. Reopen the *mingw64* shell (to load the new env variables introduced in 2.) and locate the environment variable MSMPI_BIN typing:
    ```bash
    user@machine MINGW64 ~
    > printenv | grep "WIN\|MSMPI"
    ```
    You should see something like this:
    ```
    USERDOMAIN=VIRTUALWIN8
    MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\
    COMPUTERNAME=VIRTUALWIN8
    USERDOMAIN_ROAMINGPROFILE=VIRTUALWIN8
    MSMPI_BIN=C:\Program Files\Microsoft MPI\Bin\
    MSMPI_LIB32=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\
    MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\
    WINDIR=C:\Windows
    ```
    If not, please verify you MSMPI installation on 2. Take note of MSMPI_BIN content path (here it is `C:\Program Files\Microsoft MPI\Bin\`).

4. In order to use *mpiexec.exe* inside the *mingw64* shell, the content path of MSMPI_BIN needs to be added to the local PATH of the shell. Sometimes problems of compatibility can arise in a direct export, because of blank spaces inside native Windows path. In such cases is strongly recommended to use **cygpath** (already available in the shell) to convert this path in a more compact and unix compliant format. To do so, it is sufficient to transform the original Windows path in dos short format:
    ```bash
    user@machine MINGW64 ~
    > cygpath -ms 'C:\Program Files\Microsoft MPI\Bin\'
    result will be: C:/PROGRA~1/MICROS~2/Bin/
    ```
    then transform the last result in a unix compliant path:
    ```bash
    user@machine MINGW64 ~
    > cygpath -u 'C:/PROGRA~1/MICROS~2/Bin/'
    result will be: /c/PROGRA~1/MICROS~2/Bin/
    ```
    The last result can be safely exported in the local *mingw64* shell PATH environment	variable with:
    ```bash
    export PATH=$PATH:/c/PROGRA~1/MICROS~2/Bin/
    ```
    and make *mpiexec.exe* available in the shell. It is recommended to add it permanently to the *~/.bashrc* conf file of the shell.

5. Before moving on, test MSMPI with a simple MPI Hello World code. Here an example in c++:
    ```c++
    #include <mpi.h>
    #include <iostream>
    int main(int argc,char** argv){
            MPI_Init(&argc,&argv);
            int rank;
            int nprocs;
            MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
            std::cout << "I'm " << rank << " of " << nprocs << std::endl;
            MPI_Finalize();
    }
    ```
    compile it using MSMPI language relative wrapper
    ```bash
    user@machine MINGW64 ~
    > mpicxx -o mpi_hello mpi_hello.cpp
    ```
    and launch your test sample:
    ```bash
    user@machine MINGW64 ~
    >  mpiexec -n 2 ./mpi_hello
    ```
    You should get
    ```
    I'm 0 of 2
    I'm 1 of 2
    ```

#### bitpit dependencies
For the following subsections, open a *mingw64* shell.

_**Beware**_: some of these deps can be available also as native Windows applications/libraries. Anyway we choose to use those provided by msys2/mingw64 repositories to ensure the maximum coherence during the installation flow. If you choose to use different versions/distributions, let us know your experience. We will be glad to use your tips to enrich the current guide.

**__CMAKE__**

Install *cmake* and its GUI *cmake-gui*  with:
```bash
user@machine MINGW64 ~
> pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-qt5
```
If you are short in disk space, you can avoid installing *mingw-w64-x86_64-qt5* (~1.22 GB), but you cannot use cmake-gui for setting cmake variables in building configuration. Without gui, you have to consider a large use of cmake "*-D*" flag or to manually modify *CMakeCache.txt* file, in order to set the values of the cmake variables.

**__LAPACK__**

Install *lapack* with:
```bash
user@machine MINGW64 ~
> pacman -S mingw-w64-x86_64-lapack
```

**__LIBMXL2__**

Install *libxml2* with:
```bash
user@machine MINGW64 ~
> pacman -S mingw-w64-x86_64-libxml2
```

**__PETSc__**

First, __*PETSc* sources__ are needed. Download *PETSc v3.12.5* from https://www.mcs.anl.gov/petsc/download/index.html. This is the most recent version of PETSc tested with bitpit. If you want to have a try with different versions, please let us know your experience.

- Untar the downloaded *PETSc* archive and enter the *PETSc* sources folder from the *mingw64* shell.

- Get track of *msys2* python2.exe and ar.exe binaries, both available in '/usr/bin/' inside the *mingw64* shell (if not check subsection 3 of **Install MSYS2(64bit) and MINGW64** chapter of the guide). These ones are the only that work properly with *PETSc* configure scripts. *mingw64* python and ar versions create problems.

- Check if *mpiexec.exe* is available in the shell (type "mpiex" and complete with tab). If not please check again  the section 4 of **MS-MPI** chapter.

- We configure *PETSc* using a *Python2* script. Edit a new *Python* file, e.g. __*myConf.py*__, and add these lines:
    ```python
     configure_options = [
     '--with-ar=/usr/bin/ar' ,
     '--with-shared-libraries=0',
     '--with-debugging=0',
     '--with-visibility=0',
     '--prefix=/c/msys64/mingw64/petsc',
     'FOPTFLAGS=-O3 -fno-range-check',
     'COPTFLAGS=-O3',
     'CXXOPTFLAGS=-O3'
     ]

     if __name__ == '__main__':
        import sys,os
        sys.path.insert(0,os.path.abspath('config'))
        import configure
        configure.petsc_configure(configure_options)
    ```
    The prefix can be any path on your Windows machine, provided that it is specified as unix compliant - blank space free path (use **cygpath**, an example is reported in **MS-MPI** chapter, section 4). You can customize your *PETSc* configuration following the instructions at https://www.mcs.anl.gov/petsc/documentation/installation.html. The *Python* script above gives you an optimized installation of *PETSc* library that works with *bitpit*.

- Therefore, launch:
    ```bash
    user@machine MINGW64 ~
    > /usr/bin/python2 myConf.py
    ```
    and follow the PETSc instructions to complete the installation ( copy and paste exactly the command lines PETSc suggests, up to the 'make...test' part).

- Finally, export __*PETSC_DIR*__ and __*PETSC_ARCH*__ environment variables on    *mingw64* local shell by adding  to *~/.bashrc* file : `export PETSC_DIR=#prefix# PETSC_ARCH=""`, where #prefix# is the path specified with --prefix in the python configuration script.


#### bitpit
- Download *bitpit* master archive at https://github.com/optimad/bitpit/archive/master.zip or git clone it using SSH or HTTPS.

- Once the *MinGW64* development environment is ready, you can install __*bitpit*__ following the  installation instructions you can find at https://github.com/optimad/bitpit/blob/master/INSTALL.md or locally in the *INSTALL.md* file in *bitpit* root folder.

- In order to compile a Windows native version of *bitpit* you have to:

  - specify the __*MinGW Makefiles Generator*__ to cmake, e.g. in your bitpit build folder
  ```bash
  user@machine MINGW64 ~
  > cmake -G "MinGW Makefiles" ..
  ```
  or if *cmake-gui* is used, enable a __*MinGW Makefiles*__ project.

  - use __*mingw32-make*__ in place of standard *make* command to build *bitpit*.

- **__Warning__ 1**: if *bitpit cmake* configuration exits with PETSc error, verify the values of *PETSC_DIR* and *PETSC_ARCH* *cmake* variables and possibly set them to the correct values, i.e. the values the relative environment variables have. If you used exactly *myConf.py* script above to configure *PETSc*, their values are
   ```bash
   PETSC_DIR=/c/msys64/mingw64/petsc
   PETSC_ARCH=
   ```
