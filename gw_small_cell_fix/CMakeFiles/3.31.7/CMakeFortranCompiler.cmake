set(CMAKE_Fortran_COMPILER "/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "14.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/bin/gcc-ar")
set(CMAKE_RANLIB "/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/bin/ranlib")
set(CMAKE_LINKER "/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/bin/ld")
set(CMAKE_Fortran_COMPILER_LINKER "/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/bin/ld")
set(CMAKE_Fortran_COMPILER_LINKER_ID "GNU")
set(CMAKE_Fortran_COMPILER_LINKER_VERSION 2.42)
set(CMAKE_Fortran_COMPILER_LINKER_FRONTEND_VARIANT GNU)
set(CMAKE_Fortran_COMPILER_RANLIB "/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/bin/gcc-ranlib")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED TRUE)
set(CMAKE_LINKER_PUSHPOP_STATE_SUPPORTED TRUE)
set(CMAKE_Fortran_LINKER_PUSHPOP_STATE_SUPPORTED TRUE)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib/gcc/x86_64-pc-linux-gnu/14.2.0/finclude;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/dbcsr-2.8.0/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/sirius-7.7.0/include/sirius;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/SpLA-1.6.1/include/spla;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/SpFFT-1.1.1/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/pugixml-1.15/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/spglib-2.5.0/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libvdwxc-0.4.0/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/hdf5-1.14.6/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/COSMA-2.7.0/include;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libxc-7.0.0/include;/opt/software/pc2/EB-SW/software/Python/3.13.1-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libffi/3.4.5-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/SQLite/3.47.2-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/Tcl/8.6.16-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libreadline/8.2-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/GSL/2.8-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/libarchive/3.7.7-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/zstd/1.5.6-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/lz4/1.10.0-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/cURL/8.11.1-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libpsl/0.21.5-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libunistring/1.3-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libiconv/1.18-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libidn2/2.3.7-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/bzip2/1.0.8-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/ncurses/6.5-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/FFTW.MPI/3.3.10-gompi-2025a/include;/opt/software/pc2/EB-SW/software/FFTW/3.3.10-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/FlexiBLAS/3.4.5-GCC-14.2.0/include/flexiblas;/opt/software/pc2/EB-SW/software/FlexiBLAS/3.4.5-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/AOCL-BLAS/5.1-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/OpenBLAS/0.3.29-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/OpenMPI/5.0.7-GCC-14.2.0/include;/opt/software/pc2/EB-SW/software/UCC/1.3.0-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/PRRTE/3.0.8-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/PMIx/5.0.6-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libfabric/2.0.0-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/UCX/1.18.0-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libevent/2.1.12-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/OpenSSL/3/include;/opt/software/pc2/EB-SW/software/hwloc/2.11.2-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libpciaccess/0.18.1-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/libxml2/2.13.4-GCCcore-14.2.0/include/libxml2;/opt/software/pc2/EB-SW/software/libxml2/2.13.4-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/XZ/5.6.3-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/numactl/2.0.19-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/zlib/1.3.1-GCCcore-14.2.0/include;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/include;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/lib64;/opt/software/pc2/EB-SW/software/Python/3.13.1-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libffi/3.4.5-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/SQLite/3.47.2-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/Tcl/8.6.16-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libreadline/8.2-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/GSL/2.8-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libarchive/3.7.7-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/zstd/1.5.6-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/lz4/1.10.0-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/cURL/8.11.1-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libpsl/0.21.5-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libunistring/1.3-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libiconv/1.18-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libidn2/2.3.7-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/bzip2/1.0.8-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/ncurses/6.5-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/ScaLAPACK/2.2.2-gompi-2025a-fb/lib64;/opt/software/pc2/EB-SW/software/FFTW.MPI/3.3.10-gompi-2025a/lib64;/opt/software/pc2/EB-SW/software/FFTW/3.3.10-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/FlexiBLAS/3.4.5-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/AOCL-BLAS/5.1-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/OpenBLAS/0.3.29-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/OpenMPI/5.0.7-GCC-14.2.0/lib64;/opt/software/pc2/EB-SW/software/UCC/1.3.0-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/PRRTE/3.0.8-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/PMIx/5.0.6-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libfabric/2.0.0-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/UCX/1.18.0-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libevent/2.1.12-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/OpenSSL/3/lib64;/opt/software/pc2/EB-SW/software/hwloc/2.11.2-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libpciaccess/0.18.1-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/libxml2/2.13.4-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/XZ/5.6.3-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/numactl/2.0.19-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/zlib/1.3.1-GCCcore-14.2.0/lib64;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib64;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib/gcc/x86_64-pc-linux-gnu/14.2.0;/lib64;/usr/lib64;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/dbcsr-2.8.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/sirius-7.7.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/SpLA-1.6.1/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/SpFFT-1.1.1/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/pugixml-1.15/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libvori-220621/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/spglib-2.5.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libvdwxc-0.4.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/hdf5-1.14.6/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/plumed-2.9.3/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/elpa-2024.05.001/cpu/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/COSMA-2.7.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libxsmm-e0c4a2389afba36c453233ad7de07bd92c715bec/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libxc-7.0.0/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/tools/toolchain/install/libint-v2.6.0-cp2k-lmax-7/lib;/opt/software/pc2/EB-SW/software/CP2K/2025.2-foss-2025a-gcc-openmpi-aocl/lib;/opt/software/pc2/EB-SW/software/Python/3.13.1-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libffi/3.4.5-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/SQLite/3.47.2-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/Tcl/8.6.16-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libreadline/8.2-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/GSL/2.8-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/libarchive/3.7.7-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/zstd/1.5.6-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/lz4/1.10.0-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/cURL/8.11.1-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libpsl/0.21.5-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libunistring/1.3-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libiconv/1.18-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libidn2/2.3.7-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/bzip2/1.0.8-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/ncurses/6.5-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/ScaLAPACK/2.2.2-gompi-2025a-fb/lib;/opt/software/pc2/EB-SW/software/FFTW.MPI/3.3.10-gompi-2025a/lib;/opt/software/pc2/EB-SW/software/FFTW/3.3.10-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/FlexiBLAS/3.4.5-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/AOCL-BLAS/5.1-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/OpenBLAS/0.3.29-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/OpenMPI/5.0.7-GCC-14.2.0/lib;/opt/software/pc2/EB-SW/software/UCC/1.3.0-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/PRRTE/3.0.8-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/PMIx/5.0.6-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libfabric/2.0.0-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/UCX/1.18.0-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libevent/2.1.12-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/hwloc/2.11.2-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libpciaccess/0.18.1-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/libxml2/2.13.4-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/XZ/5.6.3-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/numactl/2.0.19-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/binutils/2.42-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/zlib/1.3.1-GCCcore-14.2.0/lib;/opt/software/pc2/EB-SW/software/GCCcore/14.2.0/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
