#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

cmake_minimum_required(VERSION 3.22)
include(CMakeFindDependencyMacro)

if(NOT TARGET cp2k::cp2k)
  set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules"
                        ${CMAKE_MODULE_PATH})

  # Store CXX compiler id used in MKL package
  set(CP2K_CXX_COMPILER_ID GNU)
  if(NOT ${CMAKE_CXX_COMPILER_ID})
    set(CMAKE_CXX_COMPILER_ID ${CP2K_CXX_COMPILER_ID})
  endif()

  set(CP2K_BLAS_VENDOR OpenBLAS)
  set(CP2K_SCALAPACK_VENDOR auto)
  set(CP2K_BLAS_LINK_LIBRARIES /opt/software/pc2/EB-SW/software/OpenBLAS/0.3.29-GCC-14.2.0/lib/libopenblas.so)
  set(CP2K_LAPACK_LINK_LIBRARIES /opt/software/pc2/EB-SW/software/OpenBLAS/0.3.29-GCC-14.2.0/lib/libopenblas.so)
  set(CP2K_SCALAPACK_LINK_LIBRARIES /opt/software/pc2/EB-SW/software/ScaLAPACK/2.2.2-gompi-2025a-fb/lib/libscalapack.so)
  set(CP2K_CONFIG_PACKAGE ON)
  find_dependency(Lapack REQUIRED)

  # Define LAPACK and BLAS TARGETS

  if(ON)
    set(CP2K_USE_MPI ON)
    find_dependency(SCALAPACK REQUIRED)
  endif()
  unset(CP2K_CONFIG_PACKAGE)

  set(cp2k_VERSION 2026.1)

  find_dependency(DBCSR 2.8 REQUIRED)

  if(OFF)
    set(CP2K_USE_LIBXSMM OFF)
    find_dependency(LibXSMM REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_HIP OFF)
    # Find HIP
    find_dependency(hipfft REQUIRED IMPORTED CONFIG)
    find_dependency(hipblas REQUIRED IMPORTED CONFIG)
  endif()

  if(OFF)
    set(CP2K_USE_CUDA OFF)
    find_dependency(CUDAToolkit REQUIRED)
  endif()
  if(OFF)
    set(CP2K_USE_ELPA OFF)
    find_dependency(Elpa REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_DLAF OFF)
    find_dependency(DLAF REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_LIBXC OFF)
    find_dependency(LibXC 7 REQUIRED CONFIG)
  endif()

  if(OFF)
    set(CP2K_USE_COSMA OFF)
    find_dependency(cosma REQUIRED)
  endif()

  if(ON)
    set(CP2K_USE_MPI ON)
    find_dependency(MPI REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_FFTW3 OFF)
    find_dependency(Fftw REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_LIBINT2 OFF)
    find_dependency(Libint2 REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_SPGLIB OFF)
    find_dependency(LibSPG REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_SPLA OFF)
    find_dependency(SPLA REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_SIRIUS OFF)
    find_dependency(sirius REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_PLUMED OFF)
    find_dependency(Plumed REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_TORCH )
    find_dependency(Torch REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_HDF5 OFF)
    find_dependency(HDF5 REQUIRED COMPONENTS C Fortran)
  endif()

  if(OFF)
    set(CP2K_USE_DFTD4 OFF)
    find_dependency(dftd4 REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_DEEPMD OFF)
    find_dependency(DeepMD REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_PEXSI OFF)
    find_dependency(PEXSI REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_ACE OFF)
    find_dependency(ACE REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_ACE OFF)
    find_dependency(libsmeagol REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_VORI )
    find_dependency(LibVORI REQUIRED)
  endif()

  if(OFF)
    set(CP2K_USE_GREENX OFF)
    find_dependency(greenX REQUIRED)
  endif()

  include("${CMAKE_CURRENT_LIST_DIR}/cp2kTargets.cmake")

  # Clean-up module path
  list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")

endif()
