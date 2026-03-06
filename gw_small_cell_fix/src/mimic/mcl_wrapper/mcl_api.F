!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper module for MiMiC Communication Library (MCL) routines
!> \par History
!>       05.2025 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mcl_api

   USE kinds, ONLY: int_4, int_8, sp, dp
#if defined(__MIMIC)
   USE mcl, ONLY: mcl_finalize, mcl_get_api_version, mcl_get_program_id, mcl_receive, mcl_send
#endif

#include "../../base/base_uses.f90"

   IMPLICIT NONE



   PRIVATE

   PUBLIC :: mcl_finalize, mcl_get_api_version, mcl_get_program_id, mcl_receive, mcl_send

#if !defined(__MIMIC)

      INTERFACE mcl_send
         MODULE PROCEDURE mcl_send_char
               MODULE PROCEDURE mcl_send_i4_0d
               MODULE PROCEDURE mcl_send_i4_1d
               MODULE PROCEDURE mcl_send_i4_2d
               MODULE PROCEDURE mcl_send_i4_3d
               MODULE PROCEDURE mcl_send_i8_0d
               MODULE PROCEDURE mcl_send_i8_1d
               MODULE PROCEDURE mcl_send_i8_2d
               MODULE PROCEDURE mcl_send_i8_3d
               MODULE PROCEDURE mcl_send_sp_0d
               MODULE PROCEDURE mcl_send_sp_1d
               MODULE PROCEDURE mcl_send_sp_2d
               MODULE PROCEDURE mcl_send_sp_3d
               MODULE PROCEDURE mcl_send_dp_0d
               MODULE PROCEDURE mcl_send_dp_1d
               MODULE PROCEDURE mcl_send_dp_2d
               MODULE PROCEDURE mcl_send_dp_3d
      END INTERFACE
      INTERFACE mcl_receive
         MODULE PROCEDURE mcl_receive_char
               MODULE PROCEDURE mcl_receive_i4_0d
               MODULE PROCEDURE mcl_receive_i4_1d
               MODULE PROCEDURE mcl_receive_i4_2d
               MODULE PROCEDURE mcl_receive_i4_3d
               MODULE PROCEDURE mcl_receive_i8_0d
               MODULE PROCEDURE mcl_receive_i8_1d
               MODULE PROCEDURE mcl_receive_i8_2d
               MODULE PROCEDURE mcl_receive_i8_3d
               MODULE PROCEDURE mcl_receive_sp_0d
               MODULE PROCEDURE mcl_receive_sp_1d
               MODULE PROCEDURE mcl_receive_sp_2d
               MODULE PROCEDURE mcl_receive_sp_3d
               MODULE PROCEDURE mcl_receive_dp_0d
               MODULE PROCEDURE mcl_receive_dp_1d
               MODULE PROCEDURE mcl_receive_dp_2d
               MODULE PROCEDURE mcl_receive_dp_3d
      END INTERFACE

CONTAINS

   SUBROUTINE mcl_finalize()
      CPABORT("This CP2K executable has not been linked against MiMiC.")
   END SUBROUTINE mcl_finalize

   SUBROUTINE mcl_get_api_version(version)
      INTEGER, DIMENSION(3) :: version, dummy_version
      CPABORT("This CP2K executable has not been linked against MiMiC.")
      dummy_version = version
   END SUBROUTINE mcl_get_api_version

   SUBROUTINE mcl_get_program_id(id)
      INTEGER :: id, dummy_id
      CPABORT("This CP2K executable has not been linked against MiMiC.")
      dummy_id = id
   END SUBROUTINE mcl_get_program_id

      SUBROUTINE mcl_send_char(buffer, length, tag, remote_id, err)
         CHARACTER(LEN=*), TARGET :: buffer
         INTEGER :: length, tag, remote_id, dummy_integer
         INTEGER, OPTIONAL, INTENT(OUT) :: err
         CHARACTER(LEN=:), ALLOCATABLE :: dummy_char
         CPABORT("This CP2K executable has not been linked against MiMiC.")
         dummy_integer = length
         dummy_integer = tag
         dummy_integer = remote_id
         err = dummy_integer
         dummy_char = TRIM(buffer)
      END SUBROUTINE mcl_send_char

            SUBROUTINE mcl_send_i4_0d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i4_0d
            SUBROUTINE mcl_send_i4_1d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i4_1d
            SUBROUTINE mcl_send_i4_2d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i4_2d
            SUBROUTINE mcl_send_i4_3d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i4_3d
            SUBROUTINE mcl_send_i8_0d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i8_0d
            SUBROUTINE mcl_send_i8_1d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i8_1d
            SUBROUTINE mcl_send_i8_2d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i8_2d
            SUBROUTINE mcl_send_i8_3d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_i8_3d
            SUBROUTINE mcl_send_sp_0d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_sp_0d
            SUBROUTINE mcl_send_sp_1d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_sp_1d
            SUBROUTINE mcl_send_sp_2d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_sp_2d
            SUBROUTINE mcl_send_sp_3d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_sp_3d
            SUBROUTINE mcl_send_dp_0d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_dp_0d
            SUBROUTINE mcl_send_dp_1d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_dp_1d
            SUBROUTINE mcl_send_dp_2d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_dp_2d
            SUBROUTINE mcl_send_dp_3d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_send_dp_3d
      SUBROUTINE mcl_receive_char(buffer, length, tag, remote_id, err)
         CHARACTER(LEN=*), TARGET :: buffer
         INTEGER :: length, tag, remote_id, dummy_integer
         INTEGER, OPTIONAL, INTENT(OUT) :: err
         CHARACTER(LEN=:), ALLOCATABLE :: dummy_char
         CPABORT("This CP2K executable has not been linked against MiMiC.")
         dummy_integer = length
         dummy_integer = tag
         dummy_integer = remote_id
         err = dummy_integer
         dummy_char = TRIM(buffer)
      END SUBROUTINE mcl_receive_char

            SUBROUTINE mcl_receive_i4_0d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i4_0d
            SUBROUTINE mcl_receive_i4_1d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i4_1d
            SUBROUTINE mcl_receive_i4_2d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i4_2d
            SUBROUTINE mcl_receive_i4_3d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_4), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_4), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i4_3d
            SUBROUTINE mcl_receive_i8_0d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i8_0d
            SUBROUTINE mcl_receive_i8_1d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i8_1d
            SUBROUTINE mcl_receive_i8_2d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i8_2d
            SUBROUTINE mcl_receive_i8_3d (buffer, length, tag, remote_id, err)
               INTEGER(KIND=int_8), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               INTEGER(KIND=int_8), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_i8_3d
            SUBROUTINE mcl_receive_sp_0d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_sp_0d
            SUBROUTINE mcl_receive_sp_1d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_sp_1d
            SUBROUTINE mcl_receive_sp_2d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_sp_2d
            SUBROUTINE mcl_receive_sp_3d (buffer, length, tag, remote_id, err)
               REAL(KIND=sp), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=sp), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_sp_3d
            SUBROUTINE mcl_receive_dp_0d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_dp_0d
            SUBROUTINE mcl_receive_dp_1d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_dp_1d
            SUBROUTINE mcl_receive_dp_2d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_dp_2d
            SUBROUTINE mcl_receive_dp_3d (buffer, length, tag, remote_id, err)
               REAL(KIND=dp), DIMENSION(:,:,:), TARGET :: buffer
               INTEGER :: length, tag, remote_id, dummy_integer
               INTEGER, OPTIONAL, INTENT(OUT) :: err
               REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: dummy_buffer
               CPABORT("This CP2K executable has not been linked against MiMiC.")
               dummy_integer = length
               dummy_integer = tag
               dummy_integer = remote_id
               err = dummy_integer
               dummy_buffer => buffer
            END SUBROUTINE mcl_receive_dp_3d

#endif

END MODULE mcl_api

