!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for allocating, copying and reshaping arrays.
!> \todo with fortran 2008 support, this should be replaced by plain ALLOCATE
!> \note in particular ALLOCATE(..., SOURCE=...) does not work in gcc 5.4.0, see also
!>       https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_allocate_wrap


   USE kinds, ONLY: dp

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: allocate_any

   INTERFACE allocate_any
         MODULE PROCEDURE allocate_1d
         MODULE PROCEDURE allocate_2d
         MODULE PROCEDURE allocate_3d
         MODULE PROCEDURE allocate_4d
         MODULE PROCEDURE allocate_5d
         MODULE PROCEDURE allocate_6d
         MODULE PROCEDURE allocate_7d
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_1d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(1), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(1), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(1)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1)))
               array(:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1)))
         END IF

      END SUBROUTINE allocate_1d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_2d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(2)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2)))
               array(:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2)))
         END IF

      END SUBROUTINE allocate_2d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_3d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(3)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3)))
               array(:,:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3)))
         END IF

      END SUBROUTINE allocate_3d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_4d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(4)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4)))
               array(:,:,:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4)))
         END IF

      END SUBROUTINE allocate_4d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_5d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:,:,:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(5)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5)))
               array(:,:,:,:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5)))
         END IF

      END SUBROUTINE allocate_5d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_6d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(6), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:,:,:,:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(6), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(6)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6)))
               array(:,:,:,:,:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6)), source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6)))
         END IF

      END SUBROUTINE allocate_6d
! **************************************************************************************************
!> \brief Allocate array according to shape_spec. Possibly assign array from source.
!> \note this does not fully replace Fortran RESHAPE intrinsic since source and target array must
!>       have same rank
!> \param array target array.
!> \param shape_spec shape of array to be allocated. If not specified, it is derived from source.
!> \param source source array to be copied to target array, must have same rank as target array.
!> \param order in which order to copy source to array (same convention as RESHAPE intrinsic).
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE allocate_7d(array, shape_spec, source, order)
         REAL(dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE, INTENT(OUT) :: array
         INTEGER, DIMENSION(7), INTENT(IN), OPTIONAL                    :: shape_spec
         REAL(dp), DIMENSION(:,:,:,:,:,:,:), INTENT(IN), OPTIONAL     :: source
         INTEGER, DIMENSION(7), INTENT(IN), OPTIONAL                    :: order
         INTEGER, DIMENSION(7)                                          :: shape_prv

         IF (PRESENT(shape_spec)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = shape_spec
            ELSE
               shape_prv = shape_spec
            END IF
         ELSEIF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               shape_prv(order) = SHAPE(source)
            ELSE
               shape_prv = SHAPE(source)
            END IF
         ELSE
            CPABORT("either source or shape_spec must be present")
         END IF

         IF (PRESENT(source)) THEN
            IF (PRESENT(order)) THEN
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6), shape_prv(7)))
               array(:,:,:,:,:,:,:) = RESHAPE(source, shape_prv, order=order)
            ELSE
               ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6), shape_prv(7)),&
                   & source=source)
            END IF
         ELSE
            ALLOCATE (array(shape_prv(1), shape_prv(2), shape_prv(3), shape_prv(4), shape_prv(5), shape_prv(6), shape_prv(7)))
         END IF

      END SUBROUTINE allocate_7d
END MODULE dbt_allocate_wrap
