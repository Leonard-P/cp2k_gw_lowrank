!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for the memory handling.
!> \par History
!>      (12.2017) remove stop_memory
!> \author Matthias Krack (25.06.1999)
! **************************************************************************************************
MODULE memory_utilities

   USE kinds, ONLY: dp, int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'memory_utilities'

   PUBLIC :: reallocate

   INTERFACE reallocate
      MODULE PROCEDURE reallocate_c1, reallocate_c2, reallocate_c3, reallocate_c4, &
         reallocate_i1, reallocate_i2, reallocate_i3, reallocate_i4, &
         reallocate_r1, reallocate_r2, reallocate_r3, reallocate_r4, &
         reallocate_r5, reallocate_s1, reallocate_l1, reallocate_8i1, &
         reallocate_8i2
   END INTERFACE

CONTAINS


! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type COMPLEX(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_c1 (p, lb1_new,ub1_new)
         COMPLEX(KIND=dp), &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

            COMPLEX(KIND=dp), &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = (0.0_dp, 0.0_dp)

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_c1
! **************************************************************************************************
!> \brief (Re)Allocate a 2D vector of type COMPLEX(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
! **************************************************************************************************
      SUBROUTINE reallocate_c2 (p, lb1_new,ub1_new,lb2_new,ub2_new)
         COMPLEX(KIND=dp), &
            DIMENSION(:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old

            COMPLEX(KIND=dp), &
            DIMENSION(:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new))
         p = (0.0_dp, 0.0_dp)

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_c2
! **************************************************************************************************
!> \brief (Re)Allocate a 3D vector of type COMPLEX(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
! **************************************************************************************************
      SUBROUTINE reallocate_c3 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
         COMPLEX(KIND=dp), &
            DIMENSION(:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old

            COMPLEX(KIND=dp), &
            DIMENSION(:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new))
         p = (0.0_dp, 0.0_dp)

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_c3
! **************************************************************************************************
!> \brief (Re)Allocate a 4D vector of type COMPLEX(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
!> \param lb4_new new lower bound for dimension 4
!> \param ub4_new new upper bound for dimension 4
! **************************************************************************************************
      SUBROUTINE reallocate_c4 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
         COMPLEX(KIND=dp), &
            DIMENSION(:,:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old
            INTEGER :: lb4, lb4_old, ub4, ub4_old

            COMPLEX(KIND=dp), &
            DIMENSION(:,:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
               lb4_old = LBOUND(p, 4)
               ub4_old = UBOUND(p, 4)
               lb4 = MAX(lb4_new, lb4_old)
               ub4 = MIN(ub4_new, ub4_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new))
         p = (0.0_dp, 0.0_dp)

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_c4
! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type INTEGER with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_i1 (p, lb1_new,ub1_new)
         INTEGER, &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

            INTEGER, &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_i1
! **************************************************************************************************
!> \brief (Re)Allocate a 2D vector of type INTEGER with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
! **************************************************************************************************
      SUBROUTINE reallocate_i2 (p, lb1_new,ub1_new,lb2_new,ub2_new)
         INTEGER, &
            DIMENSION(:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old

            INTEGER, &
            DIMENSION(:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_i2
! **************************************************************************************************
!> \brief (Re)Allocate a 3D vector of type INTEGER with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
! **************************************************************************************************
      SUBROUTINE reallocate_i3 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
         INTEGER, &
            DIMENSION(:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old

            INTEGER, &
            DIMENSION(:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_i3
! **************************************************************************************************
!> \brief (Re)Allocate a 4D vector of type INTEGER with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
!> \param lb4_new new lower bound for dimension 4
!> \param ub4_new new upper bound for dimension 4
! **************************************************************************************************
      SUBROUTINE reallocate_i4 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
         INTEGER, &
            DIMENSION(:,:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old
            INTEGER :: lb4, lb4_old, ub4, ub4_old

            INTEGER, &
            DIMENSION(:,:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
               lb4_old = LBOUND(p, 4)
               ub4_old = UBOUND(p, 4)
               lb4 = MAX(lb4_new, lb4_old)
               ub4 = MIN(ub4_new, ub4_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_i4
! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type INTEGER(KIND=int_8) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_8i1 (p, lb1_new,ub1_new)
         INTEGER(KIND=int_8), &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

            INTEGER(KIND=int_8), &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_8i1
! **************************************************************************************************
!> \brief (Re)Allocate a 2D vector of type INTEGER(KIND=int_8) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
! **************************************************************************************************
      SUBROUTINE reallocate_8i2 (p, lb1_new,ub1_new,lb2_new,ub2_new)
         INTEGER(KIND=int_8), &
            DIMENSION(:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old

            INTEGER(KIND=int_8), &
            DIMENSION(:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new))
         p = 0

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_8i2
! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type REAL(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_r1 (p, lb1_new,ub1_new)
         REAL(KIND=dp), &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

            REAL(KIND=dp), &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = 0.0_dp

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_r1
! **************************************************************************************************
!> \brief (Re)Allocate a 2D vector of type REAL(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
! **************************************************************************************************
      SUBROUTINE reallocate_r2 (p, lb1_new,ub1_new,lb2_new,ub2_new)
         REAL(KIND=dp), &
            DIMENSION(:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old

            REAL(KIND=dp), &
            DIMENSION(:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new))
         p = 0.0_dp

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2) = work(lb1:ub1,lb2:ub2)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_r2
! **************************************************************************************************
!> \brief (Re)Allocate a 3D vector of type REAL(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
! **************************************************************************************************
      SUBROUTINE reallocate_r3 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new)
         REAL(KIND=dp), &
            DIMENSION(:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old

            REAL(KIND=dp), &
            DIMENSION(:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new))
         p = 0.0_dp

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3) = work(lb1:ub1,lb2:ub2,lb3:ub3)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_r3
! **************************************************************************************************
!> \brief (Re)Allocate a 4D vector of type REAL(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
!> \param lb4_new new lower bound for dimension 4
!> \param ub4_new new upper bound for dimension 4
! **************************************************************************************************
      SUBROUTINE reallocate_r4 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new)
         REAL(KIND=dp), &
            DIMENSION(:,:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old
            INTEGER :: lb4, lb4_old, ub4, ub4_old

            REAL(KIND=dp), &
            DIMENSION(:,:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
               lb4_old = LBOUND(p, 4)
               ub4_old = UBOUND(p, 4)
               lb4 = MAX(lb4_new, lb4_old)
               ub4 = MIN(ub4_new, ub4_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new))
         p = 0.0_dp

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_r4
! **************************************************************************************************
!> \brief (Re)Allocate a 5D vector of type REAL(KIND=dp) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
!> \param lb2_new new lower bound for dimension 2
!> \param ub2_new new upper bound for dimension 2
!> \param lb3_new new lower bound for dimension 3
!> \param ub3_new new upper bound for dimension 3
!> \param lb4_new new lower bound for dimension 4
!> \param ub4_new new upper bound for dimension 4
!> \param lb5_new new lower bound for dimension 5
!> \param ub5_new new upper bound for dimension 5
! **************************************************************************************************
      SUBROUTINE reallocate_r5 (p, lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new)
         REAL(KIND=dp), &
            DIMENSION(:,:,:,:,:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new,lb2_new,ub2_new,lb3_new,ub3_new,lb4_new,ub4_new,lb5_new,ub5_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old
            INTEGER :: lb2, lb2_old, ub2, ub2_old
            INTEGER :: lb3, lb3_old, ub3, ub3_old
            INTEGER :: lb4, lb4_old, ub4, ub4_old
            INTEGER :: lb5, lb5_old, ub5, ub5_old

            REAL(KIND=dp), &
            DIMENSION(:,:,:,:,:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
               lb2_old = LBOUND(p, 2)
               ub2_old = UBOUND(p, 2)
               lb2 = MAX(lb2_new, lb2_old)
               ub2 = MIN(ub2_new, ub2_old)
               lb3_old = LBOUND(p, 3)
               ub3_old = UBOUND(p, 3)
               lb3 = MAX(lb3_new, lb3_old)
               ub3 = MIN(ub3_new, ub3_old)
               lb4_old = LBOUND(p, 4)
               ub4_old = UBOUND(p, 4)
               lb4 = MAX(lb4_new, lb4_old)
               ub4 = MIN(ub4_new, ub4_old)
               lb5_old = LBOUND(p, 5)
               ub5_old = UBOUND(p, 5)
               lb5 = MAX(lb5_new, lb5_old)
               ub5 = MIN(ub5_new, ub5_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new,lb2_new:ub2_new,lb3_new:ub3_new,lb4_new:ub4_new,lb5_new:ub5_new))
         p = 0.0_dp

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5) = work(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4,lb5:ub5)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_r5
! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type LOGICAL with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_l1 (p, lb1_new,ub1_new)
         LOGICAL, &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

            LOGICAL, &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = .FALSE.

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_l1
! **************************************************************************************************
!> \brief (Re)Allocate a 1D vector of type CHARACTER(LEN=*) with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
!> \param lb1_new new lower bound for dimension 1
!> \param ub1_new new upper bound for dimension 1
! **************************************************************************************************
      SUBROUTINE reallocate_s1 (p, lb1_new,ub1_new)
         CHARACTER(LEN=*), &
            DIMENSION(:), &
            POINTER, INTENT(INOUT) :: p

         INTEGER, INTENT(IN) :: &
            lb1_new,ub1_new

            INTEGER :: lb1, lb1_old, ub1, ub1_old

         CHARACTER(LEN=LEN(p)), &
            DIMENSION(:), &
            POINTER :: work

         NULLIFY (work)

         IF (ASSOCIATED(p)) THEN
               lb1_old = LBOUND(p, 1)
               ub1_old = UBOUND(p, 1)
               lb1 = MAX(lb1_new, lb1_old)
               ub1 = MIN(ub1_new, ub1_old)
            work => p
         END IF

         ALLOCATE (p(lb1_new:ub1_new))
         p = ""

         IF (ASSOCIATED(work)) THEN
            p(lb1:ub1) = work(lb1:ub1)
            DEALLOCATE (work)
         END IF

      END SUBROUTINE reallocate_s1

END MODULE memory_utilities
