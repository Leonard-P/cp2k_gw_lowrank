!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routine for sorting an array
!> \note
!>      CP2K:
!>      Please use the interface defined in util.F for calling sort().
!>
!>      DBCSR:
!>      Please use the interface defined in dbcsr_toollib.F for calling sort().
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_array_sort

   USE kinds, ONLY: sp, dp, int_4, int_8

#include "../base/base_uses.f90"


   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cp_1d_s_sort, cp_1d_r_sort, cp_1d_i4_sort, cp_1d_i8_sort

CONTAINS


! **************************************************************************************************
!> \brief Sorts an array inplace using a combination of merge- and bubble-sort.
!>        It also returns the indices, which the elements had before the sort.
!> \param arr the array to sort
!> \param n length of array
!> \param indices returns elements-indices before the sort
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE cp_1d_s_sort(arr, n, indices)
         INTEGER, INTENT(IN)                      :: n
         REAL(kind=sp), DIMENSION(1:n), INTENT(INOUT)  :: arr
         integer, DIMENSION(1:n), INTENT(INOUT)   :: indices

         INTEGER                                  :: i
         REAL(kind=sp), ALLOCATABLE                    :: tmp_arr(:)
         INTEGER, ALLOCATABLE                     :: tmp_idx(:)

         IF (n > 1) THEN
            ! scratch space used during the merge step
            ALLOCATE (tmp_arr((n + 1)/2), tmp_idx((n + 1)/2))

            indices = [(i, i=1, n)]

            CALL cp_1d_s_sort_low(arr(1:n), indices, tmp_arr, tmp_idx)

            DEALLOCATE (tmp_arr, tmp_idx)
         ELSE IF (n > 0) THEN ! can be a frequent case in cp2k
            indices(1) = 1
         END IF

      END SUBROUTINE cp_1d_s_sort

! **************************************************************************************************
!> \brief The actual sort routing. Only cp_1d_s_sort and itself should call this.
!> \param arr the array to sort
!> \param indices elements-indices before the sort
!> \param tmp_arr scratch space
!> \param tmp_idx scratch space
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE cp_1d_s_sort_low(arr, indices, tmp_arr, tmp_idx)
         REAL(kind=sp), DIMENSION(:), INTENT(INOUT)                 :: arr
         INTEGER, DIMENSION(size(arr)), INTENT(INOUT)          :: indices
         REAL(kind=sp), DIMENSION((size(arr) + 1)/2), INTENT(INOUT)   :: tmp_arr
         INTEGER, DIMENSION((size(arr) + 1)/2), INTENT(INOUT)    :: tmp_idx
         REAL(kind=sp)                                              :: a
         INTEGER                                               :: t, m, i, j, k
         LOGICAL                                               :: swapped
         ! a,t:  used during swaping of elements in arr and indices

         swapped = .TRUE.

         ! If only a few elements are left we switch to bubble-sort for efficiency.
         IF (size(arr) <= 7) THEN ! 7 seems to be a good choice for the moment
            DO j = size(arr) - 1, 1, -1
               swapped = .FALSE.
               DO i = 1, j
                  IF (cp_1d_s_less_than(arr(i + 1), arr(i))) THEN
                     ! swap arr(i) with arr(i+1)
                     a = arr(i)
                     arr(i) = arr(i + 1)
                     arr(i + 1) = a
                     ! swap indices(i) with indices(i+1)
                     t = indices(i)
                     indices(i) = indices(i + 1)
                     indices(i + 1) = t
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT
            END DO
            RETURN
         END IF

         ! split list in half and recursively sort both sublists
         m = (size(arr) + 1)/2 ! index where we're going to divide the list in two
         CALL cp_1d_s_sort_low(arr(1:m), indices(1:m), tmp_arr, tmp_idx)
         CALL cp_1d_s_sort_low(arr(m + 1:), indices(m + 1:), tmp_arr, tmp_idx)

         ! Check for a special case: Can we just concate the two sorted sublists?
         ! This leads to O(n) scaling if the input is already sorted.
         IF (cp_1d_s_less_than(arr(m + 1), arr(m))) THEN
            ! ...no - let's merge the two sorted sublists arr(:m) and arr(m+1:)
            ! Merge will be performed directly in arr. Need backup of first sublist.
            tmp_arr(1:m) = arr(1:m)
            tmp_idx(1:m) = indices(1:m)
            i = 1; ! number of elemens consumed from 1st sublist
            j = 1; ! number of elemens consumed from 2nd sublist
            k = 1; ! number of elemens already merged

            DO WHILE (i <= m .and. j <= size(arr) - m)
               IF (cp_1d_s_less_than(arr(m + j), tmp_arr(i))) THEN
                  arr(k) = arr(m + j)
                  indices(k) = indices(m + j)
                  j = j + 1
               ELSE
                  arr(k) = tmp_arr(i)
                  indices(k) = tmp_idx(i)
                  i = i + 1
               END IF
               k = k + 1
            END DO

            ! One of the two sublist is now empty.
            ! Copy possibly remaining tail of 1st sublist
            DO WHILE (i <= m)
               arr(k) = tmp_arr(i)
               indices(k) = tmp_idx(i)
               i = i + 1
               k = k + 1
            END DO
            ! The possibly remaining tail of 2nd sublist is already at the right spot.
         END IF

      END SUBROUTINE cp_1d_s_sort_low


   PURE FUNCTION cp_1d_s_less_than(a, b) RESULT(res)
      REAL(kind=sp), INTENT(IN) :: a, b
      LOGICAL                   :: res
      res = a < b
   END FUNCTION cp_1d_s_less_than


! **************************************************************************************************
!> \brief Sorts an array inplace using a combination of merge- and bubble-sort.
!>        It also returns the indices, which the elements had before the sort.
!> \param arr the array to sort
!> \param n length of array
!> \param indices returns elements-indices before the sort
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE cp_1d_r_sort(arr, n, indices)
         INTEGER, INTENT(IN)                      :: n
         REAL(kind=dp), DIMENSION(1:n), INTENT(INOUT)  :: arr
         integer, DIMENSION(1:n), INTENT(INOUT)   :: indices

         INTEGER                                  :: i
         REAL(kind=dp), ALLOCATABLE                    :: tmp_arr(:)
         INTEGER, ALLOCATABLE                     :: tmp_idx(:)

         IF (n > 1) THEN
            ! scratch space used during the merge step
            ALLOCATE (tmp_arr((n + 1)/2), tmp_idx((n + 1)/2))

            indices = [(i, i=1, n)]

            CALL cp_1d_r_sort_low(arr(1:n), indices, tmp_arr, tmp_idx)

            DEALLOCATE (tmp_arr, tmp_idx)
         ELSE IF (n > 0) THEN ! can be a frequent case in cp2k
            indices(1) = 1
         END IF

      END SUBROUTINE cp_1d_r_sort

! **************************************************************************************************
!> \brief The actual sort routing. Only cp_1d_r_sort and itself should call this.
!> \param arr the array to sort
!> \param indices elements-indices before the sort
!> \param tmp_arr scratch space
!> \param tmp_idx scratch space
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE cp_1d_r_sort_low(arr, indices, tmp_arr, tmp_idx)
         REAL(kind=dp), DIMENSION(:), INTENT(INOUT)                 :: arr
         INTEGER, DIMENSION(size(arr)), INTENT(INOUT)          :: indices
         REAL(kind=dp), DIMENSION((size(arr) + 1)/2), INTENT(INOUT)   :: tmp_arr
         INTEGER, DIMENSION((size(arr) + 1)/2), INTENT(INOUT)    :: tmp_idx
         REAL(kind=dp)                                              :: a
         INTEGER                                               :: t, m, i, j, k
         LOGICAL                                               :: swapped
         ! a,t:  used during swaping of elements in arr and indices

         swapped = .TRUE.

         ! If only a few elements are left we switch to bubble-sort for efficiency.
         IF (size(arr) <= 7) THEN ! 7 seems to be a good choice for the moment
            DO j = size(arr) - 1, 1, -1
               swapped = .FALSE.
               DO i = 1, j
                  IF (cp_1d_r_less_than(arr(i + 1), arr(i))) THEN
                     ! swap arr(i) with arr(i+1)
                     a = arr(i)
                     arr(i) = arr(i + 1)
                     arr(i + 1) = a
                     ! swap indices(i) with indices(i+1)
                     t = indices(i)
                     indices(i) = indices(i + 1)
                     indices(i + 1) = t
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT
            END DO
            RETURN
         END IF

         ! split list in half and recursively sort both sublists
         m = (size(arr) + 1)/2 ! index where we're going to divide the list in two
         CALL cp_1d_r_sort_low(arr(1:m), indices(1:m), tmp_arr, tmp_idx)
         CALL cp_1d_r_sort_low(arr(m + 1:), indices(m + 1:), tmp_arr, tmp_idx)

         ! Check for a special case: Can we just concate the two sorted sublists?
         ! This leads to O(n) scaling if the input is already sorted.
         IF (cp_1d_r_less_than(arr(m + 1), arr(m))) THEN
            ! ...no - let's merge the two sorted sublists arr(:m) and arr(m+1:)
            ! Merge will be performed directly in arr. Need backup of first sublist.
            tmp_arr(1:m) = arr(1:m)
            tmp_idx(1:m) = indices(1:m)
            i = 1; ! number of elemens consumed from 1st sublist
            j = 1; ! number of elemens consumed from 2nd sublist
            k = 1; ! number of elemens already merged

            DO WHILE (i <= m .and. j <= size(arr) - m)
               IF (cp_1d_r_less_than(arr(m + j), tmp_arr(i))) THEN
                  arr(k) = arr(m + j)
                  indices(k) = indices(m + j)
                  j = j + 1
               ELSE
                  arr(k) = tmp_arr(i)
                  indices(k) = tmp_idx(i)
                  i = i + 1
               END IF
               k = k + 1
            END DO

            ! One of the two sublist is now empty.
            ! Copy possibly remaining tail of 1st sublist
            DO WHILE (i <= m)
               arr(k) = tmp_arr(i)
               indices(k) = tmp_idx(i)
               i = i + 1
               k = k + 1
            END DO
            ! The possibly remaining tail of 2nd sublist is already at the right spot.
         END IF

      END SUBROUTINE cp_1d_r_sort_low


   PURE FUNCTION cp_1d_r_less_than(a, b) RESULT(res)
      REAL(kind=dp), INTENT(IN) :: a, b
      LOGICAL                   :: res
      res = a < b
   END FUNCTION cp_1d_r_less_than


! **************************************************************************************************
!> \brief Sorts an array inplace using a combination of merge- and bubble-sort.
!>        It also returns the indices, which the elements had before the sort.
!> \param arr the array to sort
!> \param n length of array
!> \param indices returns elements-indices before the sort
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE cp_1d_i4_sort(arr, n, indices)
         INTEGER, INTENT(IN)                      :: n
         INTEGER(kind=int_4), DIMENSION(1:n), INTENT(INOUT)  :: arr
         integer, DIMENSION(1:n), INTENT(INOUT)   :: indices

         INTEGER                                  :: i
         INTEGER(kind=int_4), ALLOCATABLE                    :: tmp_arr(:)
         INTEGER, ALLOCATABLE                     :: tmp_idx(:)

         IF (n > 1) THEN
            ! scratch space used during the merge step
            ALLOCATE (tmp_arr((n + 1)/2), tmp_idx((n + 1)/2))

            indices = [(i, i=1, n)]

            CALL cp_1d_i4_sort_low(arr(1:n), indices, tmp_arr, tmp_idx)

            DEALLOCATE (tmp_arr, tmp_idx)
         ELSE IF (n > 0) THEN ! can be a frequent case in cp2k
            indices(1) = 1
         END IF

      END SUBROUTINE cp_1d_i4_sort

! **************************************************************************************************
!> \brief The actual sort routing. Only cp_1d_i4_sort and itself should call this.
!> \param arr the array to sort
!> \param indices elements-indices before the sort
!> \param tmp_arr scratch space
!> \param tmp_idx scratch space
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE cp_1d_i4_sort_low(arr, indices, tmp_arr, tmp_idx)
         INTEGER(kind=int_4), DIMENSION(:), INTENT(INOUT)                 :: arr
         INTEGER, DIMENSION(size(arr)), INTENT(INOUT)          :: indices
         INTEGER(kind=int_4), DIMENSION((size(arr) + 1)/2), INTENT(INOUT)   :: tmp_arr
         INTEGER, DIMENSION((size(arr) + 1)/2), INTENT(INOUT)    :: tmp_idx
         INTEGER(kind=int_4)                                              :: a
         INTEGER                                               :: t, m, i, j, k
         LOGICAL                                               :: swapped
         ! a,t:  used during swaping of elements in arr and indices

         swapped = .TRUE.

         ! If only a few elements are left we switch to bubble-sort for efficiency.
         IF (size(arr) <= 7) THEN ! 7 seems to be a good choice for the moment
            DO j = size(arr) - 1, 1, -1
               swapped = .FALSE.
               DO i = 1, j
                  IF (cp_1d_i4_less_than(arr(i + 1), arr(i))) THEN
                     ! swap arr(i) with arr(i+1)
                     a = arr(i)
                     arr(i) = arr(i + 1)
                     arr(i + 1) = a
                     ! swap indices(i) with indices(i+1)
                     t = indices(i)
                     indices(i) = indices(i + 1)
                     indices(i + 1) = t
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT
            END DO
            RETURN
         END IF

         ! split list in half and recursively sort both sublists
         m = (size(arr) + 1)/2 ! index where we're going to divide the list in two
         CALL cp_1d_i4_sort_low(arr(1:m), indices(1:m), tmp_arr, tmp_idx)
         CALL cp_1d_i4_sort_low(arr(m + 1:), indices(m + 1:), tmp_arr, tmp_idx)

         ! Check for a special case: Can we just concate the two sorted sublists?
         ! This leads to O(n) scaling if the input is already sorted.
         IF (cp_1d_i4_less_than(arr(m + 1), arr(m))) THEN
            ! ...no - let's merge the two sorted sublists arr(:m) and arr(m+1:)
            ! Merge will be performed directly in arr. Need backup of first sublist.
            tmp_arr(1:m) = arr(1:m)
            tmp_idx(1:m) = indices(1:m)
            i = 1; ! number of elemens consumed from 1st sublist
            j = 1; ! number of elemens consumed from 2nd sublist
            k = 1; ! number of elemens already merged

            DO WHILE (i <= m .and. j <= size(arr) - m)
               IF (cp_1d_i4_less_than(arr(m + j), tmp_arr(i))) THEN
                  arr(k) = arr(m + j)
                  indices(k) = indices(m + j)
                  j = j + 1
               ELSE
                  arr(k) = tmp_arr(i)
                  indices(k) = tmp_idx(i)
                  i = i + 1
               END IF
               k = k + 1
            END DO

            ! One of the two sublist is now empty.
            ! Copy possibly remaining tail of 1st sublist
            DO WHILE (i <= m)
               arr(k) = tmp_arr(i)
               indices(k) = tmp_idx(i)
               i = i + 1
               k = k + 1
            END DO
            ! The possibly remaining tail of 2nd sublist is already at the right spot.
         END IF

      END SUBROUTINE cp_1d_i4_sort_low


   PURE FUNCTION cp_1d_i4_less_than(a, b) RESULT(res)
      INTEGER(kind=int_4), INTENT(IN) :: a, b
      LOGICAL                         :: res
      res = a < b
   END FUNCTION cp_1d_i4_less_than


! **************************************************************************************************
!> \brief Sorts an array inplace using a combination of merge- and bubble-sort.
!>        It also returns the indices, which the elements had before the sort.
!> \param arr the array to sort
!> \param n length of array
!> \param indices returns elements-indices before the sort
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE cp_1d_i8_sort(arr, n, indices)
         INTEGER, INTENT(IN)                      :: n
         INTEGER(kind=int_8), DIMENSION(1:n), INTENT(INOUT)  :: arr
         integer, DIMENSION(1:n), INTENT(INOUT)   :: indices

         INTEGER                                  :: i
         INTEGER(kind=int_8), ALLOCATABLE                    :: tmp_arr(:)
         INTEGER, ALLOCATABLE                     :: tmp_idx(:)

         IF (n > 1) THEN
            ! scratch space used during the merge step
            ALLOCATE (tmp_arr((n + 1)/2), tmp_idx((n + 1)/2))

            indices = [(i, i=1, n)]

            CALL cp_1d_i8_sort_low(arr(1:n), indices, tmp_arr, tmp_idx)

            DEALLOCATE (tmp_arr, tmp_idx)
         ELSE IF (n > 0) THEN ! can be a frequent case in cp2k
            indices(1) = 1
         END IF

      END SUBROUTINE cp_1d_i8_sort

! **************************************************************************************************
!> \brief The actual sort routing. Only cp_1d_i8_sort and itself should call this.
!> \param arr the array to sort
!> \param indices elements-indices before the sort
!> \param tmp_arr scratch space
!> \param tmp_idx scratch space
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE cp_1d_i8_sort_low(arr, indices, tmp_arr, tmp_idx)
         INTEGER(kind=int_8), DIMENSION(:), INTENT(INOUT)                 :: arr
         INTEGER, DIMENSION(size(arr)), INTENT(INOUT)          :: indices
         INTEGER(kind=int_8), DIMENSION((size(arr) + 1)/2), INTENT(INOUT)   :: tmp_arr
         INTEGER, DIMENSION((size(arr) + 1)/2), INTENT(INOUT)    :: tmp_idx
         INTEGER(kind=int_8)                                              :: a
         INTEGER                                               :: t, m, i, j, k
         LOGICAL                                               :: swapped
         ! a,t:  used during swaping of elements in arr and indices

         swapped = .TRUE.

         ! If only a few elements are left we switch to bubble-sort for efficiency.
         IF (size(arr) <= 7) THEN ! 7 seems to be a good choice for the moment
            DO j = size(arr) - 1, 1, -1
               swapped = .FALSE.
               DO i = 1, j
                  IF (cp_1d_i8_less_than(arr(i + 1), arr(i))) THEN
                     ! swap arr(i) with arr(i+1)
                     a = arr(i)
                     arr(i) = arr(i + 1)
                     arr(i + 1) = a
                     ! swap indices(i) with indices(i+1)
                     t = indices(i)
                     indices(i) = indices(i + 1)
                     indices(i + 1) = t
                     swapped = .TRUE.
                  END IF
               END DO
               IF (.NOT. swapped) EXIT
            END DO
            RETURN
         END IF

         ! split list in half and recursively sort both sublists
         m = (size(arr) + 1)/2 ! index where we're going to divide the list in two
         CALL cp_1d_i8_sort_low(arr(1:m), indices(1:m), tmp_arr, tmp_idx)
         CALL cp_1d_i8_sort_low(arr(m + 1:), indices(m + 1:), tmp_arr, tmp_idx)

         ! Check for a special case: Can we just concate the two sorted sublists?
         ! This leads to O(n) scaling if the input is already sorted.
         IF (cp_1d_i8_less_than(arr(m + 1), arr(m))) THEN
            ! ...no - let's merge the two sorted sublists arr(:m) and arr(m+1:)
            ! Merge will be performed directly in arr. Need backup of first sublist.
            tmp_arr(1:m) = arr(1:m)
            tmp_idx(1:m) = indices(1:m)
            i = 1; ! number of elemens consumed from 1st sublist
            j = 1; ! number of elemens consumed from 2nd sublist
            k = 1; ! number of elemens already merged

            DO WHILE (i <= m .and. j <= size(arr) - m)
               IF (cp_1d_i8_less_than(arr(m + j), tmp_arr(i))) THEN
                  arr(k) = arr(m + j)
                  indices(k) = indices(m + j)
                  j = j + 1
               ELSE
                  arr(k) = tmp_arr(i)
                  indices(k) = tmp_idx(i)
                  i = i + 1
               END IF
               k = k + 1
            END DO

            ! One of the two sublist is now empty.
            ! Copy possibly remaining tail of 1st sublist
            DO WHILE (i <= m)
               arr(k) = tmp_arr(i)
               indices(k) = tmp_idx(i)
               i = i + 1
               k = k + 1
            END DO
            ! The possibly remaining tail of 2nd sublist is already at the right spot.
         END IF

      END SUBROUTINE cp_1d_i8_sort_low


   PURE FUNCTION cp_1d_i8_less_than(a, b) RESULT(res)
      INTEGER(kind=int_8), INTENT(IN) :: a, b
      LOGICAL                         :: res
      res = a < b
   END FUNCTION cp_1d_i8_less_than

END MODULE cp_array_sort
