!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Representation of arbitrary number of 1d integer arrays with arbitrary sizes.
!>        This is needed for generic handling of dimension-specific tensor quantities
!>        (such as block index).
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_array_list_methods



   USE dbt_index, ONLY: dbt_inverse_order
   USE dbt_allocate_wrap, ONLY: allocate_any

#include "../base/base_uses.f90"
#if defined(__LIBXSMM)
#  include "libxsmm_version.h"
#endif

#if CPVERSION_CHECK(1, 11, <=, LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
   USE libxsmm, ONLY: libxsmm_diff
#  define PURE_ARRAY_EQ
#else
#  define PURE_ARRAY_EQ PURE
#endif

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_array_list_methods'

   PUBLIC  :: &
      array_eq_i, &
      array_list, &
      array_offsets, &
      array_sublist, &
      create_array_list, &
      destroy_array_list, &
      get_array_elements, &
      get_arrays, &
      get_ith_array, &
      number_of_arrays, &
      reorder_arrays, &
      sizes_of_arrays, &
      sum_of_arrays, &
      check_equal

   TYPE array_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: col_data
      INTEGER, DIMENSION(:), ALLOCATABLE :: ptr
   END TYPE

   INTERFACE get_ith_array
      MODULE PROCEDURE allocate_and_get_ith_array
      MODULE PROCEDURE get_ith_array
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief number of arrays stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION number_of_arrays(list)
      TYPE(array_list), INTENT(IN) :: list
      INTEGER                      :: number_of_arrays

      number_of_arrays = SIZE(list%ptr) - 1

   END FUNCTION number_of_arrays

! **************************************************************************************************
!> \brief Get an element for each array.
!> \param indices element index for each array
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_array_elements(list, indices)
      TYPE(array_list), INTENT(IN)                           :: list
      INTEGER, DIMENSION(number_of_arrays(list)), INTENT(IN) :: indices
      INTEGER, DIMENSION(number_of_arrays(list))             :: get_array_elements

      INTEGER                                                :: i, ind

      DO i = 1, SIZE(indices)
         ind = indices(i) + list%ptr(i) - 1
         get_array_elements(i) = list%col_data(ind)
      END DO

   END FUNCTION get_array_elements

! **************************************************************************************************
!> \brief collects any number of arrays of different sizes into a single array (list%col_data),
!>        storing the indices that start a new array (list%ptr).
!> \param list list of arrays
!> \param ndata number of arrays
!> \param data arrays 1 and 2
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE create_array_list(list, ndata, data_1, data_2, data_3, data_4)
      TYPE(array_list), INTENT(OUT)               :: list
      INTEGER, INTENT(IN)                         :: ndata
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: data_1, data_2, data_3, data_4
      INTEGER                                     :: ptr, size_all

      size_all = 0

         IF (ndata >= 1) THEN
            CPASSERT(PRESENT(data_1))
            size_all = size_all + SIZE(data_1)
         END IF
         IF (ndata >= 2) THEN
            CPASSERT(PRESENT(data_2))
            size_all = size_all + SIZE(data_2)
         END IF
         IF (ndata >= 3) THEN
            CPASSERT(PRESENT(data_3))
            size_all = size_all + SIZE(data_3)
         END IF
         IF (ndata >= 4) THEN
            CPASSERT(PRESENT(data_4))
            size_all = size_all + SIZE(data_4)
         END IF

      ALLOCATE (list%ptr(ndata + 1))
      ALLOCATE (list%col_data(size_all))

      ptr = 1
      list%ptr(1) = ptr

         IF (ndata >= 1) THEN
            list%col_data(ptr:ptr + SIZE(data_1) - 1) = data_1 (:)
            ptr = ptr + SIZE(data_1)
            list%ptr(2) = ptr
         END IF
         IF (ndata >= 2) THEN
            list%col_data(ptr:ptr + SIZE(data_2) - 1) = data_2 (:)
            ptr = ptr + SIZE(data_2)
            list%ptr(3) = ptr
         END IF
         IF (ndata >= 3) THEN
            list%col_data(ptr:ptr + SIZE(data_3) - 1) = data_3 (:)
            ptr = ptr + SIZE(data_3)
            list%ptr(4) = ptr
         END IF
         IF (ndata >= 4) THEN
            list%col_data(ptr:ptr + SIZE(data_4) - 1) = data_4 (:)
            ptr = ptr + SIZE(data_4)
            list%ptr(5) = ptr
         END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief extract a subset of arrays
!> \param list list of arrays
!> \param i_selected array numbers to retrieve
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION array_sublist(list, i_selected)
      TYPE(array_list), INTENT(IN)                           :: list
      INTEGER, DIMENSION(:), INTENT(IN)                      :: i_selected
      TYPE(array_list)                                       :: array_sublist
      INTEGER :: ndata
      INTEGER, ALLOCATABLE, DIMENSION(:) :: data_1, data_2, data_3, data_4

      ndata = SIZE(i_selected)

         IF (ndata == 1) THEN
            CALL get_arrays(list, data_1, i_selected=i_selected)
            CALL create_array_list(array_sublist, ndata, data_1)
         END IF
         IF (ndata == 2) THEN
            CALL get_arrays(list, data_1, data_2, i_selected=i_selected)
            CALL create_array_list(array_sublist, ndata, data_1, data_2)
         END IF
         IF (ndata == 3) THEN
            CALL get_arrays(list, data_1, data_2, data_3, i_selected=i_selected)
            CALL create_array_list(array_sublist, ndata, data_1, data_2, data_3)
         END IF
         IF (ndata == 4) THEN
            CALL get_arrays(list, data_1, data_2, data_3, data_4, i_selected=i_selected)
            CALL create_array_list(array_sublist, ndata, data_1, data_2, data_3, data_4)
         END IF
   END FUNCTION

! **************************************************************************************************
!> \brief destroy array list.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE destroy_array_list(list)
      TYPE(array_list), INTENT(INOUT) :: list

      DEALLOCATE (list%ptr, list%col_data)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Get all arrays contained in list
!> \param data arrays 1 and 2
!> \param i_selected array numbers to retrieve (if not present, all arrays are returned)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE get_arrays(list, data_1, data_2, data_3, data_4, i_selected)
      !! Get all arrays contained in list
      TYPE(array_list), INTENT(IN)                       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: data_1, data_2, data_3, data_4
      INTEGER, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: i_selected
      INTEGER                                            :: i, ndata
      INTEGER, DIMENSION(number_of_arrays(list))         :: o

      o(:) = 0
      IF (PRESENT(i_selected)) THEN
         ndata = SIZE(i_selected)
         o(1:ndata) = i_selected(:)
      ELSE
         ndata = number_of_arrays(list)
         o(1:ndata) = (/(i, i=1, ndata)/)
      END IF

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
            IF (ndata > 0) THEN
               ALLOCATE (data_1, source=col_data(ptr(o(1)):ptr(o(1) + 1) - 1))
            END IF
            IF (ndata > 1) THEN
               ALLOCATE (data_2, source=col_data(ptr(o(2)):ptr(o(2) + 1) - 1))
            END IF
            IF (ndata > 2) THEN
               ALLOCATE (data_3, source=col_data(ptr(o(3)):ptr(o(3) + 1) - 1))
            END IF
            IF (ndata > 3) THEN
               ALLOCATE (data_4, source=col_data(ptr(o(4)):ptr(o(4) + 1) - 1))
            END IF
      END ASSOCIATE

   END SUBROUTINE get_arrays

! **************************************************************************************************
!> \brief get ith array
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE get_ith_array(list, i, array_size, array)
      TYPE(array_list), INTENT(IN)                    :: list
      INTEGER, INTENT(IN)                             :: i
      INTEGER, INTENT(IN)                             :: array_size
      INTEGER, DIMENSION(array_size), INTENT(OUT)     :: array

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
         CPASSERT(i <= number_of_arrays(list))

         array(:) = col_data(ptr(i):ptr(i + 1) - 1)

      END ASSOCIATE

   END SUBROUTINE

! **************************************************************************************************
!> \brief get ith array
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE allocate_and_get_ith_array(list, i, array)
      TYPE(array_list), INTENT(IN)                    :: list
      INTEGER, INTENT(IN)                             :: i
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: array

      ASSOCIATE (ptr => list%ptr, col_data => list%col_data)
         CPASSERT(i <= number_of_arrays(list))

         ALLOCATE (array, source=col_data(ptr(i):ptr(i + 1) - 1))
      END ASSOCIATE
   END SUBROUTINE

! **************************************************************************************************
!> \brief sizes of arrays stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION sizes_of_arrays(list)
      TYPE(array_list), INTENT(IN)       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sizes_of_arrays

      INTEGER                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sizes_of_arrays(num_data))
      DO i_data = 1, num_data
         sizes_of_arrays(i_data) = list%ptr(i_data + 1) - list%ptr(i_data)
      END DO
   END FUNCTION sizes_of_arrays

! **************************************************************************************************
!> \brief sum of all elements for each array stored in list
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION sum_of_arrays(list)
      TYPE(array_list), INTENT(IN)       :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sum_of_arrays

      INTEGER                            :: i_data, num_data

      num_data = number_of_arrays(list)
      ALLOCATE (sum_of_arrays(num_data))
      DO i_data = 1, num_data
         sum_of_arrays(i_data) = SUM(list%col_data(list%ptr(i_data):list%ptr(i_data + 1) - 1))
      END DO

   END FUNCTION sum_of_arrays

! **************************************************************************************************
!> \brief partial sums of array elements.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE array_offsets(list_in, list_out)
      TYPE(array_list), INTENT(IN)  :: list_in
      TYPE(array_list), INTENT(OUT) :: list_out

      INTEGER                       :: i_data, i_ptr, num_data, partial_sum

      num_data = number_of_arrays(list_in)
      ALLOCATE (list_out%ptr, source=list_in%ptr)
      ALLOCATE (list_out%col_data(SIZE(list_in%col_data)))
      DO i_data = 1, num_data
         partial_sum = 1
         DO i_ptr = list_out%ptr(i_data), list_out%ptr(i_data + 1) - 1
            list_out%col_data(i_ptr) = partial_sum
            partial_sum = partial_sum + list_in%col_data(i_ptr)
         END DO
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief reorder array list.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE reorder_arrays(list_in, list_out, order)
      TYPE(array_list), INTENT(IN)                     :: list_in
      TYPE(array_list), INTENT(OUT)                    :: list_out
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: data_1, data_2, data_3, data_4
      INTEGER, DIMENSION(number_of_arrays(list_in)), &
         INTENT(IN)                                    :: order

         IF (number_of_arrays(list_in) == 2) THEN
            CALL get_arrays(list_in, data_1, data_2, i_selected=dbt_inverse_order(order))
            CALL create_array_list(list_out, number_of_arrays(list_in), &
                                   data_1, data_2)
         END IF
         IF (number_of_arrays(list_in) == 3) THEN
            CALL get_arrays(list_in, data_1, data_2, data_3, i_selected=dbt_inverse_order(order))
            CALL create_array_list(list_out, number_of_arrays(list_in), &
                                   data_1, data_2, data_3)
         END IF
         IF (number_of_arrays(list_in) == 4) THEN
            CALL get_arrays(list_in, data_1, data_2, data_3, data_4, i_selected=dbt_inverse_order(order))
            CALL create_array_list(list_out, number_of_arrays(list_in), &
                                   data_1, data_2, data_3, data_4)
         END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief check whether two array lists are equal
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION check_equal(list1, list2)
      TYPE(array_list), INTENT(IN)  :: list1, list2
      LOGICAL :: check_equal

      check_equal = array_eq_i(list1%col_data, list2%col_data) .AND. array_eq_i(list1%ptr, list2%ptr)
   END FUNCTION

! **************************************************************************************************
!> \brief check whether two arrays are equal
!> \author Patrick Seewald
! **************************************************************************************************
   PURE_ARRAY_EQ FUNCTION array_eq_i(arr1, arr2)
      INTEGER, INTENT(IN), DIMENSION(:) :: arr1
      INTEGER, INTENT(IN), DIMENSION(:) :: arr2
      LOGICAL                           :: array_eq_i

#if CPVERSION_CHECK(1, 11, <=, LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
      array_eq_i = .NOT. libxsmm_diff(arr1, arr2)
#else
      array_eq_i = .FALSE.
      IF (SIZE(arr1) == SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)
#endif
   END FUNCTION

END MODULE dbt_array_list_methods
