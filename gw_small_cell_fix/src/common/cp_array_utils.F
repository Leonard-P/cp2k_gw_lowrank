!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!


! **************************************************************************************************
!> \brief various utilities that regard array of different kinds:
!>      output, allocation,...
!>      maybe it is not a good idea mixing output and memeory utils...
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_array_utils
   USE machine, ONLY: m_flush
   USE cp_log_handling, ONLY: cp_to_string

   USE kinds, ONLY: dp, int_4

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN = 'cp_array_utils'

      PUBLIC :: cp_1d_r_p_type, &
                cp_2d_r_p_type, &
                cp_3d_r_p_type, &
                cp_1d_r_cp_type, &
                cp_2d_r_cp_type, &
                cp_3d_r_cp_type, &
                cp_1d_r_guarantee_size, &
                cp_1d_r_write, &
                cp_2d_r_write, &
                cp_2d_r_guarantee_size, &
                cp_1d_r_bsearch
      PUBLIC :: cp_1d_c_p_type, &
                cp_2d_c_p_type, &
                cp_3d_c_p_type, &
                cp_1d_c_cp_type, &
                cp_2d_c_cp_type, &
                cp_3d_c_cp_type, &
                cp_1d_c_guarantee_size, &
                cp_1d_c_write, &
                cp_2d_c_write, &
                cp_2d_c_guarantee_size, &
                cp_1d_c_bsearch
      PUBLIC :: cp_1d_i_p_type, &
                cp_2d_i_p_type, &
                cp_3d_i_p_type, &
                cp_1d_i_cp_type, &
                cp_2d_i_cp_type, &
                cp_3d_i_cp_type, &
                cp_1d_i_guarantee_size, &
                cp_1d_i_write, &
                cp_2d_i_write, &
                cp_2d_i_guarantee_size, &
                cp_1d_i_bsearch
      PUBLIC :: cp_1d_logical_p_type, &
                cp_2d_logical_p_type, &
                cp_3d_logical_p_type, &
                cp_1d_logical_cp_type, &
                cp_2d_logical_cp_type, &
                cp_3d_logical_cp_type, &
                cp_1d_logical_guarantee_size, &
                cp_1d_logical_write, &
                cp_2d_logical_write, &
                cp_2d_logical_guarantee_size, &
                cp_1d_logical_bsearch

   ! generic interfaces
   PUBLIC :: cp_guarantee_size

   INTERFACE cp_guarantee_size
         MODULE PROCEDURE cp_1d_r_guarantee_size, &
            cp_2d_r_guarantee_size
         MODULE PROCEDURE cp_1d_c_guarantee_size, &
            cp_2d_c_guarantee_size
         MODULE PROCEDURE cp_1d_i_guarantee_size, &
            cp_2d_i_guarantee_size
         MODULE PROCEDURE cp_1d_logical_guarantee_size, &
            cp_2d_logical_guarantee_size
   END INTERFACE

!***


! **************************************************************************************************
!> \brief represent a pointer to a 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_r_p_type
         REAL(kind=dp), dimension(:), pointer :: array => NULL()
      end type cp_1d_r_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_r_p_type
         REAL(kind=dp), dimension(:, :), pointer :: array => NULL()
      end type cp_2d_r_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_r_p_type
         REAL(kind=dp), dimension(:, :, :), pointer :: array => NULL()
      end type cp_3d_r_p_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_r_cp_type
         REAL(kind=dp), dimension(:), contiguous, pointer :: array => NULL()
      end type cp_1d_r_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_r_cp_type
         REAL(kind=dp), dimension(:, :), contiguous, pointer :: array => NULL()
      end type cp_2d_r_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_r_cp_type
         REAL(kind=dp), dimension(:, :, :), contiguous, pointer :: array => NULL()
      end type cp_3d_r_cp_type


! **************************************************************************************************
!> \brief represent a pointer to a 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_c_p_type
         COMPLEX(KIND=dp), dimension(:), pointer :: array => NULL()
      end type cp_1d_c_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_c_p_type
         COMPLEX(KIND=dp), dimension(:, :), pointer :: array => NULL()
      end type cp_2d_c_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_c_p_type
         COMPLEX(KIND=dp), dimension(:, :, :), pointer :: array => NULL()
      end type cp_3d_c_p_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_c_cp_type
         COMPLEX(KIND=dp), dimension(:), contiguous, pointer :: array => NULL()
      end type cp_1d_c_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_c_cp_type
         COMPLEX(KIND=dp), dimension(:, :), contiguous, pointer :: array => NULL()
      end type cp_2d_c_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_c_cp_type
         COMPLEX(KIND=dp), dimension(:, :, :), contiguous, pointer :: array => NULL()
      end type cp_3d_c_cp_type


! **************************************************************************************************
!> \brief represent a pointer to a 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_i_p_type
         INTEGER(kind=int_4), dimension(:), pointer :: array => NULL()
      end type cp_1d_i_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_i_p_type
         INTEGER(kind=int_4), dimension(:, :), pointer :: array => NULL()
      end type cp_2d_i_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_i_p_type
         INTEGER(kind=int_4), dimension(:, :, :), pointer :: array => NULL()
      end type cp_3d_i_p_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_i_cp_type
         INTEGER(kind=int_4), dimension(:), contiguous, pointer :: array => NULL()
      end type cp_1d_i_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_i_cp_type
         INTEGER(kind=int_4), dimension(:, :), contiguous, pointer :: array => NULL()
      end type cp_2d_i_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_i_cp_type
         INTEGER(kind=int_4), dimension(:, :, :), contiguous, pointer :: array => NULL()
      end type cp_3d_i_cp_type


! **************************************************************************************************
!> \brief represent a pointer to a 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_logical_p_type
         logical, dimension(:), pointer :: array => NULL()
      end type cp_1d_logical_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_logical_p_type
         logical, dimension(:, :), pointer :: array => NULL()
      end type cp_2d_logical_p_type

! **************************************************************************************************
!> \brief represent a pointer to a 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_logical_p_type
         logical, dimension(:, :, :), pointer :: array => NULL()
      end type cp_3d_logical_p_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 1d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_1d_logical_cp_type
         logical, dimension(:), contiguous, pointer :: array => NULL()
      end type cp_1d_logical_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 2d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_2d_logical_cp_type
         logical, dimension(:, :), contiguous, pointer :: array => NULL()
      end type cp_2d_logical_cp_type

! **************************************************************************************************
!> \brief represent a pointer to a contiguous 3d array
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      type cp_3d_logical_cp_type
         logical, dimension(:, :, :), contiguous, pointer :: array => NULL()
      end type cp_3d_logical_cp_type


CONTAINS

! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated paretized list
! **************************************************************************************************
      SUBROUTINE cp_1d_r_write(array, unit_nr, el_format)
         REAL(kind=dp), INTENT(in) :: array(:)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(es14.6)"

         WRITE (unit=unit_nr, fmt="('( ')", advance="no", iostat=iostat)
         CPASSERT(iostat == 0)
         IF (PRESENT(el_format)) THEN
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(SIZE(array))
         ELSE
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(SIZE(array))
         END IF
         WRITE (unit=unit_nr, fmt="(' )')")
         call m_flush(unit_nr)

      END SUBROUTINE cp_1d_r_write

! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated parentized list
! **************************************************************************************************
      SUBROUTINE cp_2d_r_write(array, unit_nr, el_format)
         REAL(kind=dp), INTENT(in) :: array(:, :)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(es14.6)"
         CHARACTER(len=200) :: fmtstr
         CHARACTER(len=25) :: nRiga

         nRiga = cp_to_string(SIZE(array, 2))
         DO i = 1, SIZE(array, 1)
            IF (PRESENT(el_format)) THEN
               fmtstr = '(" ",'//nRiga//el_format//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            ELSE
               fmtstr = '(" ",'//nRiga//defaultFormat//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            END IF
            CPASSERT(iostat == 0)
         END DO
         call m_flush(unit_nr)
      END SUBROUTINE cp_2d_r_write

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n the wanted size
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_1d_r_guarantee_size(array, n)
         REAL(kind=dp), POINTER :: array(:)
         INTEGER, INTENT(in) :: n

         CPASSERT(n >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array) /= n) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n))
         END IF
      END SUBROUTINE cp_1d_r_guarantee_size

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n_rows the wanted number of rows
!> \param n_cols the wanted number of cols
!> \par History
!>      5.2001 first version [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_2d_r_guarantee_size(array, n_rows, n_cols)
         REAL(kind=dp), POINTER :: array(:, :)
         INTEGER, INTENT(in) :: n_rows, n_cols

         CPASSERT(n_cols >= 0)
         CPASSERT(n_rows >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array, 1) /= n_rows .OR. SIZE(array, 2) /= n_cols) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n_rows, n_cols))
         END IF
      END SUBROUTINE cp_2d_r_guarantee_size

! **************************************************************************************************
!> \brief returns the index at which the element el should be inserted in the
!>      array to keep it ordered (array(i)>=el).
!>      If the element is bigger than all the elements in the array returns
!>      the last index+1.
!> \param array the array to search
!> \param el the element to look for
!> \param l_index the lower index for binary search (defaults to 1)
!> \param u_index the upper index for binary search (defaults to size(array))
!> \return ...
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the array should be ordered in growing order
! **************************************************************************************************
      FUNCTION cp_1d_r_bsearch(array, el, l_index, u_index) &
         result(res)
         REAL(kind=dp), intent(in) :: array(:)
         REAL(kind=dp), intent(in) :: el
         INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
         integer :: res

         INTEGER :: lindex, uindex, aindex

         lindex = 1
         uindex = size(array)
         if (present(l_index)) lindex = l_index
         if (present(u_index)) uindex = u_index
         DO WHILE (lindex <= uindex)
            aindex = (lindex + uindex)/2
            IF (      array(aindex) < el) THEN
               lindex = aindex + 1
            ELSE
               uindex = aindex - 1
            END IF
         END DO
         res = lindex
      END FUNCTION cp_1d_r_bsearch
! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated paretized list
! **************************************************************************************************
      SUBROUTINE cp_1d_c_write(array, unit_nr, el_format)
         COMPLEX(KIND=dp), INTENT(in) :: array(:)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(es14.6)"

         WRITE (unit=unit_nr, fmt="('( ')", advance="no", iostat=iostat)
         CPASSERT(iostat == 0)
         IF (PRESENT(el_format)) THEN
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(SIZE(array))
         ELSE
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(SIZE(array))
         END IF
         WRITE (unit=unit_nr, fmt="(' )')")
         call m_flush(unit_nr)

      END SUBROUTINE cp_1d_c_write

! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated parentized list
! **************************************************************************************************
      SUBROUTINE cp_2d_c_write(array, unit_nr, el_format)
         COMPLEX(KIND=dp), INTENT(in) :: array(:, :)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(es14.6)"
         CHARACTER(len=200) :: fmtstr
         CHARACTER(len=25) :: nRiga

         nRiga = cp_to_string(SIZE(array, 2))
         DO i = 1, SIZE(array, 1)
            IF (PRESENT(el_format)) THEN
               fmtstr = '(" ",'//nRiga//el_format//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            ELSE
               fmtstr = '(" ",'//nRiga//defaultFormat//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            END IF
            CPASSERT(iostat == 0)
         END DO
         call m_flush(unit_nr)
      END SUBROUTINE cp_2d_c_write

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n the wanted size
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_1d_c_guarantee_size(array, n)
         COMPLEX(KIND=dp), POINTER :: array(:)
         INTEGER, INTENT(in) :: n

         CPASSERT(n >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array) /= n) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n))
         END IF
      END SUBROUTINE cp_1d_c_guarantee_size

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n_rows the wanted number of rows
!> \param n_cols the wanted number of cols
!> \par History
!>      5.2001 first version [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_2d_c_guarantee_size(array, n_rows, n_cols)
         COMPLEX(KIND=dp), POINTER :: array(:, :)
         INTEGER, INTENT(in) :: n_rows, n_cols

         CPASSERT(n_cols >= 0)
         CPASSERT(n_rows >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array, 1) /= n_rows .OR. SIZE(array, 2) /= n_cols) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n_rows, n_cols))
         END IF
      END SUBROUTINE cp_2d_c_guarantee_size

! **************************************************************************************************
!> \brief returns the index at which the element el should be inserted in the
!>      array to keep it ordered (array(i)>=el).
!>      If the element is bigger than all the elements in the array returns
!>      the last index+1.
!> \param array the array to search
!> \param el the element to look for
!> \param l_index the lower index for binary search (defaults to 1)
!> \param u_index the upper index for binary search (defaults to size(array))
!> \return ...
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the array should be ordered in growing order
! **************************************************************************************************
      FUNCTION cp_1d_c_bsearch(array, el, l_index, u_index) &
         result(res)
         COMPLEX(KIND=dp), intent(in) :: array(:)
         COMPLEX(KIND=dp), intent(in) :: el
         INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
         integer :: res

         INTEGER :: lindex, uindex, aindex

         lindex = 1
         uindex = size(array)
         if (present(l_index)) lindex = l_index
         if (present(u_index)) uindex = u_index
         DO WHILE (lindex <= uindex)
            aindex = (lindex + uindex)/2
            IF (REAL(array(aindex)) < REAL(el) .OR. (ABS(REAL(array(aindex))-REAL(el)) < epsilon(MAX(ABS(REAL(array(aindex))),&
                & ABS(REAL(el)))) .and. AIMAG(array(aindex)) <  AIMAG(el))) THEN
               lindex = aindex + 1
            ELSE
               uindex = aindex - 1
            END IF
         END DO
         res = lindex
      END FUNCTION cp_1d_c_bsearch
! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated paretized list
! **************************************************************************************************
      SUBROUTINE cp_1d_i_write(array, unit_nr, el_format)
         INTEGER(kind=int_4), INTENT(in) :: array(:)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(i6)"

         WRITE (unit=unit_nr, fmt="('( ')", advance="no", iostat=iostat)
         CPASSERT(iostat == 0)
         IF (PRESENT(el_format)) THEN
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(SIZE(array))
         ELSE
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(SIZE(array))
         END IF
         WRITE (unit=unit_nr, fmt="(' )')")
         call m_flush(unit_nr)

      END SUBROUTINE cp_1d_i_write

! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated parentized list
! **************************************************************************************************
      SUBROUTINE cp_2d_i_write(array, unit_nr, el_format)
         INTEGER(kind=int_4), INTENT(in) :: array(:, :)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(i6)"
         CHARACTER(len=200) :: fmtstr
         CHARACTER(len=25) :: nRiga

         nRiga = cp_to_string(SIZE(array, 2))
         DO i = 1, SIZE(array, 1)
            IF (PRESENT(el_format)) THEN
               fmtstr = '(" ",'//nRiga//el_format//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            ELSE
               fmtstr = '(" ",'//nRiga//defaultFormat//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            END IF
            CPASSERT(iostat == 0)
         END DO
         call m_flush(unit_nr)
      END SUBROUTINE cp_2d_i_write

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n the wanted size
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_1d_i_guarantee_size(array, n)
         INTEGER(kind=int_4), POINTER :: array(:)
         INTEGER, INTENT(in) :: n

         CPASSERT(n >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array) /= n) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n))
         END IF
      END SUBROUTINE cp_1d_i_guarantee_size

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n_rows the wanted number of rows
!> \param n_cols the wanted number of cols
!> \par History
!>      5.2001 first version [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_2d_i_guarantee_size(array, n_rows, n_cols)
         INTEGER(kind=int_4), POINTER :: array(:, :)
         INTEGER, INTENT(in) :: n_rows, n_cols

         CPASSERT(n_cols >= 0)
         CPASSERT(n_rows >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array, 1) /= n_rows .OR. SIZE(array, 2) /= n_cols) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n_rows, n_cols))
         END IF
      END SUBROUTINE cp_2d_i_guarantee_size

! **************************************************************************************************
!> \brief returns the index at which the element el should be inserted in the
!>      array to keep it ordered (array(i)>=el).
!>      If the element is bigger than all the elements in the array returns
!>      the last index+1.
!> \param array the array to search
!> \param el the element to look for
!> \param l_index the lower index for binary search (defaults to 1)
!> \param u_index the upper index for binary search (defaults to size(array))
!> \return ...
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the array should be ordered in growing order
! **************************************************************************************************
      FUNCTION cp_1d_i_bsearch(array, el, l_index, u_index) &
         result(res)
         INTEGER(kind=int_4), intent(in) :: array(:)
         INTEGER(kind=int_4), intent(in) :: el
         INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
         integer :: res

         INTEGER :: lindex, uindex, aindex

         lindex = 1
         uindex = size(array)
         if (present(l_index)) lindex = l_index
         if (present(u_index)) uindex = u_index
         DO WHILE (lindex <= uindex)
            aindex = (lindex + uindex)/2
            IF (      array(aindex) < el) THEN
               lindex = aindex + 1
            ELSE
               uindex = aindex - 1
            END IF
         END DO
         res = lindex
      END FUNCTION cp_1d_i_bsearch
! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated paretized list
! **************************************************************************************************
      SUBROUTINE cp_1d_logical_write(array, unit_nr, el_format)
         logical, INTENT(in) :: array(:)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(l1)"

         WRITE (unit=unit_nr, fmt="('( ')", advance="no", iostat=iostat)
         CPASSERT(iostat == 0)
         IF (PRESENT(el_format)) THEN
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=el_format, advance="no") array(SIZE(array))
         ELSE
            DO i = 1, SIZE(array) - 1
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(i)
               IF (MOD(i, 5) == 0) THEN  ! only a few elements per line
                  WRITE (unit=unit_nr, fmt="(',')")
               ELSE
                  WRITE (unit=unit_nr, fmt="(',')", advance="no")
               END IF
            END DO
            IF (SIZE(array) > 0) &
               WRITE (unit=unit_nr, fmt=defaultFormat, advance="no") array(SIZE(array))
         END IF
         WRITE (unit=unit_nr, fmt="(' )')")
         call m_flush(unit_nr)

      END SUBROUTINE cp_1d_logical_write

! **************************************************************************************************
!> \brief writes an array to the given unit
!> \param array the array to write
!> \param unit_nr the unit to write to (defaults to the standard out)
!> \param el_format the format of a single element
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      maybe I will move to a comma separated parentized list
! **************************************************************************************************
      SUBROUTINE cp_2d_logical_write(array, unit_nr, el_format)
         logical, INTENT(in) :: array(:, :)
         INTEGER, INTENT(in) :: unit_nr
         CHARACTER(len=*), INTENT(in), OPTIONAL :: el_format

         INTEGER :: iostat, i
         CHARACTER(len=*), PARAMETER :: defaultFormat = "(l1)"
         CHARACTER(len=200) :: fmtstr
         CHARACTER(len=25) :: nRiga

         nRiga = cp_to_string(SIZE(array, 2))
         DO i = 1, SIZE(array, 1)
            IF (PRESENT(el_format)) THEN
               fmtstr = '(" ",'//nRiga//el_format//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            ELSE
               fmtstr = '(" ",'//nRiga//defaultFormat//')'
               WRITE (unit=unit_nr, fmt=fmtstr, iostat=iostat) array(i, :)
            END IF
            CPASSERT(iostat == 0)
         END DO
         call m_flush(unit_nr)
      END SUBROUTINE cp_2d_logical_write

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n the wanted size
!> \par History
!>      12.2001 first version [fawzi]
!>      3.2002 templatized [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_1d_logical_guarantee_size(array, n)
         logical, POINTER :: array(:)
         INTEGER, INTENT(in) :: n

         CPASSERT(n >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array) /= n) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n))
         END IF
      END SUBROUTINE cp_1d_logical_guarantee_size

! **************************************************************************************************
!> \brief If the size of the array is changes reallocate it.
!>      Issues a warning when the size changes (but not on allocation
!>      and deallocation).
!>
!>      The data is NOT preserved (if you want to preserve the data see
!>      the realloc in the module memory_utilities)
!> \param array the array to reallocate if necessary
!> \param n_rows the wanted number of rows
!> \param n_cols the wanted number of cols
!> \par History
!>      5.2001 first version [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this is a different behaviour than the realloc in the module
!>      memory_utilities. It is quite low level
! **************************************************************************************************
      SUBROUTINE cp_2d_logical_guarantee_size(array, n_rows, n_cols)
         logical, POINTER :: array(:, :)
         INTEGER, INTENT(in) :: n_rows, n_cols

         CPASSERT(n_cols >= 0)
         CPASSERT(n_rows >= 0)
         IF (ASSOCIATED(array)) THEN
            IF (SIZE(array, 1) /= n_rows .OR. SIZE(array, 2) /= n_cols) THEN
               CPWARN('size has changed')
               DEALLOCATE (array)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            ALLOCATE (array(n_rows, n_cols))
         END IF
      END SUBROUTINE cp_2d_logical_guarantee_size

! **************************************************************************************************
!> \brief returns the index at which the element el should be inserted in the
!>      array to keep it ordered (array(i)>=el).
!>      If the element is bigger than all the elements in the array returns
!>      the last index+1.
!> \param array the array to search
!> \param el the element to look for
!> \param l_index the lower index for binary search (defaults to 1)
!> \param u_index the upper index for binary search (defaults to size(array))
!> \return ...
!> \par History
!>      06.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the array should be ordered in growing order
! **************************************************************************************************
      FUNCTION cp_1d_logical_bsearch(array, el, l_index, u_index) &
         result(res)
         logical, intent(in) :: array(:)
         logical, intent(in) :: el
         INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
         integer :: res

         INTEGER :: lindex, uindex, aindex

         lindex = 1
         uindex = size(array)
         if (present(l_index)) lindex = l_index
         if (present(u_index)) uindex = u_index
         DO WHILE (lindex <= uindex)
            aindex = (lindex + uindex)/2
            IF (      (.not. array(aindex)) .and. el) THEN
               lindex = aindex + 1
            ELSE
               uindex = aindex - 1
            END IF
         END DO
         res = lindex
      END FUNCTION cp_1d_logical_bsearch

END MODULE cp_array_utils
