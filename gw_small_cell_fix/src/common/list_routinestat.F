!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief An array-based list which grows on demand.
!>        When the internal array is full, a new array of twice the size will be
!>        allocated and the items are copied over.
!>
!>        This list can also be used as a stack.
!>        Have look at list_push(), list_pop() and list_peek().
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! **************************************************************************************************

MODULE list_routinestat
   USE timings_base_type, ONLY: callstack_entry_type, routine_report_type, routine_stat_type


#include "../base/base_uses.f90"
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: list_routinestat_type, list_routinestat_init, list_routinestat_push, list_routinestat_pop, list_routinestat_peek,&
          & list_routinestat_insert, list_routinestat_set, list_routinestat_get, list_routinestat_del, list_routinestat_clear,&
          & list_routinestat_size, list_routinestat_destroy, list_routinestat_isready

!this is an internal type
      TYPE private_item_type_routinestat
         PRIVATE
         TYPE(routine_stat_type), POINTER                   :: value => NULL()
      END TYPE private_item_type_routinestat

!this is an internal type
      TYPE private_item_p_type_routinestat
         PRIVATE
         TYPE(private_item_type_routinestat), POINTER :: p => Null()
      END TYPE private_item_p_type_routinestat

! this is the public type, which holds a list-instance
      TYPE list_routinestat_type
         PRIVATE
         TYPE(private_item_p_type_routinestat), DIMENSION(:), POINTER   :: arr => Null()
         INTEGER                                       :: size = -1
      END TYPE list_routinestat_type

      CONTAINS

! **************************************************************************************************
!> \brief Test if the given list has been initialized.
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_routinestat_isready(list) RESULT(res)
         TYPE(list_routinestat_type), intent(in)  :: list
         LOGICAL                                     :: res
         res = ASSOCIATED(list%arr)
      END FUNCTION list_routinestat_isready

! **************************************************************************************************
!> \brief Allocates the internal data-structures of the given list.
!>        This has to be called before any of the other routines.
!>        For deallocation call list_[valuetype]_destroy.
!> \param list ...
!> \param initial_capacity The initial size of the internal array (default=11).
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_init(list, initial_capacity)
         TYPE(list_routinestat_type), intent(inout)  :: list
         INTEGER, INTENT(in), OPTIONAL               :: initial_capacity
         INTEGER                                     :: stat
         INTEGER                                     :: initial_capacity_

         initial_capacity_ = 11
         If (PRESENT(initial_capacity)) initial_capacity_ = initial_capacity

         IF (initial_capacity_ < 0) &
            CPABORT("list_routinestat_create: initial_capacity < 0")

         IF (ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_create: list is already initialized.")

         ALLOCATE (list%arr(initial_capacity_), stat=stat)
         IF (stat /= 0) &
            CPABORT("list_routinestat_init: allocation failed")

         list%size = 0
      END SUBROUTINE list_routinestat_init

! **************************************************************************************************
!> \brief Deallocated the internal data-structures of the given list.
!>        Caution: If the stored values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_destroy(list)
         TYPE(list_routinestat_type), intent(inout)  :: list
         INTEGER :: i
         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_destroy: list is not initialized.")

         do i = 1, list%size
            deallocate (list%arr(i)%p)
         end do
         deallocate (list%arr)
         list%size = -1
      END SUBROUTINE list_routinestat_destroy

! **************************************************************************************************
!> \brief Assings the given value to the given position in the list.
!>        Thereby, the former value at that position gets overwritten.
!>        If the position is out of bounds, the program stops.
!> \param list ...
!> \param value ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1.
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_set(list, value, pos)
         TYPE(list_routinestat_type), intent(inout)  :: list
         TYPE(routine_stat_type), POINTER, intent(in) :: value
         INTEGER, intent(in) :: pos
         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_set: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_routinestat_set: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_routinestat_set: pos > size")
         list%arr(pos)%p%value =>value
      END SUBROUTINE list_routinestat_set

! **************************************************************************************************
!> \brief Appends the given value at the end of the list.
!> \param list ...
!> \param value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_push(list, value)
         TYPE(list_routinestat_type), intent(inout)  :: list
         TYPE(routine_stat_type), POINTER, intent(in)                  :: value
         INTEGER                                     :: stat

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_push: list is not initialized.")
         if (list%size == size(list%arr)) &
            call change_capacity_routinestat (list, 2*size(list%arr) + 1)

         list%size = list%size + 1
         ALLOCATE (list%arr(list%size)%p, stat=stat)
         IF (stat /= 0) &
            CPABORT("list_routinestat_push: allocation failed")
         list%arr(list%size)%p%value =>value
      END SUBROUTINE list_routinestat_push

! **************************************************************************************************
!> \brief Inserts the given value at the givenn position within the list.
!>        Values which lay behind the insertion-position move one position up.
!> \param list ...
!> \param value ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+2 .
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_insert(list, value, pos)
         TYPE(list_routinestat_type), intent(inout)  :: list
         TYPE(routine_stat_type), POINTER, intent(in) :: value
         INTEGER, intent(in) :: pos
         INTEGER :: i, stat

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_insert: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_routinestat_insert: pos < 1")
         IF (pos > list%size + 1) &
            CPABORT("list_routinestat_insert: pos > size+1")

         if (list%size == size(list%arr)) &
            call change_capacity_routinestat (list, 2*size(list%arr) + 1)

         list%size = list%size + 1
         do i = list%size, pos + 1, -1
            list%arr(i)%p => list%arr(i - 1)%p
         end do

         ALLOCATE (list%arr(pos)%p, stat=stat)
         IF (stat /= 0) &
            CPABORT("list_routinestat_insert: allocation failed.")
         list%arr(pos)%p%value =>value
      END SUBROUTINE list_routinestat_insert

! **************************************************************************************************
!> \brief Returns the last element in the list.
!>    Is equivalent to: list_routinestat_get(list, list_routinestat_size(list))
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_routinestat_peek(list) RESULT(value)
         TYPE(list_routinestat_type), intent(inout)  :: list
         TYPE(routine_stat_type), POINTER  :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_peek: list is not initialized.")
         IF (list%size < 1) &
            CPABORT("list_routinestat_peek: list is empty.")

         value =>list%arr(list%size)%p%value
      END FUNCTION list_routinestat_peek

! **************************************************************************************************
!> \brief Returns the last element in the list and removes it.
!>        Is equivialent to:
!>        value = list_routinestat_get(list, list_routinestat_size(list))
!>            call list_routinestat_del(list, list_routinestat_size(list))
!>
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_routinestat_pop(list) RESULT(value)
         TYPE(list_routinestat_type), intent(inout)  :: list
         TYPE(routine_stat_type), POINTER  :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_pop: list is not initialized.")
         IF (list%size < 1) &
            CPABORT("list_routinestat_pop: list is empty.")

         value =>list%arr(list%size)%p%value
         deallocate (list%arr(list%size)%p)
         list%size = list%size - 1
      END FUNCTION list_routinestat_pop

! **************************************************************************************************
!> \brief Removes all values from the list. The list itself is not deallocated.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_clear(list)
         TYPE(list_routinestat_type), intent(inout)  :: list
         INTEGER :: i

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_clear: list is not initialized.")

         do i = 1, list%size
            deallocate (list%arr(i)%p)
         end do
         list%size = 0
      END SUBROUTINE list_routinestat_clear

!
! **************************************************************************************************
!> \brief Returns the value at the given position from the list.
!> \param list ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1 .
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_routinestat_get(list, pos) RESULT(value)
         TYPE(list_routinestat_type), intent(in)  :: list
         INTEGER, intent(in) :: pos
         TYPE(routine_stat_type), POINTER :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_get: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_routinestat_get: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_routinestat_get: pos > size")

         value =>list%arr(pos)%p%value

      END FUNCTION list_routinestat_get

! **************************************************************************************************
!> \brief Removes the value at the given position from the list.
!> \param list ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1 .
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_routinestat_del(list, pos)
         TYPE(list_routinestat_type), intent(inout)  :: list
         INTEGER, intent(in) :: pos
         INTEGER :: i

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_del: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_routinestat_det: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_routinestat_det: pos > size")

         deallocate (list%arr(pos)%p)
         do i = pos, list%size - 1
            list%arr(i)%p => list%arr(i + 1)%p
         end do

         list%size = list%size - 1

      END SUBROUTINE list_routinestat_del

! **************************************************************************************************
!> \brief Returns the current size of the list.
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_routinestat_size(list) RESULT(size)
         TYPE(list_routinestat_type), intent(in)  :: list
         INTEGER :: size

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_routinestat_size: list is not initialized.")

         size = list%size
      END FUNCTION list_routinestat_size

! **************************************************************************************************
!> \brief Internal routine for changing the size of the internal array.
!> \param list ...
!> \param new_capacity ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE change_capacity_routinestat (list, new_capacity)
         TYPE(list_routinestat_type), intent(inout)  :: list
         INTEGER, intent(in) :: new_capacity
         INTEGER :: i, new_cap, stat
         TYPE(private_item_p_type_routinestat), DIMENSION(:), POINTER :: old_arr

         new_cap = new_capacity
         IF (new_cap < 0) &
            CPABORT("list_routinestat_change_capacity: new_capacity < 0")
         IF (new_cap < list%size) &
            CPABORT("list_routinestat_change_capacity: new_capacity < size")
         IF (new_cap > HUGE(i)) THEN
            IF (size(list%arr) == HUGE(i)) &
               CPABORT("list_routinestat_change_capacity: list has reached integer limit.")
            new_cap = HUGE(i) ! grow as far as possible
         END IF

         old_arr => list%arr
         allocate (list%arr(new_cap), stat=stat)
         IF (stat /= 0) &
            CPABORT("list_routinestat_change_capacity: allocation failed")

         do i = 1, list%size
            allocate (list%arr(i)%p, stat=stat)
            IF (stat /= 0) &
               CPABORT("list_routinestat_change_capacity: allocation failed")
            list%arr(i)%p%value =>old_arr(i)%p%value
            deallocate (old_arr(i)%p)
         end do
         deallocate (old_arr)

      END SUBROUTINE change_capacity_routinestat

END MODULE list_routinestat
