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

MODULE list_callstackentry
   USE timings_base_type, ONLY: callstack_entry_type, routine_report_type, routine_stat_type


#include "../base/base_uses.f90"
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: list_callstackentry_type, list_callstackentry_init, list_callstackentry_push, list_callstackentry_pop,&
          & list_callstackentry_peek, list_callstackentry_insert, list_callstackentry_set, list_callstackentry_get,&
          & list_callstackentry_del, list_callstackentry_clear, list_callstackentry_size, list_callstackentry_destroy,&
          & list_callstackentry_isready

!this is an internal type
      TYPE private_item_type_callstackentry
         PRIVATE
         TYPE(callstack_entry_type)                   :: value = callstack_entry_type()
      END TYPE private_item_type_callstackentry

!this is an internal type
      TYPE private_item_p_type_callstackentry
         PRIVATE
         TYPE(private_item_type_callstackentry), POINTER :: p => Null()
      END TYPE private_item_p_type_callstackentry

! this is the public type, which holds a list-instance
      TYPE list_callstackentry_type
         PRIVATE
         TYPE(private_item_p_type_callstackentry), DIMENSION(:), POINTER   :: arr => Null()
         INTEGER                                       :: size = -1
      END TYPE list_callstackentry_type

      CONTAINS

! **************************************************************************************************
!> \brief Test if the given list has been initialized.
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_callstackentry_isready(list) RESULT(res)
         TYPE(list_callstackentry_type), intent(in)  :: list
         LOGICAL                                     :: res
         res = ASSOCIATED(list%arr)
      END FUNCTION list_callstackentry_isready

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
      SUBROUTINE list_callstackentry_init(list, initial_capacity)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         INTEGER, INTENT(in), OPTIONAL               :: initial_capacity
         INTEGER                                     :: stat
         INTEGER                                     :: initial_capacity_

         initial_capacity_ = 11
         If (PRESENT(initial_capacity)) initial_capacity_ = initial_capacity

         IF (initial_capacity_ < 0) &
            CPABORT("list_callstackentry_create: initial_capacity < 0")

         IF (ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_create: list is already initialized.")

         ALLOCATE (list%arr(initial_capacity_), stat=stat)
         IF (stat /= 0) &
            CPABORT("list_callstackentry_init: allocation failed")

         list%size = 0
      END SUBROUTINE list_callstackentry_init

! **************************************************************************************************
!> \brief Deallocated the internal data-structures of the given list.
!>        Caution: If the stored values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_callstackentry_destroy(list)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         INTEGER :: i
         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_destroy: list is not initialized.")

         do i = 1, list%size
            deallocate (list%arr(i)%p)
         end do
         deallocate (list%arr)
         list%size = -1
      END SUBROUTINE list_callstackentry_destroy

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
      SUBROUTINE list_callstackentry_set(list, value, pos)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         TYPE(callstack_entry_type), intent(in) :: value
         INTEGER, intent(in) :: pos
         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_set: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_callstackentry_set: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_callstackentry_set: pos > size")
         list%arr(pos)%p%value =value
      END SUBROUTINE list_callstackentry_set

! **************************************************************************************************
!> \brief Appends the given value at the end of the list.
!> \param list ...
!> \param value ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_callstackentry_push(list, value)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         TYPE(callstack_entry_type), intent(in)                  :: value
         INTEGER                                     :: stat

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_push: list is not initialized.")
         if (list%size == size(list%arr)) &
            call change_capacity_callstackentry (list, 2*size(list%arr) + 1)

         list%size = list%size + 1
         ALLOCATE (list%arr(list%size)%p, stat=stat)
         IF (stat /= 0) &
            CPABORT("list_callstackentry_push: allocation failed")
         list%arr(list%size)%p%value =value
      END SUBROUTINE list_callstackentry_push

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
      SUBROUTINE list_callstackentry_insert(list, value, pos)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         TYPE(callstack_entry_type), intent(in) :: value
         INTEGER, intent(in) :: pos
         INTEGER :: i, stat

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_insert: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_callstackentry_insert: pos < 1")
         IF (pos > list%size + 1) &
            CPABORT("list_callstackentry_insert: pos > size+1")

         if (list%size == size(list%arr)) &
            call change_capacity_callstackentry (list, 2*size(list%arr) + 1)

         list%size = list%size + 1
         do i = list%size, pos + 1, -1
            list%arr(i)%p => list%arr(i - 1)%p
         end do

         ALLOCATE (list%arr(pos)%p, stat=stat)
         IF (stat /= 0) &
            CPABORT("list_callstackentry_insert: allocation failed.")
         list%arr(pos)%p%value =value
      END SUBROUTINE list_callstackentry_insert

! **************************************************************************************************
!> \brief Returns the last element in the list.
!>    Is equivalent to: list_callstackentry_get(list, list_callstackentry_size(list))
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_callstackentry_peek(list) RESULT(value)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         TYPE(callstack_entry_type)  :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_peek: list is not initialized.")
         IF (list%size < 1) &
            CPABORT("list_callstackentry_peek: list is empty.")

         value =list%arr(list%size)%p%value
      END FUNCTION list_callstackentry_peek

! **************************************************************************************************
!> \brief Returns the last element in the list and removes it.
!>        Is equivialent to:
!>        value = list_callstackentry_get(list, list_callstackentry_size(list))
!>            call list_callstackentry_del(list, list_callstackentry_size(list))
!>
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_callstackentry_pop(list) RESULT(value)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         TYPE(callstack_entry_type)  :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_pop: list is not initialized.")
         IF (list%size < 1) &
            CPABORT("list_callstackentry_pop: list is empty.")

         value =list%arr(list%size)%p%value
         deallocate (list%arr(list%size)%p)
         list%size = list%size - 1
      END FUNCTION list_callstackentry_pop

! **************************************************************************************************
!> \brief Removes all values from the list. The list itself is not deallocated.
!> \param list ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_callstackentry_clear(list)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         INTEGER :: i

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_clear: list is not initialized.")

         do i = 1, list%size
            deallocate (list%arr(i)%p)
         end do
         list%size = 0
      END SUBROUTINE list_callstackentry_clear

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
      FUNCTION list_callstackentry_get(list, pos) RESULT(value)
         TYPE(list_callstackentry_type), intent(in)  :: list
         INTEGER, intent(in) :: pos
         TYPE(callstack_entry_type) :: value

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_get: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_callstackentry_get: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_callstackentry_get: pos > size")

         value =list%arr(pos)%p%value

      END FUNCTION list_callstackentry_get

! **************************************************************************************************
!> \brief Removes the value at the given position from the list.
!> \param list ...
!> \param pos Position in the list - musst fulfill 0 < pos < list_size+1 .
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE list_callstackentry_del(list, pos)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         INTEGER, intent(in) :: pos
         INTEGER :: i

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_del: list is not initialized.")
         IF (pos < 1) &
            CPABORT("list_callstackentry_det: pos < 1")
         IF (pos > list%size) &
            CPABORT("list_callstackentry_det: pos > size")

         deallocate (list%arr(pos)%p)
         do i = pos, list%size - 1
            list%arr(i)%p => list%arr(i + 1)%p
         end do

         list%size = list%size - 1

      END SUBROUTINE list_callstackentry_del

! **************************************************************************************************
!> \brief Returns the current size of the list.
!> \param list ...
!> \return ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION list_callstackentry_size(list) RESULT(size)
         TYPE(list_callstackentry_type), intent(in)  :: list
         INTEGER :: size

         IF (.not. ASSOCIATED(list%arr)) &
            CPABORT("list_callstackentry_size: list is not initialized.")

         size = list%size
      END FUNCTION list_callstackentry_size

! **************************************************************************************************
!> \brief Internal routine for changing the size of the internal array.
!> \param list ...
!> \param new_capacity ...
!> \par History
!>      12.2012 created [ole]
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE change_capacity_callstackentry (list, new_capacity)
         TYPE(list_callstackentry_type), intent(inout)  :: list
         INTEGER, intent(in) :: new_capacity
         INTEGER :: i, new_cap, stat
         TYPE(private_item_p_type_callstackentry), DIMENSION(:), POINTER :: old_arr

         new_cap = new_capacity
         IF (new_cap < 0) &
            CPABORT("list_callstackentry_change_capacity: new_capacity < 0")
         IF (new_cap < list%size) &
            CPABORT("list_callstackentry_change_capacity: new_capacity < size")
         IF (new_cap > HUGE(i)) THEN
            IF (size(list%arr) == HUGE(i)) &
               CPABORT("list_callstackentry_change_capacity: list has reached integer limit.")
            new_cap = HUGE(i) ! grow as far as possible
         END IF

         old_arr => list%arr
         allocate (list%arr(new_cap), stat=stat)
         IF (stat /= 0) &
            CPABORT("list_callstackentry_change_capacity: allocation failed")

         do i = 1, list%size
            allocate (list%arr(i)%p, stat=stat)
            IF (stat /= 0) &
               CPABORT("list_callstackentry_change_capacity: allocation failed")
            list%arr(i)%p%value =old_arr(i)%p%value
            deallocate (old_arr(i)%p)
         end do
         deallocate (old_arr)

      END SUBROUTINE change_capacity_callstackentry

END MODULE list_callstackentry
