!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!


MODULE cp_linked_list_pw
   USE kinds, ONLY: dp
   USE cp_array_utils, ONLY: cp_3d_r_cp_type, cp_3d_c_cp_type, cp_1d_r_cp_type, cp_1d_c_cp_type
   USE realspace_grid_types, ONLY: realspace_grid_type, realspace_grid_p_type
#include "../base/base_uses.f90"


! **************************************************************************************************
!> \brief describes a generic linked list template.
!>      Linked list are supposed to always use pointers to the nodes for
!>      basically everything: a pointer to the node is a list, an element of
!>      the list, an iterator between the elment of the list.
!>      An empty list is represented by an unassociated pointer.
!> \par History
!>      1.2002 created
!>      4.2002 changed most of the methods, by making access to the list
!>             always through pointers (identifying pointer, list and iterators)
!>      6.2004 removed %initialized from list elements
!> \author Fawzi Mohamed
! **************************************************************************************************
      IMPLICIT NONE
      PRIVATE

      LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.


!API type
         PUBLIC :: cp_sll_3d_r_type, cp_sll_3d_r_p_type
         PUBLIC :: cp_sll_3d_c_type, cp_sll_3d_c_p_type
         PUBLIC :: cp_sll_1d_r_type, cp_sll_1d_r_p_type
         PUBLIC :: cp_sll_1d_c_type, cp_sll_1d_c_p_type
         PUBLIC :: cp_sll_rs_type, cp_sll_rs_p_type

!API common methods
      PUBLIC :: cp_create, cp_dealloc, cp_next
!API special get
      PUBLIC :: cp_get_first_el, cp_get_rest, cp_get_length, cp_get_element_at, cp_to_array
!API special set
      PUBLIC :: cp_set_element_at
!API structure manipulation
      PUBLIC :: cp_insert, cp_remove_first_el, cp_remove_el, cp_remove_all

      !API underlying routines

         PUBLIC :: cp_sll_3d_r_create, cp_sll_3d_r_dealloc, &
                   cp_sll_3d_r_next, &
                   cp_sll_3d_r_get_length, &
                   cp_sll_3d_r_get_el_at, cp_sll_3d_r_set_el_at, &
                   cp_sll_3d_r_get_first_el, cp_sll_3d_r_get_rest, &
                   cp_sll_3d_r_insert_el, cp_sll_3d_r_insert_el_at, &
                   cp_sll_3d_r_rm_first_el, cp_sll_3d_r_rm_el_at, &
                   cp_sll_3d_r_rm_all_el, &
                   cp_sll_3d_r_to_array
         PUBLIC :: cp_sll_3d_c_create, cp_sll_3d_c_dealloc, &
                   cp_sll_3d_c_next, &
                   cp_sll_3d_c_get_length, &
                   cp_sll_3d_c_get_el_at, cp_sll_3d_c_set_el_at, &
                   cp_sll_3d_c_get_first_el, cp_sll_3d_c_get_rest, &
                   cp_sll_3d_c_insert_el, cp_sll_3d_c_insert_el_at, &
                   cp_sll_3d_c_rm_first_el, cp_sll_3d_c_rm_el_at, &
                   cp_sll_3d_c_rm_all_el, &
                   cp_sll_3d_c_to_array
         PUBLIC :: cp_sll_1d_r_create, cp_sll_1d_r_dealloc, &
                   cp_sll_1d_r_next, &
                   cp_sll_1d_r_get_length, &
                   cp_sll_1d_r_get_el_at, cp_sll_1d_r_set_el_at, &
                   cp_sll_1d_r_get_first_el, cp_sll_1d_r_get_rest, &
                   cp_sll_1d_r_insert_el, cp_sll_1d_r_insert_el_at, &
                   cp_sll_1d_r_rm_first_el, cp_sll_1d_r_rm_el_at, &
                   cp_sll_1d_r_rm_all_el, &
                   cp_sll_1d_r_to_array
         PUBLIC :: cp_sll_1d_c_create, cp_sll_1d_c_dealloc, &
                   cp_sll_1d_c_next, &
                   cp_sll_1d_c_get_length, &
                   cp_sll_1d_c_get_el_at, cp_sll_1d_c_set_el_at, &
                   cp_sll_1d_c_get_first_el, cp_sll_1d_c_get_rest, &
                   cp_sll_1d_c_insert_el, cp_sll_1d_c_insert_el_at, &
                   cp_sll_1d_c_rm_first_el, cp_sll_1d_c_rm_el_at, &
                   cp_sll_1d_c_rm_all_el, &
                   cp_sll_1d_c_to_array
         PUBLIC :: cp_sll_rs_create, cp_sll_rs_dealloc, &
                   cp_sll_rs_next, &
                   cp_sll_rs_get_length, &
                   cp_sll_rs_get_el_at, cp_sll_rs_set_el_at, &
                   cp_sll_rs_get_first_el, cp_sll_rs_get_rest, &
                   cp_sll_rs_insert_el, cp_sll_rs_insert_el_at, &
                   cp_sll_rs_rm_first_el, cp_sll_rs_rm_el_at, &
                   cp_sll_rs_rm_all_el, &
                   cp_sll_rs_to_array

! creation of an object (from a pointer)
      INTERFACE cp_create
            MODULE PROCEDURE cp_sll_3d_r_create
            MODULE PROCEDURE cp_sll_3d_c_create
            MODULE PROCEDURE cp_sll_1d_r_create
            MODULE PROCEDURE cp_sll_1d_c_create
            MODULE PROCEDURE cp_sll_rs_create
      END INTERFACE
! destruction of an object (from a pointer)
      INTERFACE cp_dealloc
            MODULE PROCEDURE cp_sll_3d_r_dealloc
            MODULE PROCEDURE cp_sll_3d_c_dealloc
            MODULE PROCEDURE cp_sll_1d_r_dealloc
            MODULE PROCEDURE cp_sll_1d_c_dealloc
            MODULE PROCEDURE cp_sll_rs_dealloc
      END INTERFACE
! iterates to the next element
      INTERFACE cp_next
            MODULE PROCEDURE cp_sll_3d_r_next
            MODULE PROCEDURE cp_sll_3d_c_next
            MODULE PROCEDURE cp_sll_1d_r_next
            MODULE PROCEDURE cp_sll_1d_c_next
            MODULE PROCEDURE cp_sll_rs_next
      END INTERFACE
! returns the first element
      INTERFACE cp_get_first_el
            MODULE PROCEDURE cp_sll_3d_r_get_first_el
            MODULE PROCEDURE cp_sll_3d_c_get_first_el
            MODULE PROCEDURE cp_sll_1d_r_get_first_el
            MODULE PROCEDURE cp_sll_1d_c_get_first_el
            MODULE PROCEDURE cp_sll_rs_get_first_el
      END INTERFACE
! returns the rest of the list
      INTERFACE cp_get_rest
            MODULE PROCEDURE cp_sll_3d_r_get_rest
            MODULE PROCEDURE cp_sll_3d_c_get_rest
            MODULE PROCEDURE cp_sll_1d_r_get_rest
            MODULE PROCEDURE cp_sll_1d_c_get_rest
            MODULE PROCEDURE cp_sll_rs_get_rest
      END INTERFACE
! returns the length of the list
      INTERFACE cp_get_length
            MODULE PROCEDURE cp_sll_3d_r_get_length
            MODULE PROCEDURE cp_sll_3d_c_get_length
            MODULE PROCEDURE cp_sll_1d_r_get_length
            MODULE PROCEDURE cp_sll_1d_c_get_length
            MODULE PROCEDURE cp_sll_rs_get_length
      END INTERFACE
! returns the element at the given position
      INTERFACE cp_get_element_at
            MODULE PROCEDURE cp_sll_3d_r_get_el_at
            MODULE PROCEDURE cp_sll_3d_c_get_el_at
            MODULE PROCEDURE cp_sll_1d_r_get_el_at
            MODULE PROCEDURE cp_sll_1d_c_get_el_at
            MODULE PROCEDURE cp_sll_rs_get_el_at
      END INTERFACE
! sets the element at the given position
      INTERFACE cp_set_element_at
            MODULE PROCEDURE cp_sll_3d_r_set_el_at
            MODULE PROCEDURE cp_sll_3d_c_set_el_at
            MODULE PROCEDURE cp_sll_1d_r_set_el_at
            MODULE PROCEDURE cp_sll_1d_c_set_el_at
            MODULE PROCEDURE cp_sll_rs_set_el_at
      END INTERFACE
! inserts one element call cp_insert(list,element,...)
      INTERFACE cp_insert
            MODULE PROCEDURE cp_sll_3d_r_insert_el
            MODULE PROCEDURE cp_sll_3d_c_insert_el
            MODULE PROCEDURE cp_sll_1d_r_insert_el
            MODULE PROCEDURE cp_sll_1d_c_insert_el
            MODULE PROCEDURE cp_sll_rs_insert_el
      END INTERFACE
      INTERFACE cp_insert_at
            MODULE PROCEDURE cp_sll_3d_r_insert_el_at
            MODULE PROCEDURE cp_sll_3d_c_insert_el_at
            MODULE PROCEDURE cp_sll_1d_r_insert_el_at
            MODULE PROCEDURE cp_sll_1d_c_insert_el_at
            MODULE PROCEDURE cp_sll_rs_insert_el_at
      END INTERFACE
! removes an element
      INTERFACE cp_remove_el
            MODULE PROCEDURE cp_sll_3d_r_rm_first_el, &
               cp_sll_3d_r_rm_el_at
            MODULE PROCEDURE cp_sll_3d_c_rm_first_el, &
               cp_sll_3d_c_rm_el_at
            MODULE PROCEDURE cp_sll_1d_r_rm_first_el, &
               cp_sll_1d_r_rm_el_at
            MODULE PROCEDURE cp_sll_1d_c_rm_first_el, &
               cp_sll_1d_c_rm_el_at
            MODULE PROCEDURE cp_sll_rs_rm_first_el, &
               cp_sll_rs_rm_el_at
      END INTERFACE
! removes the first el
      INTERFACE cp_remove_first_el
            MODULE PROCEDURE cp_sll_3d_r_rm_first_el
            MODULE PROCEDURE cp_sll_3d_c_rm_first_el
            MODULE PROCEDURE cp_sll_1d_r_rm_first_el
            MODULE PROCEDURE cp_sll_1d_c_rm_first_el
            MODULE PROCEDURE cp_sll_rs_rm_first_el
      END INTERFACE
! remove all the elements
      INTERFACE cp_remove_all
            MODULE PROCEDURE cp_sll_3d_r_rm_all_el
            MODULE PROCEDURE cp_sll_3d_c_rm_all_el
            MODULE PROCEDURE cp_sll_1d_r_rm_all_el
            MODULE PROCEDURE cp_sll_1d_c_rm_all_el
            MODULE PROCEDURE cp_sll_rs_rm_all_el
      END INTERFACE
! transorms the list in array
      INTERFACE cp_to_array
            module procedure cp_sll_3d_r_to_array
            module procedure cp_sll_3d_c_to_array
            module procedure cp_sll_1d_r_to_array
            module procedure cp_sll_1d_c_to_array
            module procedure cp_sll_rs_to_array
      END INTERFACE

! **************************************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_3d_r_type
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous :: first_el  => NULL()
            TYPE(cp_sll_3d_r_type), POINTER :: rest => NULL()
         END TYPE cp_sll_3d_r_type
! **************************************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_3d_c_type
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous :: first_el  => NULL()
            TYPE(cp_sll_3d_c_type), POINTER :: rest => NULL()
         END TYPE cp_sll_3d_c_type
! **************************************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_1d_r_type
            REAL(kind=dp),dimension(:),pointer,contiguous :: first_el  => NULL()
            TYPE(cp_sll_1d_r_type), POINTER :: rest => NULL()
         END TYPE cp_sll_1d_r_type
! **************************************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_1d_c_type
            COMPLEX(kind=dp),dimension(:),pointer,contiguous :: first_el  => NULL()
            TYPE(cp_sll_1d_c_type), POINTER :: rest => NULL()
         END TYPE cp_sll_1d_c_type
! **************************************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_rs_type
            type(realspace_grid_type),pointer :: first_el  => NULL()
            TYPE(cp_sll_rs_type), POINTER :: rest => NULL()
         END TYPE cp_sll_rs_type

! **************************************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_3d_r_p_type
            TYPE(cp_sll_3d_r_type), POINTER :: list => NULL()
         END TYPE cp_sll_3d_r_p_type
! **************************************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_3d_c_p_type
            TYPE(cp_sll_3d_c_type), POINTER :: list => NULL()
         END TYPE cp_sll_3d_c_p_type
! **************************************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_1d_r_p_type
            TYPE(cp_sll_1d_r_type), POINTER :: list => NULL()
         END TYPE cp_sll_1d_r_p_type
! **************************************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_1d_c_p_type
            TYPE(cp_sll_1d_c_type), POINTER :: list => NULL()
         END TYPE cp_sll_1d_c_p_type
! **************************************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         TYPE cp_sll_rs_p_type
            TYPE(cp_sll_rs_type), POINTER :: list => NULL()
         END TYPE cp_sll_rs_p_type

      CONTAINS

       ! this comment prevents trailing whitespace

! =========== creation / destruction ========

! **************************************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_create(sll, first_el, rest)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous, OPTIONAL :: first_el
            TYPE(cp_sll_3d_r_type), POINTER, OPTIONAL :: rest

            IF (.NOT. PRESENT(first_el)) THEN
               NULLIFY (sll)
               IF (PRESENT(rest)) sll => rest
            ELSE
               ALLOCATE (sll)
               sll%first_el =>first_el
               NULLIFY (sll%rest)
               IF (PRESENT(rest)) sll%rest => rest
            END IF
         END SUBROUTINE cp_sll_3d_r_create
! **************************************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_create(sll, first_el, rest)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous, OPTIONAL :: first_el
            TYPE(cp_sll_3d_c_type), POINTER, OPTIONAL :: rest

            IF (.NOT. PRESENT(first_el)) THEN
               NULLIFY (sll)
               IF (PRESENT(rest)) sll => rest
            ELSE
               ALLOCATE (sll)
               sll%first_el =>first_el
               NULLIFY (sll%rest)
               IF (PRESENT(rest)) sll%rest => rest
            END IF
         END SUBROUTINE cp_sll_3d_c_create
! **************************************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_create(sll, first_el, rest)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:),pointer,contiguous, OPTIONAL :: first_el
            TYPE(cp_sll_1d_r_type), POINTER, OPTIONAL :: rest

            IF (.NOT. PRESENT(first_el)) THEN
               NULLIFY (sll)
               IF (PRESENT(rest)) sll => rest
            ELSE
               ALLOCATE (sll)
               sll%first_el =>first_el
               NULLIFY (sll%rest)
               IF (PRESENT(rest)) sll%rest => rest
            END IF
         END SUBROUTINE cp_sll_1d_r_create
! **************************************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_create(sll, first_el, rest)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:),pointer,contiguous, OPTIONAL :: first_el
            TYPE(cp_sll_1d_c_type), POINTER, OPTIONAL :: rest

            IF (.NOT. PRESENT(first_el)) THEN
               NULLIFY (sll)
               IF (PRESENT(rest)) sll => rest
            ELSE
               ALLOCATE (sll)
               sll%first_el =>first_el
               NULLIFY (sll%rest)
               IF (PRESENT(rest)) sll%rest => rest
            END IF
         END SUBROUTINE cp_sll_1d_c_create
! **************************************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_create(sll, first_el, rest)
            TYPE(cp_sll_rs_type), POINTER :: sll
            type(realspace_grid_type),pointer, OPTIONAL :: first_el
            TYPE(cp_sll_rs_type), POINTER, OPTIONAL :: rest

            IF (.NOT. PRESENT(first_el)) THEN
               NULLIFY (sll)
               IF (PRESENT(rest)) sll => rest
            ELSE
               ALLOCATE (sll)
               sll%first_el =>first_el
               NULLIFY (sll%rest)
               IF (PRESENT(rest)) sll%rest => rest
            END IF
         END SUBROUTINE cp_sll_rs_create

! **************************************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elements that are stored in the list
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_dealloc(sll)
            TYPE(cp_sll_3d_r_type), POINTER :: sll

            CALL cp_sll_3d_r_rm_all_el(sll)
         END SUBROUTINE cp_sll_3d_r_dealloc
! **************************************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elements that are stored in the list
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_dealloc(sll)
            TYPE(cp_sll_3d_c_type), POINTER :: sll

            CALL cp_sll_3d_c_rm_all_el(sll)
         END SUBROUTINE cp_sll_3d_c_dealloc
! **************************************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elements that are stored in the list
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_dealloc(sll)
            TYPE(cp_sll_1d_r_type), POINTER :: sll

            CALL cp_sll_1d_r_rm_all_el(sll)
         END SUBROUTINE cp_sll_1d_r_dealloc
! **************************************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elements that are stored in the list
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_dealloc(sll)
            TYPE(cp_sll_1d_c_type), POINTER :: sll

            CALL cp_sll_1d_c_rm_all_el(sll)
         END SUBROUTINE cp_sll_1d_c_dealloc
! **************************************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elements that are stored in the list
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_dealloc(sll)
            TYPE(cp_sll_rs_type), POINTER :: sll

            CALL cp_sll_rs_rm_all_el(sll)
         END SUBROUTINE cp_sll_rs_dealloc

! * low-level *

! **************************************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_dealloc_node(sll)
            TYPE(cp_sll_3d_r_type), POINTER :: sll

            DEALLOCATE (sll)
         END SUBROUTINE cp_sll_3d_r_dealloc_node
! **************************************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_dealloc_node(sll)
            TYPE(cp_sll_3d_c_type), POINTER :: sll

            DEALLOCATE (sll)
         END SUBROUTINE cp_sll_3d_c_dealloc_node
! **************************************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_dealloc_node(sll)
            TYPE(cp_sll_1d_r_type), POINTER :: sll

            DEALLOCATE (sll)
         END SUBROUTINE cp_sll_1d_r_dealloc_node
! **************************************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_dealloc_node(sll)
            TYPE(cp_sll_1d_c_type), POINTER :: sll

            DEALLOCATE (sll)
         END SUBROUTINE cp_sll_1d_c_dealloc_node
! **************************************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_dealloc_node(sll)
            TYPE(cp_sll_rs_type), POINTER :: sll

            DEALLOCATE (sll)
         END SUBROUTINE cp_sll_rs_dealloc_node

! ============= get/set ============

! **************************************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_get_first_el(sll) RESULT(res)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous :: res

            res =>sll%first_el
         END FUNCTION cp_sll_3d_r_get_first_el
! **************************************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_get_first_el(sll) RESULT(res)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous :: res

            res =>sll%first_el
         END FUNCTION cp_sll_3d_c_get_first_el
! **************************************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_get_first_el(sll) RESULT(res)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:),pointer,contiguous :: res

            res =>sll%first_el
         END FUNCTION cp_sll_1d_r_get_first_el
! **************************************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_get_first_el(sll) RESULT(res)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:),pointer,contiguous :: res

            res =>sll%first_el
         END FUNCTION cp_sll_1d_c_get_first_el
! **************************************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_rs_get_first_el(sll) RESULT(res)
            TYPE(cp_sll_rs_type), POINTER :: sll
            type(realspace_grid_type),pointer :: res

            res =>sll%first_el
         END FUNCTION cp_sll_rs_get_first_el

! **************************************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_get_rest(sll, iter) RESULT(res)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            INTEGER, OPTIONAL :: iter

            TYPE(cp_sll_3d_r_type), POINTER :: res

            INTEGER :: i

            IF (.NOT. ASSOCIATED(sll)) THEN
               NULLIFY (res)
            ELSE
               IF (PRESENT(iter)) THEN
                  res => sll
                  DO i = 1, iter
                     IF (ASSOCIATED(res%rest)) THEN
                        res => res%rest
                     ELSE
                        CPABORT("tried to go past end")
                     END IF
                  END DO
                  IF (iter == -1) THEN
                     DO
                        IF (.NOT. ASSOCIATED(res%rest)) EXIT
                        res => res%rest
                     END DO
                  END IF
               ELSE
                  res => sll%rest ! make the common case fast...
               END IF
            END IF
         END FUNCTION cp_sll_3d_r_get_rest
! **************************************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_get_rest(sll, iter) RESULT(res)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            INTEGER, OPTIONAL :: iter

            TYPE(cp_sll_3d_c_type), POINTER :: res

            INTEGER :: i

            IF (.NOT. ASSOCIATED(sll)) THEN
               NULLIFY (res)
            ELSE
               IF (PRESENT(iter)) THEN
                  res => sll
                  DO i = 1, iter
                     IF (ASSOCIATED(res%rest)) THEN
                        res => res%rest
                     ELSE
                        CPABORT("tried to go past end")
                     END IF
                  END DO
                  IF (iter == -1) THEN
                     DO
                        IF (.NOT. ASSOCIATED(res%rest)) EXIT
                        res => res%rest
                     END DO
                  END IF
               ELSE
                  res => sll%rest ! make the common case fast...
               END IF
            END IF
         END FUNCTION cp_sll_3d_c_get_rest
! **************************************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_get_rest(sll, iter) RESULT(res)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            INTEGER, OPTIONAL :: iter

            TYPE(cp_sll_1d_r_type), POINTER :: res

            INTEGER :: i

            IF (.NOT. ASSOCIATED(sll)) THEN
               NULLIFY (res)
            ELSE
               IF (PRESENT(iter)) THEN
                  res => sll
                  DO i = 1, iter
                     IF (ASSOCIATED(res%rest)) THEN
                        res => res%rest
                     ELSE
                        CPABORT("tried to go past end")
                     END IF
                  END DO
                  IF (iter == -1) THEN
                     DO
                        IF (.NOT. ASSOCIATED(res%rest)) EXIT
                        res => res%rest
                     END DO
                  END IF
               ELSE
                  res => sll%rest ! make the common case fast...
               END IF
            END IF
         END FUNCTION cp_sll_1d_r_get_rest
! **************************************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_get_rest(sll, iter) RESULT(res)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            INTEGER, OPTIONAL :: iter

            TYPE(cp_sll_1d_c_type), POINTER :: res

            INTEGER :: i

            IF (.NOT. ASSOCIATED(sll)) THEN
               NULLIFY (res)
            ELSE
               IF (PRESENT(iter)) THEN
                  res => sll
                  DO i = 1, iter
                     IF (ASSOCIATED(res%rest)) THEN
                        res => res%rest
                     ELSE
                        CPABORT("tried to go past end")
                     END IF
                  END DO
                  IF (iter == -1) THEN
                     DO
                        IF (.NOT. ASSOCIATED(res%rest)) EXIT
                        res => res%rest
                     END DO
                  END IF
               ELSE
                  res => sll%rest ! make the common case fast...
               END IF
            END IF
         END FUNCTION cp_sll_1d_c_get_rest
! **************************************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! **************************************************************************************************
         FUNCTION cp_sll_rs_get_rest(sll, iter) RESULT(res)
            TYPE(cp_sll_rs_type), POINTER :: sll
            INTEGER, OPTIONAL :: iter

            TYPE(cp_sll_rs_type), POINTER :: res

            INTEGER :: i

            IF (.NOT. ASSOCIATED(sll)) THEN
               NULLIFY (res)
            ELSE
               IF (PRESENT(iter)) THEN
                  res => sll
                  DO i = 1, iter
                     IF (ASSOCIATED(res%rest)) THEN
                        res => res%rest
                     ELSE
                        CPABORT("tried to go past end")
                     END IF
                  END DO
                  IF (iter == -1) THEN
                     DO
                        IF (.NOT. ASSOCIATED(res%rest)) EXIT
                        res => res%rest
                     END DO
                  END IF
               ELSE
                  res => sll%rest ! make the common case fast...
               END IF
            END IF
         END FUNCTION cp_sll_rs_get_rest

! **************************************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_get_length(sll) RESULT(res)
            TYPE(cp_sll_3d_r_type), pointer :: sll
            INTEGER ::res

            TYPE(cp_sll_3d_r_type), POINTER :: iterator

            res = 0
            iterator => sll
            DO
               IF (ASSOCIATED(iterator)) THEN
                  res = res + 1
                  iterator => iterator%rest
               ELSE
                  EXIT
               END IF
            END DO
         END FUNCTION cp_sll_3d_r_get_length
! **************************************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_get_length(sll) RESULT(res)
            TYPE(cp_sll_3d_c_type), pointer :: sll
            INTEGER ::res

            TYPE(cp_sll_3d_c_type), POINTER :: iterator

            res = 0
            iterator => sll
            DO
               IF (ASSOCIATED(iterator)) THEN
                  res = res + 1
                  iterator => iterator%rest
               ELSE
                  EXIT
               END IF
            END DO
         END FUNCTION cp_sll_3d_c_get_length
! **************************************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_get_length(sll) RESULT(res)
            TYPE(cp_sll_1d_r_type), pointer :: sll
            INTEGER ::res

            TYPE(cp_sll_1d_r_type), POINTER :: iterator

            res = 0
            iterator => sll
            DO
               IF (ASSOCIATED(iterator)) THEN
                  res = res + 1
                  iterator => iterator%rest
               ELSE
                  EXIT
               END IF
            END DO
         END FUNCTION cp_sll_1d_r_get_length
! **************************************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_get_length(sll) RESULT(res)
            TYPE(cp_sll_1d_c_type), pointer :: sll
            INTEGER ::res

            TYPE(cp_sll_1d_c_type), POINTER :: iterator

            res = 0
            iterator => sll
            DO
               IF (ASSOCIATED(iterator)) THEN
                  res = res + 1
                  iterator => iterator%rest
               ELSE
                  EXIT
               END IF
            END DO
         END FUNCTION cp_sll_1d_c_get_length
! **************************************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! **************************************************************************************************
         FUNCTION cp_sll_rs_get_length(sll) RESULT(res)
            TYPE(cp_sll_rs_type), pointer :: sll
            INTEGER ::res

            TYPE(cp_sll_rs_type), POINTER :: iterator

            res = 0
            iterator => sll
            DO
               IF (ASSOCIATED(iterator)) THEN
                  res = res + 1
                  iterator => iterator%rest
               ELSE
                  EXIT
               END IF
            END DO
         END FUNCTION cp_sll_rs_get_length

! **************************************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_get_el_at(sll, index) RESULT(res)
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous :: res
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_3d_r_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_3d_r_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_3d_r_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            res =>pos%first_el
         END FUNCTION cp_sll_3d_r_get_el_at
! **************************************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_get_el_at(sll, index) RESULT(res)
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous :: res
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_3d_c_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_3d_c_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_3d_c_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            res =>pos%first_el
         END FUNCTION cp_sll_3d_c_get_el_at
! **************************************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_get_el_at(sll, index) RESULT(res)
            REAL(kind=dp),dimension(:),pointer,contiguous :: res
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_1d_r_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_1d_r_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_1d_r_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            res =>pos%first_el
         END FUNCTION cp_sll_1d_r_get_el_at
! **************************************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_get_el_at(sll, index) RESULT(res)
            COMPLEX(kind=dp),dimension(:),pointer,contiguous :: res
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_1d_c_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_1d_c_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_1d_c_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            res =>pos%first_el
         END FUNCTION cp_sll_1d_c_get_el_at
! **************************************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         FUNCTION cp_sll_rs_get_el_at(sll, index) RESULT(res)
            type(realspace_grid_type),pointer :: res
            TYPE(cp_sll_rs_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_rs_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_rs_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_rs_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            res =>pos%first_el
         END FUNCTION cp_sll_rs_get_el_at

! **************************************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_set_el_at(sll, index, value)
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous :: value
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_3d_r_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_3d_r_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_3d_r_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            pos%first_el =>value
         END SUBROUTINE cp_sll_3d_r_set_el_at
! **************************************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_set_el_at(sll, index, value)
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous :: value
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_3d_c_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_3d_c_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_3d_c_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            pos%first_el =>value
         END SUBROUTINE cp_sll_3d_c_set_el_at
! **************************************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_set_el_at(sll, index, value)
            REAL(kind=dp),dimension(:),pointer,contiguous :: value
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_1d_r_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_1d_r_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_1d_r_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            pos%first_el =>value
         END SUBROUTINE cp_sll_1d_r_set_el_at
! **************************************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_set_el_at(sll, index, value)
            COMPLEX(kind=dp),dimension(:),pointer,contiguous :: value
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_1d_c_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_1d_c_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_1d_c_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            pos%first_el =>value
         END SUBROUTINE cp_sll_1d_c_set_el_at
! **************************************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_set_el_at(sll, index, value)
            type(realspace_grid_type),pointer :: value
            TYPE(cp_sll_rs_type), POINTER :: sll
            INTEGER, INTENT(in) :: index

            TYPE(cp_sll_rs_type), POINTER :: pos

            IF (index == -1) THEN
               pos => cp_sll_rs_get_rest(sll, iter=-1)
            ELSE
               pos => cp_sll_rs_get_rest(sll, iter=index - 1)
            END IF
            CPASSERT(ASSOCIATED(pos))

            pos%first_el =>value
         END SUBROUTINE cp_sll_rs_set_el_at

! * iteration *

! **************************************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_next(iterator, el_att) RESULT(res)
            TYPE(cp_sll_3d_r_type), POINTER :: iterator
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous, OPTIONAL :: el_att
            LOGICAL :: res

            IF (ASSOCIATED(iterator)) THEN
               res = .TRUE.
               if (present(el_att)) el_att =>iterator%first_el
               iterator => iterator%rest
            ELSE
               res = .FALSE.
            END IF
         END FUNCTION cp_sll_3d_r_next
! **************************************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_next(iterator, el_att) RESULT(res)
            TYPE(cp_sll_3d_c_type), POINTER :: iterator
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous, OPTIONAL :: el_att
            LOGICAL :: res

            IF (ASSOCIATED(iterator)) THEN
               res = .TRUE.
               if (present(el_att)) el_att =>iterator%first_el
               iterator => iterator%rest
            ELSE
               res = .FALSE.
            END IF
         END FUNCTION cp_sll_3d_c_next
! **************************************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_next(iterator, el_att) RESULT(res)
            TYPE(cp_sll_1d_r_type), POINTER :: iterator
            REAL(kind=dp),dimension(:),pointer,contiguous, OPTIONAL :: el_att
            LOGICAL :: res

            IF (ASSOCIATED(iterator)) THEN
               res = .TRUE.
               if (present(el_att)) el_att =>iterator%first_el
               iterator => iterator%rest
            ELSE
               res = .FALSE.
            END IF
         END FUNCTION cp_sll_1d_r_next
! **************************************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_next(iterator, el_att) RESULT(res)
            TYPE(cp_sll_1d_c_type), POINTER :: iterator
            COMPLEX(kind=dp),dimension(:),pointer,contiguous, OPTIONAL :: el_att
            LOGICAL :: res

            IF (ASSOCIATED(iterator)) THEN
               res = .TRUE.
               if (present(el_att)) el_att =>iterator%first_el
               iterator => iterator%rest
            ELSE
               res = .FALSE.
            END IF
         END FUNCTION cp_sll_1d_c_next
! **************************************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \return ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_rs_next(iterator, el_att) RESULT(res)
            TYPE(cp_sll_rs_type), POINTER :: iterator
            type(realspace_grid_type),pointer, OPTIONAL :: el_att
            LOGICAL :: res

            IF (ASSOCIATED(iterator)) THEN
               res = .TRUE.
               if (present(el_att)) el_att =>iterator%first_el
               iterator => iterator%rest
            ELSE
               res = .FALSE.
            END IF
         END FUNCTION cp_sll_rs_next

! ============ structure modifications ============

! **************************************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_insert_el(sll, el)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous:: el

            TYPE(cp_sll_3d_r_type), POINTER :: newSlot

            NULLIFY (newSlot)

            CALL cp_sll_3d_r_create(newSlot, first_el=el, &
                                                rest=sll)
            sll => newSlot
         END SUBROUTINE cp_sll_3d_r_insert_el
! **************************************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_insert_el(sll, el)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous:: el

            TYPE(cp_sll_3d_c_type), POINTER :: newSlot

            NULLIFY (newSlot)

            CALL cp_sll_3d_c_create(newSlot, first_el=el, &
                                                rest=sll)
            sll => newSlot
         END SUBROUTINE cp_sll_3d_c_insert_el
! **************************************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_insert_el(sll, el)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            REAL(kind=dp),dimension(:),pointer,contiguous:: el

            TYPE(cp_sll_1d_r_type), POINTER :: newSlot

            NULLIFY (newSlot)

            CALL cp_sll_1d_r_create(newSlot, first_el=el, &
                                                rest=sll)
            sll => newSlot
         END SUBROUTINE cp_sll_1d_r_insert_el
! **************************************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_insert_el(sll, el)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            COMPLEX(kind=dp),dimension(:),pointer,contiguous:: el

            TYPE(cp_sll_1d_c_type), POINTER :: newSlot

            NULLIFY (newSlot)

            CALL cp_sll_1d_c_create(newSlot, first_el=el, &
                                                rest=sll)
            sll => newSlot
         END SUBROUTINE cp_sll_1d_c_insert_el
! **************************************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_insert_el(sll, el)
            TYPE(cp_sll_rs_type), POINTER :: sll
            type(realspace_grid_type),pointer:: el

            TYPE(cp_sll_rs_type), POINTER :: newSlot

            NULLIFY (newSlot)

            CALL cp_sll_rs_create(newSlot, first_el=el, &
                                                rest=sll)
            sll => newSlot
         END SUBROUTINE cp_sll_rs_insert_el

! **************************************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_rm_first_el(sll)
            TYPE(cp_sll_3d_r_type), POINTER :: sll

            TYPE(cp_sll_3d_r_type), POINTER :: node_to_rm
            node_to_rm => sll

            IF (ASSOCIATED(sll)) THEN
               sll => sll%rest
               CALL cp_sll_3d_r_dealloc_node(node_to_rm)
            ELSE
               CPABORT("tried to remove first el of an empty list")
            END IF
         END SUBROUTINE cp_sll_3d_r_rm_first_el
! **************************************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_rm_first_el(sll)
            TYPE(cp_sll_3d_c_type), POINTER :: sll

            TYPE(cp_sll_3d_c_type), POINTER :: node_to_rm
            node_to_rm => sll

            IF (ASSOCIATED(sll)) THEN
               sll => sll%rest
               CALL cp_sll_3d_c_dealloc_node(node_to_rm)
            ELSE
               CPABORT("tried to remove first el of an empty list")
            END IF
         END SUBROUTINE cp_sll_3d_c_rm_first_el
! **************************************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_rm_first_el(sll)
            TYPE(cp_sll_1d_r_type), POINTER :: sll

            TYPE(cp_sll_1d_r_type), POINTER :: node_to_rm
            node_to_rm => sll

            IF (ASSOCIATED(sll)) THEN
               sll => sll%rest
               CALL cp_sll_1d_r_dealloc_node(node_to_rm)
            ELSE
               CPABORT("tried to remove first el of an empty list")
            END IF
         END SUBROUTINE cp_sll_1d_r_rm_first_el
! **************************************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_rm_first_el(sll)
            TYPE(cp_sll_1d_c_type), POINTER :: sll

            TYPE(cp_sll_1d_c_type), POINTER :: node_to_rm
            node_to_rm => sll

            IF (ASSOCIATED(sll)) THEN
               sll => sll%rest
               CALL cp_sll_1d_c_dealloc_node(node_to_rm)
            ELSE
               CPABORT("tried to remove first el of an empty list")
            END IF
         END SUBROUTINE cp_sll_1d_c_rm_first_el
! **************************************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_rm_first_el(sll)
            TYPE(cp_sll_rs_type), POINTER :: sll

            TYPE(cp_sll_rs_type), POINTER :: node_to_rm
            node_to_rm => sll

            IF (ASSOCIATED(sll)) THEN
               sll => sll%rest
               CALL cp_sll_rs_dealloc_node(node_to_rm)
            ELSE
               CPABORT("tried to remove first el of an empty list")
            END IF
         END SUBROUTINE cp_sll_rs_rm_first_el

! **************************************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_insert_el_at(sll, el, index)
            REAL(kind=dp),dimension(:,:,:),pointer,contiguous :: el
            INTEGER, INTENT(in) :: index
            TYPE(cp_sll_3d_r_type), POINTER :: sll

            TYPE(cp_sll_3d_r_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_3d_r_insert_el(sll, el)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_3d_r_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_3d_r_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_3d_r_insert_el(pos%rest, el)
            END IF
         END SUBROUTINE cp_sll_3d_r_insert_el_at
! **************************************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_insert_el_at(sll, el, index)
            COMPLEX(kind=dp),dimension(:,:,:),pointer,contiguous :: el
            INTEGER, INTENT(in) :: index
            TYPE(cp_sll_3d_c_type), POINTER :: sll

            TYPE(cp_sll_3d_c_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_3d_c_insert_el(sll, el)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_3d_c_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_3d_c_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_3d_c_insert_el(pos%rest, el)
            END IF
         END SUBROUTINE cp_sll_3d_c_insert_el_at
! **************************************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_insert_el_at(sll, el, index)
            REAL(kind=dp),dimension(:),pointer,contiguous :: el
            INTEGER, INTENT(in) :: index
            TYPE(cp_sll_1d_r_type), POINTER :: sll

            TYPE(cp_sll_1d_r_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_1d_r_insert_el(sll, el)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_1d_r_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_1d_r_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_1d_r_insert_el(pos%rest, el)
            END IF
         END SUBROUTINE cp_sll_1d_r_insert_el_at
! **************************************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_insert_el_at(sll, el, index)
            COMPLEX(kind=dp),dimension(:),pointer,contiguous :: el
            INTEGER, INTENT(in) :: index
            TYPE(cp_sll_1d_c_type), POINTER :: sll

            TYPE(cp_sll_1d_c_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_1d_c_insert_el(sll, el)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_1d_c_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_1d_c_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_1d_c_insert_el(pos%rest, el)
            END IF
         END SUBROUTINE cp_sll_1d_c_insert_el_at
! **************************************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_insert_el_at(sll, el, index)
            type(realspace_grid_type),pointer :: el
            INTEGER, INTENT(in) :: index
            TYPE(cp_sll_rs_type), POINTER :: sll

            TYPE(cp_sll_rs_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_rs_insert_el(sll, el)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_rs_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_rs_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_rs_insert_el(pos%rest, el)
            END IF
         END SUBROUTINE cp_sll_rs_insert_el_at

! **************************************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_rm_el_at(sll, index)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            INTEGER, INTENT(in)::index

            TYPE(cp_sll_3d_r_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_3d_r_rm_first_el(sll)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_3d_r_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_3d_r_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_3d_r_rm_first_el(pos%rest)
            END IF
         END SUBROUTINE cp_sll_3d_r_rm_el_at
! **************************************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_rm_el_at(sll, index)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            INTEGER, INTENT(in)::index

            TYPE(cp_sll_3d_c_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_3d_c_rm_first_el(sll)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_3d_c_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_3d_c_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_3d_c_rm_first_el(pos%rest)
            END IF
         END SUBROUTINE cp_sll_3d_c_rm_el_at
! **************************************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_rm_el_at(sll, index)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            INTEGER, INTENT(in)::index

            TYPE(cp_sll_1d_r_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_1d_r_rm_first_el(sll)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_1d_r_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_1d_r_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_1d_r_rm_first_el(pos%rest)
            END IF
         END SUBROUTINE cp_sll_1d_r_rm_el_at
! **************************************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_rm_el_at(sll, index)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            INTEGER, INTENT(in)::index

            TYPE(cp_sll_1d_c_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_1d_c_rm_first_el(sll)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_1d_c_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_1d_c_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_1d_c_rm_first_el(pos%rest)
            END IF
         END SUBROUTINE cp_sll_1d_c_rm_el_at
! **************************************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_rm_el_at(sll, index)
            TYPE(cp_sll_rs_type), POINTER :: sll
            INTEGER, INTENT(in)::index

            TYPE(cp_sll_rs_type), POINTER :: pos

            IF (index == 1) THEN
               CALL cp_sll_rs_rm_first_el(sll)
            ELSE
               IF (index == -1) THEN
                  pos => cp_sll_rs_get_rest(sll, iter=-1)
               ELSE
                  pos => cp_sll_rs_get_rest(sll, iter=index - 2)
               END IF
               CPASSERT(ASSOCIATED(pos))
               CALL cp_sll_rs_rm_first_el(pos%rest)
            END IF
         END SUBROUTINE cp_sll_rs_rm_el_at

! **************************************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_r_rm_all_el(sll)
            TYPE(cp_sll_3d_r_type), POINTER :: sll

            TYPE(cp_sll_3d_r_type), POINTER :: next_node, actual_node

            actual_node => sll
            DO
               IF (.NOT. ASSOCIATED(actual_node)) EXIT
               next_node => actual_node%rest
               CALL cp_sll_3d_r_dealloc_node(actual_node)
               actual_node => next_node
            END DO
            NULLIFY (sll)
         END SUBROUTINE cp_sll_3d_r_rm_all_el
! **************************************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_3d_c_rm_all_el(sll)
            TYPE(cp_sll_3d_c_type), POINTER :: sll

            TYPE(cp_sll_3d_c_type), POINTER :: next_node, actual_node

            actual_node => sll
            DO
               IF (.NOT. ASSOCIATED(actual_node)) EXIT
               next_node => actual_node%rest
               CALL cp_sll_3d_c_dealloc_node(actual_node)
               actual_node => next_node
            END DO
            NULLIFY (sll)
         END SUBROUTINE cp_sll_3d_c_rm_all_el
! **************************************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_r_rm_all_el(sll)
            TYPE(cp_sll_1d_r_type), POINTER :: sll

            TYPE(cp_sll_1d_r_type), POINTER :: next_node, actual_node

            actual_node => sll
            DO
               IF (.NOT. ASSOCIATED(actual_node)) EXIT
               next_node => actual_node%rest
               CALL cp_sll_1d_r_dealloc_node(actual_node)
               actual_node => next_node
            END DO
            NULLIFY (sll)
         END SUBROUTINE cp_sll_1d_r_rm_all_el
! **************************************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_1d_c_rm_all_el(sll)
            TYPE(cp_sll_1d_c_type), POINTER :: sll

            TYPE(cp_sll_1d_c_type), POINTER :: next_node, actual_node

            actual_node => sll
            DO
               IF (.NOT. ASSOCIATED(actual_node)) EXIT
               next_node => actual_node%rest
               CALL cp_sll_1d_c_dealloc_node(actual_node)
               actual_node => next_node
            END DO
            NULLIFY (sll)
         END SUBROUTINE cp_sll_1d_c_rm_all_el
! **************************************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! **************************************************************************************************
         SUBROUTINE cp_sll_rs_rm_all_el(sll)
            TYPE(cp_sll_rs_type), POINTER :: sll

            TYPE(cp_sll_rs_type), POINTER :: next_node, actual_node

            actual_node => sll
            DO
               IF (.NOT. ASSOCIATED(actual_node)) EXIT
               next_node => actual_node%rest
               CALL cp_sll_rs_dealloc_node(actual_node)
               actual_node => next_node
            END DO
            NULLIFY (sll)
         END SUBROUTINE cp_sll_rs_rm_all_el

! **************************************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to transform in array
!> \return ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_r_to_array(sll) RESULT(res)
            TYPE(cp_sll_3d_r_type), POINTER :: sll
            type(cp_3d_r_cp_type), DIMENSION(:), POINTER :: res

            INTEGER :: len, i
            LOGICAL :: ok
            TYPE(cp_sll_3d_r_type), POINTER :: iter

            len = cp_sll_3d_r_get_length(sll)
            ALLOCATE (res(len))
            iter => sll
            DO i = 1, len
               res(i) %array=>iter%first_el
               ok = cp_sll_3d_r_next(iter)
               CPASSERT(ok .OR. i == len)
            END DO
         END FUNCTION cp_sll_3d_r_to_array
! **************************************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to transform in array
!> \return ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_3d_c_to_array(sll) RESULT(res)
            TYPE(cp_sll_3d_c_type), POINTER :: sll
            type(cp_3d_c_cp_type), DIMENSION(:), POINTER :: res

            INTEGER :: len, i
            LOGICAL :: ok
            TYPE(cp_sll_3d_c_type), POINTER :: iter

            len = cp_sll_3d_c_get_length(sll)
            ALLOCATE (res(len))
            iter => sll
            DO i = 1, len
               res(i) %array=>iter%first_el
               ok = cp_sll_3d_c_next(iter)
               CPASSERT(ok .OR. i == len)
            END DO
         END FUNCTION cp_sll_3d_c_to_array
! **************************************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to transform in array
!> \return ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_r_to_array(sll) RESULT(res)
            TYPE(cp_sll_1d_r_type), POINTER :: sll
            type(cp_1d_r_cp_type), DIMENSION(:), POINTER :: res

            INTEGER :: len, i
            LOGICAL :: ok
            TYPE(cp_sll_1d_r_type), POINTER :: iter

            len = cp_sll_1d_r_get_length(sll)
            ALLOCATE (res(len))
            iter => sll
            DO i = 1, len
               res(i) %array=>iter%first_el
               ok = cp_sll_1d_r_next(iter)
               CPASSERT(ok .OR. i == len)
            END DO
         END FUNCTION cp_sll_1d_r_to_array
! **************************************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to transform in array
!> \return ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_1d_c_to_array(sll) RESULT(res)
            TYPE(cp_sll_1d_c_type), POINTER :: sll
            type(cp_1d_c_cp_type), DIMENSION(:), POINTER :: res

            INTEGER :: len, i
            LOGICAL :: ok
            TYPE(cp_sll_1d_c_type), POINTER :: iter

            len = cp_sll_1d_c_get_length(sll)
            ALLOCATE (res(len))
            iter => sll
            DO i = 1, len
               res(i) %array=>iter%first_el
               ok = cp_sll_1d_c_next(iter)
               CPASSERT(ok .OR. i == len)
            END DO
         END FUNCTION cp_sll_1d_c_to_array
! **************************************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to transform in array
!> \return ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         FUNCTION cp_sll_rs_to_array(sll) RESULT(res)
            TYPE(cp_sll_rs_type), POINTER :: sll
            type(realspace_grid_p_type), DIMENSION(:), POINTER :: res

            INTEGER :: len, i
            LOGICAL :: ok
            TYPE(cp_sll_rs_type), POINTER :: iter

            len = cp_sll_rs_get_length(sll)
            ALLOCATE (res(len))
            iter => sll
            DO i = 1, len
               res(i) %rs_grid=>iter%first_el
               ok = cp_sll_rs_next(iter)
               CPASSERT(ok .OR. i == len)
            END DO
         END FUNCTION cp_sll_rs_to_array
END MODULE cp_linked_list_pw
