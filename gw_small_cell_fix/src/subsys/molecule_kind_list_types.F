!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a simple array based list of the given type
!> \par History
!>      08.2003 created [fawzi]
!>      01.2017 ported to Fypp [Ole Schuett]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE molecule_kind_list_types
   use molecule_kind_types, only: molecule_kind_type, deallocate_molecule_kind_set

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE




      LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

      !API
      PUBLIC :: molecule_kind_list_type, molecule_kind_list_p_type, &
                molecule_kind_list_create, molecule_kind_list_retain, &
                molecule_kind_list_release

! **************************************************************************************************
!> \brief represent a list of objects
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param n_el the number of elements in the list
!> \param owns_list if the els are owned by this structure, and
!>        should be deallocated by it
!> \param list the array of object, might be oversized,
!>        only the fist n_el have some meaning
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      TYPE molecule_kind_list_type
         INTEGER :: ref_count = 0, n_els = 0
         LOGICAL :: owns_els = .FALSE.
         TYPE(molecule_kind_type), DIMENSION(:), POINTER :: els => NULL()
      END TYPE molecule_kind_list_type

! **************************************************************************************************
!> \brief represents a pointer to a list
!> \param list the pointer to the list
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      TYPE molecule_kind_list_p_type
         TYPE(molecule_kind_list_type), POINTER :: list => NULL()
      END TYPE molecule_kind_list_p_type

      CONTAINS

! **************************************************************************************************
!> \brief creates a list
!> \param list the list to allocate and initialize
!> \param els_ptr the elements to store in the list (the array is only,
!>        referenced, not copied!)
!> \param owns_els if the list takes the ownership of els_ptr and
!>        will deallocate it (defaults to true)
!> \param n_els number of elements in the list (at least one els_ptr or
!>        n_els should be given)
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE molecule_kind_list_create(list, els_ptr, &
                                      owns_els, n_els)
         TYPE(molecule_kind_list_type), POINTER, OPTIONAL :: list
         TYPE(molecule_kind_type), DIMENSION(:), POINTER, OPTIONAL :: els_ptr
         LOGICAL, INTENT(in), OPTIONAL :: owns_els
         INTEGER, INTENT(in), OPTIONAL :: n_els

         CPASSERT(PRESENT(els_ptr) .OR. PRESENT(n_els))

         ALLOCATE (list)
         list%ref_count = 1
         list%owns_els = .TRUE.
         list%n_els = 0
         IF (PRESENT(owns_els)) list%owns_els = owns_els
         NULLIFY (list%els)
         IF (PRESENT(els_ptr)) THEN
            list%els => els_ptr
            IF (ASSOCIATED(els_ptr)) THEN
               list%n_els = SIZE(els_ptr)
            END IF
         END IF
         IF (PRESENT(n_els)) list%n_els = n_els
         IF (.NOT. ASSOCIATED(list%els)) THEN
            ALLOCATE (list%els(list%n_els))
            CPASSERT(list%owns_els)
         END IF
      END SUBROUTINE molecule_kind_list_create

! **************************************************************************************************
!> \brief retains a list (see doc/ReferenceCounting.html)
!> \param list the list to retain
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE molecule_kind_list_retain(list)
         TYPE(molecule_kind_list_type), POINTER :: list

         CPASSERT(ASSOCIATED(list))
         CPASSERT(list%ref_count > 0)
         list%ref_count = list%ref_count + 1
      END SUBROUTINE molecule_kind_list_retain

! **************************************************************************************************
!> \brief releases a  list (see doc/ReferenceCounting.html)
!> \param list the list to release
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE molecule_kind_list_release(list)
         TYPE(molecule_kind_list_type), POINTER :: list

         IF (ASSOCIATED(list)) THEN
            CPASSERT(list%ref_count > 0)
            list%ref_count = list%ref_count - 1
            IF (list%ref_count == 0) THEN
               IF (list%owns_els) THEN
                  IF (ASSOCIATED(list%els)) THEN
                     CALL deallocate_molecule_kind_set(list%els)
                  END IF
               END IF
               NULLIFY (list%els)
               DEALLOCATE (list)
            END IF
         END IF
         NULLIFY (list)
      END SUBROUTINE molecule_kind_list_release


END MODULE molecule_kind_list_types
