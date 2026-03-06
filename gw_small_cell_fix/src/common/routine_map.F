!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE routine_map
   USE kinds, ONLY: default_string_length, int_4, int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE


! **************************************************************************************************
!> \brief A hash map (also known as hashtable or dictionary).
!>        Internally the hash map uses an array to holds its data.
!>        If this array reaches a load-factor of 75%, a new array with twice the
!>        size will be allocated and the items are then copied over.
!>        This ensures that the dictionary will perform operations in O(1).
!> \par History
!>      12.2012 first version [Ole Schuett]
!>      08.2019 refactored for Fypp [Ole Schuett]
!> \author Ole Schuett
! ***************************************************************************************************

      PUBLIC :: routine_map_init
      PUBLIC :: routine_map_items
      PUBLIC :: routine_map_haskey
      PUBLIC :: routine_map_set
      PUBLIC :: routine_map_get
      PUBLIC :: routine_map_size
      PUBLIC :: routine_map_destroy
      PUBLIC :: routine_map_type
      PUBLIC :: routine_map_item_type

!this is an internal type
!Calculating hashes might be expensive, therefore they are stored
!for use during change_capacity().
      TYPE private_item_type
         PRIVATE
         CHARACTER(LEN=default_string_length)                                :: key  = ""
         INTEGER(kind=int_4)                              :: value  = 0_int_4
         INTEGER(KIND=int_8)                         :: hash = 0_int_8
         TYPE(private_item_type), POINTER            :: next => Null()
      END TYPE private_item_type

!this is an internal type
      TYPE private_item_p_type
         PRIVATE
         TYPE(private_item_type), POINTER           :: p => Null()
      END TYPE private_item_p_type

! this is the public type, which holds a hash map instance
      TYPE routine_map_type
         PRIVATE
         TYPE(private_item_p_type), DIMENSION(:), POINTER      :: buckets => Null()
         INTEGER                                               :: size = -1
      END TYPE routine_map_type

! this is a public type, its returned by  routine_map_items()
      TYPE routine_map_item_type
         CHARACTER(LEN=default_string_length)        :: key  = ""
         INTEGER(kind=int_4)      :: value  = 0_int_4
      END TYPE routine_map_item_type

      CONTAINS


! **************************************************************************************************
!> \brief Allocates the internal data-structures of the given hash map.
!> \param hash_map ...
!> \param initial_capacity The initial size of the internal array (default=11).
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE routine_map_init(hash_map, initial_capacity)
         TYPE(routine_map_type), INTENT(inout)                 :: hash_map
         INTEGER, INTENT(in), OPTIONAL                             :: initial_capacity

         INTEGER :: initial_capacity_

         IF (PRESENT(initial_capacity)) THEN
            initial_capacity_ = initial_capacity
         ELSE
            initial_capacity_ = 11
         END IF

         IF (initial_capacity_ < 1) &
            CPABORT("initial_capacity < 1")

         IF (ASSOCIATED(hash_map%buckets)) &
            CPABORT("hash map is already initialized.")

         ALLOCATE (hash_map%buckets(initial_capacity_))
         hash_map%size = 0

      END SUBROUTINE routine_map_init

! **************************************************************************************************
!> \brief Test if the given hash map has been initialized.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION routine_map_isready(hash_map) RESULT(res)
         TYPE(routine_map_type), INTENT(inout)             :: hash_map
         LOGICAL                                               :: res
         res = ASSOCIATED(hash_map%buckets)
      END FUNCTION routine_map_isready

! **************************************************************************************************
!> \brief Deallocated the internal data-structures if the given hash map.
!>        Caution: If the stored keys or values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param hash_map ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE routine_map_destroy(hash_map)
         TYPE(routine_map_type), INTENT(inout)  :: hash_map
         TYPE(private_item_type), POINTER           :: item, prev_item
         INTEGER :: i

         CPASSERT(ASSOCIATED(hash_map%buckets))

         DO i = 1, size(hash_map%buckets)
            item => hash_map%buckets(i)%p
            DO WHILE (ASSOCIATED(item))
               prev_item => item
               item => item%next
               DEALLOCATE (prev_item)
            END DO
         END DO

         DEALLOCATE (hash_map%buckets)
         hash_map%size = -1
      END SUBROUTINE routine_map_destroy

! **************************************************************************************************
!> \brief Stores, and possibly overwrites, a given value under a given key.
!> \param hash_map ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE routine_map_set(hash_map, key, value)
         TYPE(routine_map_type), INTENT(inout)             :: hash_map
         CHARACTER(LEN=default_string_length), INTENT(in)                            :: key
         INTEGER(kind=int_4), INTENT(in)                            :: value
         INTEGER(KIND=int_8)                                   :: hash
         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = routine_map_hash_function(key)
         CALL routine_map_set_hashed(hash_map, key, value, hash)
      END SUBROUTINE routine_map_set

! **************************************************************************************************
!> \brief Common code used internally by routine_map_set() and routine_map_change_capacity().
!> \param hash_map ...
!> \param key ...
!> \param value ...
!> \param hash ...
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE routine_map_set_hashed(hash_map, key, value, hash)
         TYPE(routine_map_type), INTENT(inout)               :: hash_map
         CHARACTER(LEN=default_string_length), intent(in)                              :: key
         INTEGER(kind=int_4), intent(in)                              :: value
         INTEGER(KIND=int_8), intent(in)                         :: hash
         TYPE(private_item_type), POINTER                        :: item, new_item
         INTEGER(KIND=int_8)                                     :: idx

         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         ! if already in hash map just update its value
         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (routine_map_keys_equal(item%key, key)) THEN
                  item%value =value
                  RETURN
               END IF
            END IF
            item => item%next
         END DO

         ! check load-factor
         IF (4*hash_map%size > 3*size(hash_map%buckets)) THEN ! load-factor > 75%
            call routine_map_change_capacity(hash_map, 2*size(hash_map%buckets)) !double capacity
            idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1
         END IF

         ! create a new item
         allocate (new_item)
         new_item%hash = hash
         new_item%key =key
         new_item%value =value
         new_item%next => hash_map%buckets(idx)%p
         hash_map%buckets(idx)%p => new_item
         hash_map%size = hash_map%size + 1

      END SUBROUTINE routine_map_set_hashed

! **************************************************************************************************
!> \brief Internal routine for changing the hash map's capacity.
!> \param hash_map ...
!> \param new_capacity ...
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE routine_map_change_capacity(hash_map, new_capacity)
         TYPE(routine_map_type), INTENT(inout)         :: hash_map
         INTEGER, INTENT(in)                               :: new_capacity
         INTEGER                                           :: i, old_size, new_cap
         TYPE(private_item_type), POINTER                  :: item, prev_item
         TYPE(private_item_p_type), DIMENSION(:), POINTER  :: old_buckets
         new_cap = new_capacity
         ! pre checks
         IF (new_cap > HUGE(i)) THEN
            IF (size(hash_map%buckets) == HUGE(i)) RETURN ! reached maximum - stay there.
            new_cap = HUGE(i) ! grow as far as possible
         END IF
         CPASSERT(new_cap >= 1)
         CPASSERT(4*hash_map%size < 3*new_cap)

         old_size = hash_map%size
         old_buckets => hash_map%buckets
         ALLOCATE (hash_map%buckets(new_capacity))
         hash_map%size = 0
         DO i = 1, size(old_buckets)
            item => old_buckets(i)%p
            DO WHILE (ASSOCIATED(item))
               CALL routine_map_set_hashed(hash_map, item%key, item%value, item%hash)
               prev_item => item
               item => item%next
               DEALLOCATE (prev_item)
            END DO
         END DO

         DEALLOCATE (old_buckets)

         CPASSERT(old_size == hash_map%size)
      END SUBROUTINE routine_map_change_capacity

! **************************************************************************************************
!> \brief Gets a value for a given key from the hash map.
!>        If the key is not found the default_value will be returned.
!>        If the key is not found and default_value was not provided the program stops.
!> \param hash_map ...
!> \param key ...
!> \param default_value ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION routine_map_get(hash_map, key, default_value) RESULT(value)
         TYPE(routine_map_type), INTENT(in)              :: hash_map
         CHARACTER(LEN=default_string_length), INTENT(in)                            :: key
         INTEGER(kind=int_4), INTENT(in), OPTIONAL                :: default_value
         INTEGER(kind=int_4)                                      :: value
         TYPE(private_item_type), POINTER                    :: item
         INTEGER(KIND=int_8)                                 :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = routine_map_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (routine_map_keys_equal(item%key, key)) THEN
                  value =item%value
                  RETURN
               END IF
            END IF
            item => item%next
         END DO

         IF (PRESENT(default_value)) THEN
            value =default_value
            RETURN
         END IF

         CPABORT("Key not found.")
      END FUNCTION routine_map_get

! **************************************************************************************************
!> \brief Remove the value for a given key from the hash map.
!>        If the key is not found the program stops.
!> \param hash_map ...
!> \param key ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE routine_map_del(hash_map, key)
         TYPE(routine_map_type), INTENT(inout)               :: hash_map
         CHARACTER(LEN=default_string_length), INTENT(in)                                :: key
         TYPE(private_item_type), POINTER                        :: item, prev_item
         INTEGER(KIND=int_8)                                     :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = routine_map_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         prev_item => Null()
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (routine_map_keys_equal(item%key, key)) THEN
                  IF (ASSOCIATED(prev_item)) THEN
                     prev_item%next => item%next
                  ELSE
                     hash_map%buckets(idx)%p => item%next
                  END IF
                  DEALLOCATE (item)
                  hash_map%size = hash_map%size - 1
                  RETURN
               END IF
            END IF
            prev_item => item
            item => item%next
         END DO

         CPABORT("Key not found.")
      END SUBROUTINE routine_map_del

! **************************************************************************************************
!> \brief Returns the number of key/value-items currently stored in the hash map.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION routine_map_size(hash_map) RESULT(size)
         TYPE(routine_map_type), INTENT(IN)  :: hash_map
         INTEGER                                 :: size

         CPASSERT(ASSOCIATED(hash_map%buckets))
         size = hash_map%size
      END FUNCTION routine_map_size

! **************************************************************************************************
!> \brief Checks whether a given key is currently stored in the hash_map.
!> \param hash_map ...
!> \param key ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION routine_map_haskey(hash_map, key) RESULT(res)
         TYPE(routine_map_type), INTENT(IN)                :: hash_map
         CHARACTER(LEN=default_string_length), INTENT(IN)                              :: key
         LOGICAL                                               :: res
         TYPE(private_item_type), POINTER                      :: item
         INTEGER(KIND=int_8)                                   :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         res = .FALSE.
         IF (hash_map%size == 0) RETURN

         hash = routine_map_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (routine_map_keys_equal(item%key, key)) THEN
                  res = .TRUE.
                  return
               END IF
            END IF
            item => item%next
         END DO

      END FUNCTION routine_map_haskey

! **************************************************************************************************
!> \brief Returns a pointer to an array of all key/value-items stored in the hash map.
!>        Caution: The caller is responsible for deallocating targeted array after usage.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION routine_map_items(hash_map) RESULT(items)
         TYPE(routine_map_type), INTENT(IN)                 :: hash_map
         TYPE(routine_map_item_type), DIMENSION(:), POINTER :: items

         TYPE(private_item_type), POINTER  :: item
         INTEGER :: i, j

         CPASSERT(ASSOCIATED(hash_map%buckets))

         ALLOCATE (items(hash_map%size))
         j = 1
         DO i = 1, size(hash_map%buckets)
            item => hash_map%buckets(i)%p
            DO WHILE (ASSOCIATED(item))
               items(j)%key =item%key
               items(j)%value =item%value
               j = j + 1
               item => item%next
            END DO
         END DO

         CPASSERT(j == hash_map%size + 1)
      END FUNCTION routine_map_items

! **************************************************************************************************
!> \brief Copies all key/values-items from one hash map to another.
!>        Afterwards hash_map will contain all items from the from_hash_map and
!>        additionally all its previous items, which were not overwritten.
!>        The two hash maps have to be of the same type.
!> \param hash_map destination of items
!> \param from_hash_map source of items - will not be change
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE routine_map_update(hash_map, from_hash_map)
         TYPE(routine_map_type), INTENT(inout)              :: hash_map
         TYPE(routine_map_type), INTENT(in)                 :: from_hash_map
         TYPE(routine_map_item_type), DIMENSION(:), POINTER :: from_items
         INTEGER :: i

         CPASSERT(ASSOCIATED(hash_map%buckets))
         CPASSERT(ASSOCIATED(from_hash_map%buckets))

         from_items => routine_map_items(from_hash_map)
         DO i = 1, size(from_items)
            CALL routine_map_set(hash_map, from_items(i)%key, from_items(i)%value)
         END DO
         DEALLOCATE (from_items)
      END SUBROUTINE routine_map_update


! **************************************************************************************************
! This is joaat_hash from string_table.F
!
!> \brief generates the hash of a given string
!> \param key a string of any length
!> \return ...
!> \par History
!>       09.2006 created [Joost VandeVondele]
!>       12.2012 copied and adopted [ole]
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
! **************************************************************************************************
   PURE FUNCTION routine_map_hash_function(key) RESULT(hash)
      CHARACTER(LEN=*), INTENT(IN)                       :: key
      INTEGER(KIND=int_8)                                :: hash

      INTEGER(KIND=int_8), PARAMETER                     :: b32 = 2_int_8**32 - 1_int_8

      INTEGER                                            :: i

      hash = 0_int_8
      DO i = 1, LEN(key)
         hash = IAND(hash + ICHAR(key(i:i)), b32)
         hash = IAND(hash + IAND(ISHFT(hash, 10), b32), b32)
         hash = IAND(IEOR(hash, IAND(ISHFT(hash, -6), b32)), b32)
      END DO
      hash = IAND(hash + IAND(ISHFT(hash, 3), b32), b32)
      hash = IAND(IEOR(hash, IAND(ISHFT(hash, -11), b32)), b32)
      hash = IAND(hash + IAND(ISHFT(hash, 15), b32), b32)
   END FUNCTION routine_map_hash_function

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION routine_map_keys_equal(key1, key2) RESULT(res)
      CHARACTER(LEN=*), INTENT(IN)                      :: key1, key2
      LOGICAL                                            :: res

      res = (key1 == key2)
   END FUNCTION routine_map_keys_equal

END MODULE routine_map
