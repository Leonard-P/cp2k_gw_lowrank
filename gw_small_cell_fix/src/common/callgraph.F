!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE callgraph
   USE kinds, ONLY: int_4, int_8
   USE timings_base_type, ONLY: call_stat_type
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

      PUBLIC :: callgraph_init
      PUBLIC :: callgraph_items
      PUBLIC :: callgraph_haskey
      PUBLIC :: callgraph_set
      PUBLIC :: callgraph_get
      PUBLIC :: callgraph_size
      PUBLIC :: callgraph_destroy
      PUBLIC :: callgraph_type
      PUBLIC :: callgraph_item_type

!this is an internal type
!Calculating hashes might be expensive, therefore they are stored
!for use during change_capacity().
      TYPE private_item_type
         PRIVATE
         INTEGER(kind=int_4), DIMENSION(2)                                :: key  = 0_int_4
         TYPE(call_stat_type), POINTER                              :: value  => NULL()
         INTEGER(KIND=int_8)                         :: hash = 0_int_8
         TYPE(private_item_type), POINTER            :: next => Null()
      END TYPE private_item_type

!this is an internal type
      TYPE private_item_p_type
         PRIVATE
         TYPE(private_item_type), POINTER           :: p => Null()
      END TYPE private_item_p_type

! this is the public type, which holds a hash map instance
      TYPE callgraph_type
         PRIVATE
         TYPE(private_item_p_type), DIMENSION(:), POINTER      :: buckets => Null()
         INTEGER                                               :: size = -1
      END TYPE callgraph_type

! this is a public type, its returned by  callgraph_items()
      TYPE callgraph_item_type
         INTEGER(kind=int_4), DIMENSION(2)        :: key  = 0_int_4
         TYPE(call_stat_type), POINTER      :: value  => NULL()
      END TYPE callgraph_item_type

      CONTAINS


! **************************************************************************************************
!> \brief Allocates the internal data-structures of the given hash map.
!> \param hash_map ...
!> \param initial_capacity The initial size of the internal array (default=11).
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE callgraph_init(hash_map, initial_capacity)
         TYPE(callgraph_type), INTENT(inout)                 :: hash_map
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

      END SUBROUTINE callgraph_init

! **************************************************************************************************
!> \brief Test if the given hash map has been initialized.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION callgraph_isready(hash_map) RESULT(res)
         TYPE(callgraph_type), INTENT(inout)             :: hash_map
         LOGICAL                                               :: res
         res = ASSOCIATED(hash_map%buckets)
      END FUNCTION callgraph_isready

! **************************************************************************************************
!> \brief Deallocated the internal data-structures if the given hash map.
!>        Caution: If the stored keys or values are pointers, their targets will
!>                 not get deallocated by this routine.
!> \param hash_map ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE callgraph_destroy(hash_map)
         TYPE(callgraph_type), INTENT(inout)  :: hash_map
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
      END SUBROUTINE callgraph_destroy

! **************************************************************************************************
!> \brief Stores, and possibly overwrites, a given value under a given key.
!> \param hash_map ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE callgraph_set(hash_map, key, value)
         TYPE(callgraph_type), INTENT(inout)             :: hash_map
         INTEGER(kind=int_4), DIMENSION(2), INTENT(in)                            :: key
         TYPE(call_stat_type), POINTER, INTENT(in)                            :: value
         INTEGER(KIND=int_8)                                   :: hash
         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = callgraph_hash_function(key)
         CALL callgraph_set_hashed(hash_map, key, value, hash)
      END SUBROUTINE callgraph_set

! **************************************************************************************************
!> \brief Common code used internally by callgraph_set() and callgraph_change_capacity().
!> \param hash_map ...
!> \param key ...
!> \param value ...
!> \param hash ...
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE callgraph_set_hashed(hash_map, key, value, hash)
         TYPE(callgraph_type), INTENT(inout)               :: hash_map
         INTEGER(kind=int_4), DIMENSION(2), intent(in)                              :: key
         TYPE(call_stat_type), POINTER, intent(in)                              :: value
         INTEGER(KIND=int_8), intent(in)                         :: hash
         TYPE(private_item_type), POINTER                        :: item, new_item
         INTEGER(KIND=int_8)                                     :: idx

         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         ! if already in hash map just update its value
         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (callgraph_keys_equal(item%key, key)) THEN
                  item%value =>value
                  RETURN
               END IF
            END IF
            item => item%next
         END DO

         ! check load-factor
         IF (4*hash_map%size > 3*size(hash_map%buckets)) THEN ! load-factor > 75%
            call callgraph_change_capacity(hash_map, 2*size(hash_map%buckets)) !double capacity
            idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1
         END IF

         ! create a new item
         allocate (new_item)
         new_item%hash = hash
         new_item%key =key
         new_item%value =>value
         new_item%next => hash_map%buckets(idx)%p
         hash_map%buckets(idx)%p => new_item
         hash_map%size = hash_map%size + 1

      END SUBROUTINE callgraph_set_hashed

! **************************************************************************************************
!> \brief Internal routine for changing the hash map's capacity.
!> \param hash_map ...
!> \param new_capacity ...
!> \author Ole Schuett
! **************************************************************************************************
      RECURSIVE SUBROUTINE callgraph_change_capacity(hash_map, new_capacity)
         TYPE(callgraph_type), INTENT(inout)         :: hash_map
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
               CALL callgraph_set_hashed(hash_map, item%key, item%value, item%hash)
               prev_item => item
               item => item%next
               DEALLOCATE (prev_item)
            END DO
         END DO

         DEALLOCATE (old_buckets)

         CPASSERT(old_size == hash_map%size)
      END SUBROUTINE callgraph_change_capacity

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
      FUNCTION callgraph_get(hash_map, key, default_value) RESULT(value)
         TYPE(callgraph_type), INTENT(in)              :: hash_map
         INTEGER(kind=int_4), DIMENSION(2), INTENT(in)                            :: key
         TYPE(call_stat_type), POINTER, INTENT(in), OPTIONAL                :: default_value
         TYPE(call_stat_type), POINTER                                      :: value
         TYPE(private_item_type), POINTER                    :: item
         INTEGER(KIND=int_8)                                 :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = callgraph_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (callgraph_keys_equal(item%key, key)) THEN
                  value =>item%value
                  RETURN
               END IF
            END IF
            item => item%next
         END DO

         IF (PRESENT(default_value)) THEN
            value =>default_value
            RETURN
         END IF

         CPABORT("Key not found.")
      END FUNCTION callgraph_get

! **************************************************************************************************
!> \brief Remove the value for a given key from the hash map.
!>        If the key is not found the program stops.
!> \param hash_map ...
!> \param key ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE callgraph_del(hash_map, key)
         TYPE(callgraph_type), INTENT(inout)               :: hash_map
         INTEGER(kind=int_4), DIMENSION(2), INTENT(in)                                :: key
         TYPE(private_item_type), POINTER                        :: item, prev_item
         INTEGER(KIND=int_8)                                     :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         hash = callgraph_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         prev_item => Null()
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (callgraph_keys_equal(item%key, key)) THEN
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
      END SUBROUTINE callgraph_del

! **************************************************************************************************
!> \brief Returns the number of key/value-items currently stored in the hash map.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION callgraph_size(hash_map) RESULT(size)
         TYPE(callgraph_type), INTENT(IN)  :: hash_map
         INTEGER                                 :: size

         CPASSERT(ASSOCIATED(hash_map%buckets))
         size = hash_map%size
      END FUNCTION callgraph_size

! **************************************************************************************************
!> \brief Checks whether a given key is currently stored in the hash_map.
!> \param hash_map ...
!> \param key ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION callgraph_haskey(hash_map, key) RESULT(res)
         TYPE(callgraph_type), INTENT(IN)                :: hash_map
         INTEGER(kind=int_4), DIMENSION(2), INTENT(IN)                              :: key
         LOGICAL                                               :: res
         TYPE(private_item_type), POINTER                      :: item
         INTEGER(KIND=int_8)                                   :: hash, idx

         CPASSERT(ASSOCIATED(hash_map%buckets))

         res = .FALSE.
         IF (hash_map%size == 0) RETURN

         hash = callgraph_hash_function(key)
         idx = MOD(hash, INT(size(hash_map%buckets), KIND=int_8)) + 1

         item => hash_map%buckets(idx)%p
         DO WHILE (ASSOCIATED(item))
            IF (item%hash == hash) THEN
               IF (callgraph_keys_equal(item%key, key)) THEN
                  res = .TRUE.
                  return
               END IF
            END IF
            item => item%next
         END DO

      END FUNCTION callgraph_haskey

! **************************************************************************************************
!> \brief Returns a pointer to an array of all key/value-items stored in the hash map.
!>        Caution: The caller is responsible for deallocating targeted array after usage.
!> \param hash_map ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
      FUNCTION callgraph_items(hash_map) RESULT(items)
         TYPE(callgraph_type), INTENT(IN)                 :: hash_map
         TYPE(callgraph_item_type), DIMENSION(:), POINTER :: items

         TYPE(private_item_type), POINTER  :: item
         INTEGER :: i, j

         CPASSERT(ASSOCIATED(hash_map%buckets))

         ALLOCATE (items(hash_map%size))
         j = 1
         DO i = 1, size(hash_map%buckets)
            item => hash_map%buckets(i)%p
            DO WHILE (ASSOCIATED(item))
               items(j)%key =item%key
               items(j)%value =>item%value
               j = j + 1
               item => item%next
            END DO
         END DO

         CPASSERT(j == hash_map%size + 1)
      END FUNCTION callgraph_items

! **************************************************************************************************
!> \brief Copies all key/values-items from one hash map to another.
!>        Afterwards hash_map will contain all items from the from_hash_map and
!>        additionally all its previous items, which were not overwritten.
!>        The two hash maps have to be of the same type.
!> \param hash_map destination of items
!> \param from_hash_map source of items - will not be change
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE callgraph_update(hash_map, from_hash_map)
         TYPE(callgraph_type), INTENT(inout)              :: hash_map
         TYPE(callgraph_type), INTENT(in)                 :: from_hash_map
         TYPE(callgraph_item_type), DIMENSION(:), POINTER :: from_items
         INTEGER :: i

         CPASSERT(ASSOCIATED(hash_map%buckets))
         CPASSERT(ASSOCIATED(from_hash_map%buckets))

         from_items => callgraph_items(from_hash_map)
         DO i = 1, size(from_items)
            CALL callgraph_set(hash_map, from_items(i)%key, from_items(i)%value)
         END DO
         DEALLOCATE (from_items)
      END SUBROUTINE callgraph_update


! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION callgraph_hash_function(key) RESULT(hash)
      INTEGER(kind=int_4), DIMENSION(2), INTENT(in)      :: key
      INTEGER(KIND=int_8)                                :: hash

      INTEGER(kind=int_8)                                :: k1, k2
      k1 = key(1) ! cast to int_8
      k2 = key(2)
      hash = IOR(k1, ISHFT(k2, 32))
   END FUNCTION callgraph_hash_function

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION callgraph_keys_equal(key1, key2) RESULT(res)
      INTEGER(kind=int_4), DIMENSION(2), INTENT(in)      :: key1, key2
      LOGICAL                                            :: res

      res = ALL(key1 == key2)
   END FUNCTION callgraph_keys_equal

END MODULE callgraph
