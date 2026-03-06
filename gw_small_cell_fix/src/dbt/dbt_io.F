!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBT tensor Input / Output
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_io



   USE dbt_types, ONLY: &
      dbt_get_info, dbt_type, ndims_tensor, dbt_get_num_blocks, dbt_get_num_blocks_total, &
      blk_dims_tensor, dbt_get_stored_coordinates, dbt_get_nze, dbt_get_nze_total, &
      dbt_pgrid_type, dbt_nblks_total
   USE kinds, ONLY: default_string_length, int_8, dp
   USE message_passing, ONLY: mp_comm_type
   USE dbt_block, ONLY: &
      dbt_iterator_type, dbt_iterator_next_block, dbt_iterator_start, &
      dbt_iterator_blocks_left, dbt_iterator_stop, dbt_get_block
   USE dbt_tas_io, ONLY: dbt_tas_write_split_info

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_types'

   PUBLIC :: &
      dbt_write_tensor_info, &
      dbt_write_tensor_dist, &
      dbt_write_blocks, &
      dbt_write_block, &
      dbt_write_block_indices, &
      dbt_write_split_info, &
      prep_output_unit

CONTAINS

! **************************************************************************************************
!> \brief Write tensor global info: block dimensions, full dimensions and process grid dimensions
!> \param full_info Whether to print distribution and block size vectors
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_write_tensor_info(tensor, unit_nr, full_info)
      TYPE(dbt_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: unit_nr
      LOGICAL, OPTIONAL, INTENT(IN)  :: full_info
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: nblks_total, nfull_total, pdims, my_ploc, nblks_local, nfull_local

         INTEGER, DIMENSION(dbt_nblks_total(tensor, 1)) :: proc_dist_1
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 1)) :: blk_size_1
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 1)) :: blks_local_1
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 2)) :: proc_dist_2
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 2)) :: blk_size_2
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 2)) :: blks_local_2
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 3)) :: proc_dist_3
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 3)) :: blk_size_3
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 3)) :: blks_local_3
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 4)) :: proc_dist_4
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 4)) :: blk_size_4
         INTEGER, DIMENSION(dbt_nblks_total(tensor, 4)) :: blks_local_4
      CHARACTER(len=default_string_length)                     :: name
      INTEGER                                                  :: idim
      INTEGER                                                  :: iblk
      INTEGER                                                  :: unit_nr_prv

      unit_nr_prv = prep_output_unit(unit_nr)
      IF (unit_nr_prv == 0) RETURN

      CALL dbt_get_info(tensor, nblks_total, nfull_total, nblks_local, nfull_local, pdims, my_ploc, &
                        blks_local_1, blks_local_2, blks_local_3, blks_local_4, proc_dist_1, proc_dist_2, proc_dist_3,&
                            & proc_dist_4, blk_size_1, blk_size_2, blk_size_3, blk_size_4, &
                        name=name)

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, "(T2,A)") &
            "GLOBAL INFO OF "//TRIM(name)
         WRITE (unit_nr_prv, "(T4,A,1X)", advance="no") "block dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (unit_nr_prv, "(I6)", advance="no") nblks_total(idim)
         END DO
         WRITE (unit_nr_prv, "(/T4,A,1X)", advance="no") "full dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (unit_nr_prv, "(I8)", advance="no") nfull_total(idim)
         END DO
         WRITE (unit_nr_prv, "(/T4,A,1X)", advance="no") "process grid dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (unit_nr_prv, "(I6)", advance="no") pdims(idim)
         END DO
         WRITE (unit_nr_prv, *)

         IF (PRESENT(full_info)) THEN
            IF (full_info) THEN
               WRITE (unit_nr_prv, '(T4,A)', advance='no') "Block sizes:"
                  IF (ndims_tensor(tensor) >= 1) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 1, ':'
                     DO iblk = 1, SIZE(blk_size_1)
                        WRITE (unit_nr_prv, '(I2,1X)', advance='no') blk_size_1 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 2) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 2, ':'
                     DO iblk = 1, SIZE(blk_size_2)
                        WRITE (unit_nr_prv, '(I2,1X)', advance='no') blk_size_2 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 3) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 3, ':'
                     DO iblk = 1, SIZE(blk_size_3)
                        WRITE (unit_nr_prv, '(I2,1X)', advance='no') blk_size_3 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 4) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 4, ':'
                     DO iblk = 1, SIZE(blk_size_4)
                        WRITE (unit_nr_prv, '(I2,1X)', advance='no') blk_size_4 (iblk)
                     END DO
                  END IF
               WRITE (unit_nr_prv, '(/T4,A)', advance='no') "Block distribution:"
                  IF (ndims_tensor(tensor) >= 1) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 1, ':'
                     DO iblk = 1, SIZE(proc_dist_1)
                        WRITE (unit_nr_prv, '(I3,1X)', advance='no') proc_dist_1 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 2) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 2, ':'
                     DO iblk = 1, SIZE(proc_dist_2)
                        WRITE (unit_nr_prv, '(I3,1X)', advance='no') proc_dist_2 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 3) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 3, ':'
                     DO iblk = 1, SIZE(proc_dist_3)
                        WRITE (unit_nr_prv, '(I3,1X)', advance='no') proc_dist_3 (iblk)
                     END DO
                  END IF
                  IF (ndims_tensor(tensor) >= 4) THEN
                     WRITE (unit_nr_prv, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', 4, ':'
                     DO iblk = 1, SIZE(proc_dist_4)
                        WRITE (unit_nr_prv, '(I3,1X)', advance='no') proc_dist_4 (iblk)
                     END DO
                  END IF
            END IF
            WRITE (unit_nr_prv, *)
         END IF
      END IF

   END SUBROUTINE dbt_write_tensor_info

! **************************************************************************************************
!> \brief Write info on tensor distribution & load balance
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_write_tensor_dist(tensor, unit_nr)
      TYPE(dbt_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: unit_nr
      INTEGER                        :: nproc, nblock_max, nelement_max
      INTEGER(KIND=int_8)            :: nblock_sum, nelement_sum, nblock_tot
      INTEGER                        :: nblock, nelement, unit_nr_prv
      INTEGER, DIMENSION(2)          :: tmp
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: bdims
      REAL(KIND=dp)              :: occupation

      unit_nr_prv = prep_output_unit(unit_nr)
      IF (unit_nr_prv == 0) RETURN

      nproc = tensor%pgrid%mp_comm_2d%num_pe

      nblock = dbt_get_num_blocks(tensor)
      nelement = dbt_get_nze(tensor)

      nblock_sum = dbt_get_num_blocks_total(tensor)
      nelement_sum = dbt_get_nze_total(tensor)

      tmp = [nblock, nelement]
      CALL tensor%pgrid%mp_comm_2d%max(tmp)
      nblock_max = tmp(1); nelement_max = tmp(2)

      CALL blk_dims_tensor(tensor, bdims)
      nblock_tot = PRODUCT(INT(bdims, KIND=int_8))

      occupation = -1.0_dp
      IF (nblock_tot /= 0) occupation = 100.0_dp*REAL(nblock_sum, dp)/REAL(nblock_tot, dp)

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, "(T2,A)") &
            "DISTRIBUTION OF "//TRIM(tensor%name)
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Number of non-zero blocks:", nblock_sum
         WRITE (unit_nr_prv, "(T15,A,T75,F6.2)") "Percentage of non-zero blocks:", occupation
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of blocks per CPU:", (nblock_sum + nproc - 1)/nproc
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of blocks per CPU:", nblock_max
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Average number of matrix elements per CPU:", (nelement_sum + nproc - 1)/nproc
         WRITE (unit_nr_prv, "(T15,A,T68,I13)") "Maximum number of matrix elements per CPU:", nelement_max
      END IF

   END SUBROUTINE dbt_write_tensor_dist

! **************************************************************************************************
!> \brief Write all tensor blocks
!> \param io_unit_master for global output
!> \param io_unit_all for local output
!> \param write_int convert to integers (useful for testing with integer tensors)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_write_blocks(tensor, io_unit_master, io_unit_all, write_int)
      TYPE(dbt_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: io_unit_master, io_unit_all
      LOGICAL, INTENT(IN), OPTIONAL                      :: write_int
      INTEGER, DIMENSION(ndims_tensor(tensor))          :: blk_index, blk_size
         REAL(KIND=dp), ALLOCATABLE, &
            DIMENSION(:,:)                :: blk_values_2
         REAL(KIND=dp), ALLOCATABLE, &
            DIMENSION(:,:,:)                :: blk_values_3
         REAL(KIND=dp), ALLOCATABLE, &
            DIMENSION(:,:,:,:)                :: blk_values_4
      TYPE(dbt_iterator_type)                        :: iterator
      INTEGER                                            :: proc, mynode
      LOGICAL                                            :: found

      IF (io_unit_master > 0) THEN
         WRITE (io_unit_master, '(T7,A)') "(block index) @ process: (array index) value"
      END IF
      CALL dbt_iterator_start(iterator, tensor)
      DO WHILE (dbt_iterator_blocks_left(iterator))
         CALL dbt_iterator_next_block(iterator, blk_index, blk_size=blk_size)
         CALL dbt_get_stored_coordinates(tensor, blk_index, proc)
         mynode = tensor%pgrid%mp_comm_2d%mepos
         CPASSERT(proc == mynode)
            IF (ndims_tensor(tensor) == 2) THEN
               CALL dbt_get_block(tensor, blk_index, blk_values_2, found)
               CPASSERT(found)
               CALL dbt_write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, &
                                    blk_values_2=blk_values_2, write_int=write_int)
               DEALLOCATE (blk_values_2)
            END IF
            IF (ndims_tensor(tensor) == 3) THEN
               CALL dbt_get_block(tensor, blk_index, blk_values_3, found)
               CPASSERT(found)
               CALL dbt_write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, &
                                    blk_values_3=blk_values_3, write_int=write_int)
               DEALLOCATE (blk_values_3)
            END IF
            IF (ndims_tensor(tensor) == 4) THEN
               CALL dbt_get_block(tensor, blk_index, blk_values_4, found)
               CPASSERT(found)
               CALL dbt_write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, &
                                    blk_values_4=blk_values_4, write_int=write_int)
               DEALLOCATE (blk_values_4)
            END IF
      END DO
      CALL dbt_iterator_stop(iterator)
   END SUBROUTINE dbt_write_blocks

! **************************************************************************************************
!> \brief Write a tensor block
!> \param name tensor name
!> \param blk_size block size
!> \param blk_index block index
!> \param blk_values_i block values for 2 dimensions
!> \param write_int write_int convert values to integers
!> \param unit_nr unit number
!> \param proc which process am I
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_write_block(name, blk_size, blk_index, proc, unit_nr, &
                              blk_values_2, blk_values_3, blk_values_4, write_int)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_index
         REAL(KIND=dp), &
            DIMENSION(blk_size(1), blk_size(2)), &
            INTENT(IN), OPTIONAL                            :: blk_values_2
         REAL(KIND=dp), &
            DIMENSION(blk_size(1), blk_size(2), blk_size(3)), &
            INTENT(IN), OPTIONAL                            :: blk_values_3
         REAL(KIND=dp), &
            DIMENSION(blk_size(1), blk_size(2), blk_size(3), blk_size(4)), &
            INTENT(IN), OPTIONAL                            :: blk_values_4
      LOGICAL, INTENT(IN), OPTIONAL                      :: write_int
      LOGICAL                                            :: write_int_prv
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, INTENT(IN)                                :: proc
      INTEGER                                            :: i_1, i_2, i_3, i_4
      INTEGER                                            :: ndim

      IF (PRESENT(write_int)) THEN
         write_int_prv = write_int
      ELSE
         write_int_prv = .FALSE.
      END IF

      ndim = SIZE(blk_size)

      IF (unit_nr > 0) THEN
            IF (ndim == 2) THEN
                  DO i_2 = 1, blk_size(2)
                  DO i_1 = 1, blk_size(1)
                     IF (write_int_prv) THEN
                        WRITE (unit_nr, '(T7,A,T16,A,2I3,1X,A,1X,I3,A,1X,A,2I3,1X,A,1X,I20)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, ")", &
                           INT(blk_values_2 (i_1, i_2), KIND=int_8)
                     ELSE
                        WRITE (unit_nr, '(T7,A,T16,A,2I3,1X,A,1X,I3,A,1X,A,2I3,1X,A,1X,F10.5)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, ")", &
                           blk_values_2 (i_1, i_2)
                     END IF
                        END DO
                        END DO
                  END IF
            IF (ndim == 3) THEN
                  DO i_3 = 1, blk_size(3)
                  DO i_2 = 1, blk_size(2)
                  DO i_1 = 1, blk_size(1)
                     IF (write_int_prv) THEN
                        WRITE (unit_nr, '(T7,A,T16,A,3I3,1X,A,1X,I3,A,1X,A,3I3,1X,A,1X,I20)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, i_3, ")", &
                           INT(blk_values_3 (i_1, i_2, i_3), KIND=int_8)
                     ELSE
                        WRITE (unit_nr, '(T7,A,T16,A,3I3,1X,A,1X,I3,A,1X,A,3I3,1X,A,1X,F10.5)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, i_3, ")", &
                           blk_values_3 (i_1, i_2, i_3)
                     END IF
                        END DO
                        END DO
                        END DO
                  END IF
            IF (ndim == 4) THEN
                  DO i_4 = 1, blk_size(4)
                  DO i_3 = 1, blk_size(3)
                  DO i_2 = 1, blk_size(2)
                  DO i_1 = 1, blk_size(1)
                     IF (write_int_prv) THEN
                        WRITE (unit_nr, '(T7,A,T16,A,4I3,1X,A,1X,I3,A,1X,A,4I3,1X,A,1X,I20)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, i_3, i_4, ")", &
                           INT(blk_values_4 (i_1, i_2, i_3, i_4), KIND=int_8)
                     ELSE
                        WRITE (unit_nr, '(T7,A,T16,A,4I3,1X,A,1X,I3,A,1X,A,4I3,1X,A,1X,F10.5)') &
                           TRIM(name), "(", blk_index, ") @", proc, ':', &
                           "(", i_1, i_2, i_3, i_4, ")", &
                           blk_values_4 (i_1, i_2, i_3, i_4)
                     END IF
                        END DO
                        END DO
                        END DO
                        END DO
                  END IF
            END IF
         END SUBROUTINE dbt_write_block

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_write_block_indices(tensor, io_unit_master, io_unit_all)
            TYPE(dbt_type), INTENT(INOUT)                  :: tensor
            INTEGER, INTENT(IN)                                :: io_unit_master, io_unit_all
            TYPE(dbt_iterator_type)                        :: iterator
            INTEGER, DIMENSION(ndims_tensor(tensor))          :: blk_index, blk_size
            INTEGER                                            :: mynode, proc

            IF (io_unit_master > 0) THEN
               WRITE (io_unit_master, '(T7,A)') "(block index) @ process: size"
            END IF

            CALL dbt_iterator_start(iterator, tensor)
            DO WHILE (dbt_iterator_blocks_left(iterator))
               CALL dbt_iterator_next_block(iterator, blk_index, blk_size=blk_size)
               CALL dbt_get_stored_coordinates(tensor, blk_index, proc)
               mynode = tensor%pgrid%mp_comm_2d%mepos
               CPASSERT(proc == mynode)
                  IF (ndims_tensor(tensor) == 2) THEN
                     WRITE (io_unit_all, '(T7,A,T16,A,2I3,1X,A,1X,I3,A2,2I3)') &
                        TRIM(tensor%name), "blk index (", blk_index, ") @", proc, ":", blk_size
                  END IF
                  IF (ndims_tensor(tensor) == 3) THEN
                     WRITE (io_unit_all, '(T7,A,T16,A,3I3,1X,A,1X,I3,A2,3I3)') &
                        TRIM(tensor%name), "blk index (", blk_index, ") @", proc, ":", blk_size
                  END IF
                  IF (ndims_tensor(tensor) == 4) THEN
                     WRITE (io_unit_all, '(T7,A,T16,A,4I3,1X,A,1X,I3,A2,4I3)') &
                        TRIM(tensor%name), "blk index (", blk_index, ") @", proc, ":", blk_size
                  END IF
            END DO
            CALL dbt_iterator_stop(iterator)
         END SUBROUTINE dbt_write_block_indices

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_write_split_info(pgrid, unit_nr)
            TYPE(dbt_pgrid_type), INTENT(IN) :: pgrid
            INTEGER, INTENT(IN) :: unit_nr

            IF (ALLOCATED(pgrid%tas_split_info)) THEN
               CALL dbt_tas_write_split_info(pgrid%tas_split_info, unit_nr)
            END IF
         END SUBROUTINE dbt_write_split_info

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION prep_output_unit(unit_nr) RESULT(unit_nr_out)
            INTEGER, INTENT(IN), OPTIONAL :: unit_nr
            INTEGER                       :: unit_nr_out

            IF (PRESENT(unit_nr)) THEN
               unit_nr_out = unit_nr
            ELSE
               unit_nr_out = 0
            END IF

         END FUNCTION prep_output_unit

      END MODULE dbt_io
