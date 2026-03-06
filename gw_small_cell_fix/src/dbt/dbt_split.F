!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to split blocks and to convert between tensors with different block sizes.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_split



   USE dbt_allocate_wrap, ONLY: allocate_any
   USE dbt_array_list_methods, ONLY: get_ith_array
   USE dbt_block, ONLY: dbt_iterator_type, &
                        dbt_get_block, &
                        dbt_put_block, &
                        dbt_iterator_start, &
                        dbt_iterator_num_blocks, &
                        dbt_iterator_blocks_left, &
                        dbt_iterator_stop, &
                        dbt_iterator_next_block, &
                        dbt_reserve_blocks
   USE dbt_index, ONLY: dbt_get_mapping_info, &
                        dbt_inverse_order
   USE dbt_types, ONLY: dbt_create, &
                        dbt_type, &
                        ndims_tensor, &
                        dbt_distribution_type, &
                        dbt_distribution, &
                        dbt_distribution_destroy, &
                        dbt_distribution_new_expert, &
                        dbt_clear, &
                        dbt_finalize, &
                        dbt_get_num_blocks, &
                        dbt_blk_offsets, &
                        dbt_blk_sizes, &
                        ndims_matrix_row, &
                        ndims_matrix_column, &
                        dbt_filter, &
                        dbt_copy_contraction_storage
   USE kinds, ONLY: dp, dp

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_split'
   PUBLIC :: &
      dbt_make_compatible_blocks, &
      dbt_split_blocks, &
      dbt_split_blocks_generic, &
      dbt_split_copyback, &
      dbt_crop

CONTAINS

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (splitted blocks)
!> \param blk_size_i block sizes for each of the tensor dimensions
!> \param nodata don't copy data from tensor_in to tensor_out
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_split_blocks_generic(tensor_in, tensor_out, blk_size_1, blk_size_2, blk_size_3, blk_size_4, nodata)
      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: blk_size_1, blk_size_2, blk_size_3, blk_size_4
      LOGICAL, INTENT(IN), OPTIONAL                   :: nodata

      TYPE(dbt_distribution_type)                 :: dist_old, dist_split
      TYPE(dbt_iterator_type)                     :: iter
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nd_dist_split_1, nd_dist_split_2, nd_dist_split_3, nd_dist_split_4
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nd_blk_size_split_1, nd_blk_size_split_2, nd_blk_size_split_3,&
          & nd_blk_size_split_4
      INTEGER, DIMENSION(:), ALLOCATABLE              :: index_split_offset_1, index_split_offset_2, index_split_offset_3,&
          & index_split_offset_4
      INTEGER, DIMENSION(:), ALLOCATABLE              :: inblock_offset_1, inblock_offset_2, inblock_offset_3, inblock_offset_4
      INTEGER, DIMENSION(:), ALLOCATABLE              :: blk_nsplit_1, blk_nsplit_2, blk_nsplit_3, blk_nsplit_4
      INTEGER                                         :: split_blk_1, split_blk_2, split_blk_3, split_blk_4
      INTEGER :: idim, i, isplit_sum, nsplit, handle, splitsum, bcount
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: blks_to_allocate
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_d, blk_size_d, blk_size_split_d, dist_split_d
      INTEGER, DIMENSION(ndims_matrix_row(tensor_in)) :: map1_2d
      INTEGER, DIMENSION(ndims_matrix_column(tensor_in)) :: map2_2d
      INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: blk_index, blk_size, blk_offset, &
                                                     blk_shape
      INTEGER, DIMENSION(4) :: bi_split, inblock_offset
      LOGICAL :: found

         REAL(dp), DIMENSION(:,:), ALLOCATABLE :: block_2d
         REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: block_3d
         REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: block_4d
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_split_blocks_generic'

      CALL timeset(routineN, handle)

      dist_old = dbt_distribution(tensor_in)

      DO idim = 1, ndims_tensor(tensor_in)
         CALL get_ith_array(dist_old%nd_dist, idim, dist_d)
         CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

            IF (idim == 1) THEN
               ! split block index offset for each normal block index
               ALLOCATE (index_split_offset_1 (SIZE(dist_d)))
               ! how many split blocks for each normal block index
               ALLOCATE (blk_nsplit_1 (SIZE(dist_d)))
               ! data offset of split blocks w.r.t. corresponding normal block
               ALLOCATE (inblock_offset_1 (SIZE(blk_size_1)))
               ALLOCATE (blk_size_split_d, source=blk_size_1)
            END IF
            IF (idim == 2) THEN
               ! split block index offset for each normal block index
               ALLOCATE (index_split_offset_2 (SIZE(dist_d)))
               ! how many split blocks for each normal block index
               ALLOCATE (blk_nsplit_2 (SIZE(dist_d)))
               ! data offset of split blocks w.r.t. corresponding normal block
               ALLOCATE (inblock_offset_2 (SIZE(blk_size_2)))
               ALLOCATE (blk_size_split_d, source=blk_size_2)
            END IF
            IF (idim == 3) THEN
               ! split block index offset for each normal block index
               ALLOCATE (index_split_offset_3 (SIZE(dist_d)))
               ! how many split blocks for each normal block index
               ALLOCATE (blk_nsplit_3 (SIZE(dist_d)))
               ! data offset of split blocks w.r.t. corresponding normal block
               ALLOCATE (inblock_offset_3 (SIZE(blk_size_3)))
               ALLOCATE (blk_size_split_d, source=blk_size_3)
            END IF
            IF (idim == 4) THEN
               ! split block index offset for each normal block index
               ALLOCATE (index_split_offset_4 (SIZE(dist_d)))
               ! how many split blocks for each normal block index
               ALLOCATE (blk_nsplit_4 (SIZE(dist_d)))
               ! data offset of split blocks w.r.t. corresponding normal block
               ALLOCATE (inblock_offset_4 (SIZE(blk_size_4)))
               ALLOCATE (blk_size_split_d, source=blk_size_4)
            END IF

         ! distribution vector for split blocks
         ALLOCATE (dist_split_d(SIZE(blk_size_split_d)))

         isplit_sum = 0 ! counting splits
         DO i = 1, SIZE(blk_size_d)
            nsplit = 0 ! number of splits for current normal block
            splitsum = 0 ! summing split block sizes for current normal block
            DO WHILE (splitsum < blk_size_d(i))
               nsplit = nsplit + 1
               isplit_sum = isplit_sum + 1
                  IF (idim == 1) inblock_offset_1 (isplit_sum) = splitsum
                  IF (idim == 2) inblock_offset_2 (isplit_sum) = splitsum
                  IF (idim == 3) inblock_offset_3 (isplit_sum) = splitsum
                  IF (idim == 4) inblock_offset_4 (isplit_sum) = splitsum
               dist_split_d(isplit_sum) = dist_d(i)
               splitsum = splitsum + blk_size_split_d(isplit_sum)
            END DO
            CPASSERT(splitsum == blk_size_d(i))
               IF (idim == 1) THEN
                  blk_nsplit_1 (i) = nsplit
                  index_split_offset_1 (i) = isplit_sum - nsplit
               END IF
               IF (idim == 2) THEN
                  blk_nsplit_2 (i) = nsplit
                  index_split_offset_2 (i) = isplit_sum - nsplit
               END IF
               IF (idim == 3) THEN
                  blk_nsplit_3 (i) = nsplit
                  index_split_offset_3 (i) = isplit_sum - nsplit
               END IF
               IF (idim == 4) THEN
                  blk_nsplit_4 (i) = nsplit
                  index_split_offset_4 (i) = isplit_sum - nsplit
               END IF
         END DO

            IF (idim == 1) THEN
               ALLOCATE (nd_dist_split_1, source=dist_split_d)
               ALLOCATE (nd_blk_size_split_1, source=blk_size_split_d)
            END IF
            IF (idim == 2) THEN
               ALLOCATE (nd_dist_split_2, source=dist_split_d)
               ALLOCATE (nd_blk_size_split_2, source=blk_size_split_d)
            END IF
            IF (idim == 3) THEN
               ALLOCATE (nd_dist_split_3, source=dist_split_d)
               ALLOCATE (nd_blk_size_split_3, source=blk_size_split_d)
            END IF
            IF (idim == 4) THEN
               ALLOCATE (nd_dist_split_4, source=dist_split_d)
               ALLOCATE (nd_blk_size_split_4, source=blk_size_split_d)
            END IF
         DEALLOCATE (dist_split_d)
         DEALLOCATE (blk_size_split_d)

      END DO

      CALL dbt_get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)

         IF (ndims_tensor(tensor_in) == 2) THEN
            CALL dbt_distribution_new_expert(dist_split, tensor_in%pgrid, map1_2d, map2_2d, &
                                             nd_dist_split_1, nd_dist_split_2)
            CALL dbt_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                            nd_blk_size_split_1, nd_blk_size_split_2)
         END IF
         IF (ndims_tensor(tensor_in) == 3) THEN
            CALL dbt_distribution_new_expert(dist_split, tensor_in%pgrid, map1_2d, map2_2d, &
                                             nd_dist_split_1, nd_dist_split_2, nd_dist_split_3)
            CALL dbt_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                            nd_blk_size_split_1, nd_blk_size_split_2, nd_blk_size_split_3)
         END IF
         IF (ndims_tensor(tensor_in) == 4) THEN
            CALL dbt_distribution_new_expert(dist_split, tensor_in%pgrid, map1_2d, map2_2d, &
                                             nd_dist_split_1, nd_dist_split_2, nd_dist_split_3, nd_dist_split_4)
            CALL dbt_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                            nd_blk_size_split_1, nd_blk_size_split_2, nd_blk_size_split_3, nd_blk_size_split_4)
         END IF

      CALL dbt_distribution_destroy(dist_split)

      IF (PRESENT(nodata)) THEN
         IF (nodata) THEN
            CALL timestop(handle)
            RETURN
         END IF
      END IF

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor_in,tensor_out) &
!$OMP SHARED(blk_nsplit_1, blk_nsplit_2, blk_nsplit_3, blk_nsplit_4) &
!$OMP SHARED(inblock_offset_1, inblock_offset_2, inblock_offset_3, inblock_offset_4) &
!$OMP SHARED(blk_size_1, blk_size_2, blk_size_3, blk_size_4) &
!$OMP SHARED(index_split_offset_1, index_split_offset_2, index_split_offset_3, index_split_offset_4) &
!$OMP PRIVATE(iter,found,bcount,blks_to_allocate,bi_split,inblock_offset) &
!$OMP PRIVATE(blk_index,blk_size,blk_offset,blk_shape) &
!$OMP PRIVATE(block_2d,block_3d,block_4d)
      CALL dbt_iterator_start(iter, tensor_in)

      bcount = 0
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size)
            IF (ndims_tensor(tensor_in) == 2) THEN
                  DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                  DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                     bcount = bcount + 1
                        END DO
                        END DO
                  END IF
            IF (ndims_tensor(tensor_in) == 3) THEN
                  DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                  DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                  DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                     bcount = bcount + 1
                        END DO
                        END DO
                        END DO
                  END IF
            IF (ndims_tensor(tensor_in) == 4) THEN
                  DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                  DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                  DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                  DO split_blk_4 = 1, blk_nsplit_4 (blk_index(4))
                     bcount = bcount + 1
                        END DO
                        END DO
                        END DO
                        END DO
                  END IF
            END DO
            CALL dbt_iterator_stop(iter)

            ALLOCATE (blks_to_allocate(bcount, ndims_tensor(tensor_in)))

            CALL dbt_iterator_start(iter, tensor_in)

            bcount = 0
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

                  IF (ndims_tensor(tensor_in) == 2) THEN
                        DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                           bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                        DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                           bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                           bcount = bcount + 1
                           blks_to_allocate(bcount, :) = bi_split(1:ndims_tensor(tensor_in))
                              END DO
                              END DO
                        END IF
                  IF (ndims_tensor(tensor_in) == 3) THEN
                        DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                           bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                        DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                           bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                        DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                           bi_split(3) = index_split_offset_3 (blk_index(3)) + split_blk_3
                           bcount = bcount + 1
                           blks_to_allocate(bcount, :) = bi_split(1:ndims_tensor(tensor_in))
                              END DO
                              END DO
                              END DO
                        END IF
                  IF (ndims_tensor(tensor_in) == 4) THEN
                        DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                           bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                        DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                           bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                        DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                           bi_split(3) = index_split_offset_3 (blk_index(3)) + split_blk_3
                        DO split_blk_4 = 1, blk_nsplit_4 (blk_index(4))
                           bi_split(4) = index_split_offset_4 (blk_index(4)) + split_blk_4
                           bcount = bcount + 1
                           blks_to_allocate(bcount, :) = bi_split(1:ndims_tensor(tensor_in))
                              END DO
                              END DO
                              END DO
                              END DO
                        END IF
                  END DO

                  CALL dbt_iterator_stop(iter)

                  CALL dbt_reserve_blocks(tensor_out, blks_to_allocate)

                  CALL dbt_iterator_start(iter, tensor_in)

                  DO WHILE (dbt_iterator_blocks_left(iter))
                     CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)
                        IF (ndims_tensor(tensor_in) == 2) THEN
                           CALL dbt_get_block(tensor_in, blk_index, block_2d, found)
                           CPASSERT(found)
                        END IF
                        IF (ndims_tensor(tensor_in) == 3) THEN
                           CALL dbt_get_block(tensor_in, blk_index, block_3d, found)
                           CPASSERT(found)
                        END IF
                        IF (ndims_tensor(tensor_in) == 4) THEN
                           CALL dbt_get_block(tensor_in, blk_index, block_4d, found)
                           CPASSERT(found)
                        END IF
                        IF (ndims_tensor(tensor_in) == 2) THEN
                              DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                                 ! split block index
                                 bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                                 blk_shape(1) = blk_size_1 (bi_split(1))
                              DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                                 ! split block index
                                 bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                                 blk_shape(2) = blk_size_2 (bi_split(2))

                                    inblock_offset(1) = inblock_offset_1 (bi_split(1))
                                    inblock_offset(2) = inblock_offset_2 (bi_split(2))
                                 CALL dbt_put_block(tensor_out, bi_split(1:2), &
                                                    blk_shape, &
                                                    block_2d( &
                                      inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) +&
                                          & 1:inblock_offset(2) + blk_shape(2)))

                                    END DO
                                    END DO

                                 DEALLOCATE (block_2d)
                              END IF
                        IF (ndims_tensor(tensor_in) == 3) THEN
                              DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                                 ! split block index
                                 bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                                 blk_shape(1) = blk_size_1 (bi_split(1))
                              DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                                 ! split block index
                                 bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                                 blk_shape(2) = blk_size_2 (bi_split(2))
                              DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                                 ! split block index
                                 bi_split(3) = index_split_offset_3 (blk_index(3)) + split_blk_3
                                 blk_shape(3) = blk_size_3 (bi_split(3))

                                    inblock_offset(1) = inblock_offset_1 (bi_split(1))
                                    inblock_offset(2) = inblock_offset_2 (bi_split(2))
                                    inblock_offset(3) = inblock_offset_3 (bi_split(3))
                                 CALL dbt_put_block(tensor_out, bi_split(1:3), &
                                                    blk_shape, &
                                                    block_3d( &
                                      inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) +&
                                          & 1:inblock_offset(2) + blk_shape(2), inblock_offset(3) + 1:inblock_offset(3) +&
                                          & blk_shape(3)))

                                    END DO
                                    END DO
                                    END DO

                                 DEALLOCATE (block_3d)
                              END IF
                        IF (ndims_tensor(tensor_in) == 4) THEN
                              DO split_blk_1 = 1, blk_nsplit_1 (blk_index(1))
                                 ! split block index
                                 bi_split(1) = index_split_offset_1 (blk_index(1)) + split_blk_1
                                 blk_shape(1) = blk_size_1 (bi_split(1))
                              DO split_blk_2 = 1, blk_nsplit_2 (blk_index(2))
                                 ! split block index
                                 bi_split(2) = index_split_offset_2 (blk_index(2)) + split_blk_2
                                 blk_shape(2) = blk_size_2 (bi_split(2))
                              DO split_blk_3 = 1, blk_nsplit_3 (blk_index(3))
                                 ! split block index
                                 bi_split(3) = index_split_offset_3 (blk_index(3)) + split_blk_3
                                 blk_shape(3) = blk_size_3 (bi_split(3))
                              DO split_blk_4 = 1, blk_nsplit_4 (blk_index(4))
                                 ! split block index
                                 bi_split(4) = index_split_offset_4 (blk_index(4)) + split_blk_4
                                 blk_shape(4) = blk_size_4 (bi_split(4))

                                    inblock_offset(1) = inblock_offset_1 (bi_split(1))
                                    inblock_offset(2) = inblock_offset_2 (bi_split(2))
                                    inblock_offset(3) = inblock_offset_3 (bi_split(3))
                                    inblock_offset(4) = inblock_offset_4 (bi_split(4))
                                 CALL dbt_put_block(tensor_out, bi_split(1:4), &
                                                    blk_shape, &
                                                    block_4d( &
                                      inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) +&
                                          & 1:inblock_offset(2) + blk_shape(2), inblock_offset(3) + 1:inblock_offset(3) +&
                                          & blk_shape(3), inblock_offset(4) + 1:inblock_offset(4) + blk_shape(4)))

                                    END DO
                                    END DO
                                    END DO
                                    END DO

                                 DEALLOCATE (block_4d)
                              END IF
                        END DO
                        CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

                        CALL dbt_finalize(tensor_out)

                        ! remove blocks that are exactly 0, these can occur if a cropping operation was performed before splitting
                        CALL dbt_filter(tensor_out, TINY(0.0_dp))

                        CALL timestop(handle)

                     END SUBROUTINE

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks of maximum size PRODUCT(block_sizes).
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (split blocks)
!> \param block_sizes block sizes for each of the tensor dimensions
!> \param nodata don't copy data from tensor_in to tensor_out
!> \author Patrick Seewald
! **************************************************************************************************
                     SUBROUTINE dbt_split_blocks(tensor_in, tensor_out, block_sizes, nodata)

                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
                        TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
                        INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
                           INTENT(IN)                                   :: block_sizes
                        LOGICAL, INTENT(IN), OPTIONAL                   :: nodata

                        INTEGER, DIMENSION(:), ALLOCATABLE              :: nd_blk_size_split_1, nd_blk_size_split_2,&
                            & nd_blk_size_split_3, nd_blk_size_split_4
                        INTEGER :: idim, i, isplit_sum, blk_remainder, nsplit, isplit
                        INTEGER, DIMENSION(:), ALLOCATABLE :: blk_size_d, blk_size_split_d

                        DO idim = 1, ndims_tensor(tensor_in)
                           CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

                           isplit_sum = 0
                           DO i = 1, SIZE(blk_size_d)
                              nsplit = (blk_size_d(i) + block_sizes(idim) - 1)/block_sizes(idim)
                              isplit_sum = isplit_sum + nsplit
                           END DO

                           ALLOCATE (blk_size_split_d(isplit_sum))

                           isplit_sum = 0
                           DO i = 1, SIZE(blk_size_d)
                              nsplit = (blk_size_d(i) + block_sizes(idim) - 1)/block_sizes(idim)
                              blk_remainder = blk_size_d(i)
                              DO isplit = 1, nsplit
                                 isplit_sum = isplit_sum + 1
                                 blk_size_split_d(isplit_sum) = MIN(block_sizes(idim), blk_remainder)
                                 blk_remainder = blk_remainder - block_sizes(idim)
                              END DO

                           END DO

                              IF (idim == 1) THEN
                                 ALLOCATE (nd_blk_size_split_1, source=blk_size_split_d)
                              END IF
                              IF (idim == 2) THEN
                                 ALLOCATE (nd_blk_size_split_2, source=blk_size_split_d)
                              END IF
                              IF (idim == 3) THEN
                                 ALLOCATE (nd_blk_size_split_3, source=blk_size_split_d)
                              END IF
                              IF (idim == 4) THEN
                                 ALLOCATE (nd_blk_size_split_4, source=blk_size_split_d)
                              END IF
                           DEALLOCATE (blk_size_split_d)
                        END DO

                           IF (ndims_tensor(tensor_in) == 2) CALL dbt_split_blocks_generic(tensor_in, tensor_out, &
                                                                                      nd_blk_size_split_1, nd_blk_size_split_2, &
                                                                                                  nodata=nodata)
                           IF (ndims_tensor(tensor_in) == 3) CALL dbt_split_blocks_generic(tensor_in, tensor_out, &
                                                                                      nd_blk_size_split_1, nd_blk_size_split_2,&
                                                                                          & nd_blk_size_split_3, &
                                                                                                  nodata=nodata)
                           IF (ndims_tensor(tensor_in) == 4) CALL dbt_split_blocks_generic(tensor_in, tensor_out, &
                                                                                      nd_blk_size_split_1, nd_blk_size_split_2,&
                                                                                          & nd_blk_size_split_3, nd_blk_size_split_&
                                                                                          &4, &
                                                                                                  nodata=nodata)

                     END SUBROUTINE

! **************************************************************************************************
!> \brief Copy tensor with split blocks to tensor with original block sizes.
!> \param tensor_split_in tensor with smaller blocks
!> \param tensor_out original tensor
!> \author Patrick Seewald
! **************************************************************************************************
                     SUBROUTINE dbt_split_copyback(tensor_split_in, tensor_out, summation)
                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_split_in
                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_out
                        LOGICAL, INTENT(IN), OPTIONAL                   :: summation
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: first_split_d, last_split_d
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: blk_size_split_d, blk_size_d
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: last_split_1, last_split_2, last_split_3, last_split_4, &
                                                                           first_split_1, first_split_2, first_split_3,&
                                                                               & first_split_4, &
                                                                           split_1, split_2, split_3, split_4
                     INTEGER, DIMENSION(:), ALLOCATABLE              :: inblock_offset_1, inblock_offset_2, inblock_offset_3,&
                         & inblock_offset_4, blk_size_split_1, blk_size_split_2, blk_size_split_3, blk_size_split_4
                        INTEGER, DIMENSION(:, :), ALLOCATABLE            :: blks_to_allocate
                        INTEGER                                         :: idim, iblk, bcount
                        INTEGER                                         :: iblk_1, iblk_2, iblk_3, iblk_4, isplit_sum, splitsum
                        TYPE(dbt_iterator_type)                     :: iter
                        INTEGER, DIMENSION(ndims_tensor(tensor_out)) :: blk_index, blk_size, blk_offset, blk_shape, blk_index_n
                        LOGICAL                                         :: found

                        INTEGER, DIMENSION(4)                  :: inblock_offset
                        INTEGER                                            :: handle
                        CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_split_copyback'
                           REAL(dp), DIMENSION(:,:), ALLOCATABLE :: block_2d
                           REAL(dp), DIMENSION(:,:), ALLOCATABLE :: block_split_2d
                           REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: block_3d
                           REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: block_split_3d
                           REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: block_4d
                           REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: block_split_4d

                        CALL timeset(routineN, handle)
                        CPASSERT(tensor_out%valid)
                        IF (PRESENT(summation)) THEN
                           IF (.NOT. summation) CALL dbt_clear(tensor_out)
                        ELSE
                           CALL dbt_clear(tensor_out)
                        END IF

                        DO idim = 1, ndims_tensor(tensor_split_in)
                           CALL get_ith_array(tensor_split_in%blk_sizes, idim, blk_size_split_d)
                           CALL get_ith_array(tensor_out%blk_sizes, idim, blk_size_d)

                              IF (idim == 1) THEN
                                 ! data offset of split blocks w.r.t. corresponding normal block
                                 ALLOCATE (inblock_offset_1 (SIZE(blk_size_split_d)))
                                 ! normal block index for each split block
                                 ALLOCATE (split_1 (SIZE(blk_size_split_d)))
                              END IF
                              IF (idim == 2) THEN
                                 ! data offset of split blocks w.r.t. corresponding normal block
                                 ALLOCATE (inblock_offset_2 (SIZE(blk_size_split_d)))
                                 ! normal block index for each split block
                                 ALLOCATE (split_2 (SIZE(blk_size_split_d)))
                              END IF
                              IF (idim == 3) THEN
                                 ! data offset of split blocks w.r.t. corresponding normal block
                                 ALLOCATE (inblock_offset_3 (SIZE(blk_size_split_d)))
                                 ! normal block index for each split block
                                 ALLOCATE (split_3 (SIZE(blk_size_split_d)))
                              END IF
                              IF (idim == 4) THEN
                                 ! data offset of split blocks w.r.t. corresponding normal block
                                 ALLOCATE (inblock_offset_4 (SIZE(blk_size_split_d)))
                                 ! normal block index for each split block
                                 ALLOCATE (split_4 (SIZE(blk_size_split_d)))
                              END IF

                           ALLOCATE (last_split_d(SIZE(blk_size_d)))
                           ALLOCATE (first_split_d(SIZE(blk_size_d)))
                           first_split_d(1) = 1
                           isplit_sum = 0
                           DO iblk = 1, SIZE(blk_size_d)
                              splitsum = 0
                              IF (iblk > 1) first_split_d(iblk) = last_split_d(iblk - 1) + 1
                              DO WHILE (splitsum < blk_size_d(iblk))
                                 isplit_sum = isplit_sum + 1
                                    IF (idim == 1) THEN
                                       inblock_offset_1 (isplit_sum) = splitsum
                                       split_1 (isplit_sum) = iblk
                                    END IF
                                    IF (idim == 2) THEN
                                       inblock_offset_2 (isplit_sum) = splitsum
                                       split_2 (isplit_sum) = iblk
                                    END IF
                                    IF (idim == 3) THEN
                                       inblock_offset_3 (isplit_sum) = splitsum
                                       split_3 (isplit_sum) = iblk
                                    END IF
                                    IF (idim == 4) THEN
                                       inblock_offset_4 (isplit_sum) = splitsum
                                       split_4 (isplit_sum) = iblk
                                    END IF
                                 splitsum = splitsum + blk_size_split_d(isplit_sum)
                              END DO
                              CPASSERT(splitsum == blk_size_d(iblk))
                              last_split_d(iblk) = isplit_sum
                           END DO
                              IF (idim == 1) THEN
                                 ALLOCATE (first_split_1, source=first_split_d)
                                 ALLOCATE (last_split_1, source=last_split_d)
                                 ALLOCATE (blk_size_split_1, source=blk_size_split_d)
                              END IF
                              IF (idim == 2) THEN
                                 ALLOCATE (first_split_2, source=first_split_d)
                                 ALLOCATE (last_split_2, source=last_split_d)
                                 ALLOCATE (blk_size_split_2, source=blk_size_split_d)
                              END IF
                              IF (idim == 3) THEN
                                 ALLOCATE (first_split_3, source=first_split_d)
                                 ALLOCATE (last_split_3, source=last_split_d)
                                 ALLOCATE (blk_size_split_3, source=blk_size_split_d)
                              END IF
                              IF (idim == 4) THEN
                                 ALLOCATE (first_split_4, source=first_split_d)
                                 ALLOCATE (last_split_4, source=last_split_d)
                                 ALLOCATE (blk_size_split_4, source=blk_size_split_d)
                              END IF
                           DEALLOCATE (first_split_d, last_split_d)
                           DEALLOCATE (blk_size_split_d, blk_size_d)
                        END DO

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor_split_in,tensor_out,summation) &
!$OMP SHARED(split_1, split_2, split_3, split_4) &
!$OMP SHARED(first_split_1, first_split_2, first_split_3, first_split_4) &
!$OMP SHARED(last_split_1, last_split_2, last_split_3, last_split_4) &
!$OMP SHARED(inblock_offset_1, inblock_offset_2, inblock_offset_3, inblock_offset_4) &
!$OMP PRIVATE(iter,blks_to_allocate,bcount,blk_index_n) &
!$OMP PRIVATE(blk_index,blk_size,blk_shape,blk_offset,inblock_offset,found) &
!$OMP PRIVATE(block_2d,block_3d,block_4d,block_split_2d,block_split_3d,block_split_4d)
                        CALL dbt_iterator_start(iter, tensor_split_in)
                        ALLOCATE (blks_to_allocate(dbt_iterator_num_blocks(iter), ndims_tensor(tensor_split_in)))
                        bcount = 0
                        DO WHILE (dbt_iterator_blocks_left(iter))
                           CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size)
                              IF (ndims_tensor(tensor_out) == 2) THEN
                                    blk_index_n(1) = split_1 (blk_index(1))
                                    blk_index_n(2) = split_2 (blk_index(2))
                              END IF
                              IF (ndims_tensor(tensor_out) == 3) THEN
                                    blk_index_n(1) = split_1 (blk_index(1))
                                    blk_index_n(2) = split_2 (blk_index(2))
                                    blk_index_n(3) = split_3 (blk_index(3))
                              END IF
                              IF (ndims_tensor(tensor_out) == 4) THEN
                                    blk_index_n(1) = split_1 (blk_index(1))
                                    blk_index_n(2) = split_2 (blk_index(2))
                                    blk_index_n(3) = split_3 (blk_index(3))
                                    blk_index_n(4) = split_4 (blk_index(4))
                              END IF
                           blks_to_allocate(bcount + 1, :) = blk_index_n
                           bcount = bcount + 1
                        END DO
                        CALL dbt_iterator_stop(iter)
                        CALL dbt_reserve_blocks(tensor_out, blks_to_allocate)

                        CALL dbt_iterator_start(iter, tensor_out)
                        DO WHILE (dbt_iterator_blocks_left(iter))
                           CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

                              IF (ndims_tensor(tensor_out) == 2) THEN
                                 CALL allocate_any(block_2d, blk_size)
                                 block_2d = 0.0_dp
                            DO iblk_1 = first_split_1 (blk_index(1)), last_split_1 (blk_index(1))
                            DO iblk_2 = first_split_2 (blk_index(2)), last_split_2 (blk_index(2))
                                          inblock_offset(1) = inblock_offset_1 (iblk_1)
                                          inblock_offset(2) = inblock_offset_2 (iblk_2)

                            CALL dbt_get_block(tensor_split_in, [iblk_1, iblk_2], &
                                                          block_split_2d, found)
                                       IF (found) THEN
                                          blk_shape(1:2) = SHAPE(block_split_2d)
                                          block_2d( &
                     inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) + 1:inblock_offset(2) +&
                         & blk_shape(2)) = &
                                             block_split_2d
                                       END IF

                                          END DO
                                          END DO
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_2d, summation=summation)
                                       DEALLOCATE (block_2d)
                                    END IF
                              IF (ndims_tensor(tensor_out) == 3) THEN
                                 CALL allocate_any(block_3d, blk_size)
                                 block_3d = 0.0_dp
                            DO iblk_1 = first_split_1 (blk_index(1)), last_split_1 (blk_index(1))
                            DO iblk_2 = first_split_2 (blk_index(2)), last_split_2 (blk_index(2))
                            DO iblk_3 = first_split_3 (blk_index(3)), last_split_3 (blk_index(3))
                                          inblock_offset(1) = inblock_offset_1 (iblk_1)
                                          inblock_offset(2) = inblock_offset_2 (iblk_2)
                                          inblock_offset(3) = inblock_offset_3 (iblk_3)

                            CALL dbt_get_block(tensor_split_in, [iblk_1, iblk_2, iblk_3], &
                                                          block_split_3d, found)
                                       IF (found) THEN
                                          blk_shape(1:3) = SHAPE(block_split_3d)
                                          block_3d( &
                     inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) + 1:inblock_offset(2) +&
                         & blk_shape(2), inblock_offset(3) + 1:inblock_offset(3) + blk_shape(3)) = &
                                             block_split_3d
                                       END IF

                                          END DO
                                          END DO
                                          END DO
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_3d, summation=summation)
                                       DEALLOCATE (block_3d)
                                    END IF
                              IF (ndims_tensor(tensor_out) == 4) THEN
                                 CALL allocate_any(block_4d, blk_size)
                                 block_4d = 0.0_dp
                            DO iblk_1 = first_split_1 (blk_index(1)), last_split_1 (blk_index(1))
                            DO iblk_2 = first_split_2 (blk_index(2)), last_split_2 (blk_index(2))
                            DO iblk_3 = first_split_3 (blk_index(3)), last_split_3 (blk_index(3))
                            DO iblk_4 = first_split_4 (blk_index(4)), last_split_4 (blk_index(4))
                                          inblock_offset(1) = inblock_offset_1 (iblk_1)
                                          inblock_offset(2) = inblock_offset_2 (iblk_2)
                                          inblock_offset(3) = inblock_offset_3 (iblk_3)
                                          inblock_offset(4) = inblock_offset_4 (iblk_4)

                            CALL dbt_get_block(tensor_split_in, [iblk_1, iblk_2, iblk_3, iblk_4], &
                                                          block_split_4d, found)
                                       IF (found) THEN
                                          blk_shape(1:4) = SHAPE(block_split_4d)
                                          block_4d( &
                     inblock_offset(1) + 1:inblock_offset(1) + blk_shape(1), inblock_offset(2) + 1:inblock_offset(2) +&
                         & blk_shape(2), inblock_offset(3) + 1:inblock_offset(3) + blk_shape(3), inblock_offset(4) +&
                         & 1:inblock_offset(4) + blk_shape(4)) = &
                                             block_split_4d
                                       END IF

                                          END DO
                                          END DO
                                          END DO
                                          END DO
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_4d, summation=summation)
                                       DEALLOCATE (block_4d)
                                    END IF
                              END DO
                              CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

                              CALL timestop(handle)

                           END SUBROUTINE

! **************************************************************************************************
!> \brief split two tensors with same total sizes but different block sizes such that they have
!>        equal block sizes
!> \param move_data memory optimization: transfer data s.t. tensor1 and tensor2 may be empty on return
!> \param tensor1_split tensor 1 with split blocks
!> \param tensor2_split tensor 2 with split blocks
!> \param nodata1 don't copy data of tensor 1
!> \param nodata2 don't copy data of tensor 2
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
                           SUBROUTINE dbt_make_compatible_blocks(tensor1, tensor2, tensor1_split, tensor2_split, &
                                                                 order, nodata1, nodata2, move_data)
                              TYPE(dbt_type), INTENT(INOUT) :: tensor1, tensor2
                              TYPE(dbt_type), INTENT(OUT)   :: tensor1_split, tensor2_split
                              INTEGER, DIMENSION(ndims_tensor(tensor1)), &
                                 INTENT(IN), OPTIONAL                        :: order
                              LOGICAL, INTENT(IN), OPTIONAL     :: nodata1, nodata2, move_data
                          INTEGER, DIMENSION(:), ALLOCATABLE  :: blk_size_split_1_1, blk_size_split_1_2, blk_size_split_1_3,&
                              & blk_size_split_1_4, blk_size_split_2_1, blk_size_split_2_2, blk_size_split_2_3,&
                              & blk_size_split_2_4, &
                                                                     blk_size_d_1, blk_size_d_2, blk_size_d_split
                              INTEGER :: size_sum_1, size_sum_2, size_sum, bind_1, bind_2, isplit, bs, idim, i
                              LOGICAL :: move_prv, nodata1_prv, nodata2_prv
                              INTEGER, DIMENSION(ndims_tensor(tensor1)) :: order_prv

                              IF (PRESENT(move_data)) THEN
                                 move_prv = move_data
                              ELSE
                                 move_prv = .FALSE.
                              END IF

                              IF (PRESENT(nodata1)) THEN
                                 nodata1_prv = nodata1
                              ELSE
                                 nodata1_prv = .FALSE.
                              END IF
                              IF (PRESENT(nodata2)) THEN
                                 nodata2_prv = nodata2
                              ELSE
                                 nodata2_prv = .FALSE.
                              END IF

                              IF (PRESENT(order)) THEN
                                 order_prv(:) = dbt_inverse_order(order)
                              ELSE
                                 order_prv(:) = (/(i, i=1, ndims_tensor(tensor1))/)
                              END IF

                              DO idim = 1, ndims_tensor(tensor2)
                                 CALL get_ith_array(tensor1%blk_sizes, order_prv(idim), blk_size_d_1)
                                 CALL get_ith_array(tensor2%blk_sizes, idim, blk_size_d_2)
                                 ALLOCATE (blk_size_d_split(SIZE(blk_size_d_1) + SIZE(blk_size_d_2)))
                                 size_sum_1 = 0
                                 size_sum_2 = 0
                                 size_sum = 0
                                 bind_1 = 0
                                 bind_2 = 0
                                 isplit = 0

                                 DO WHILE (bind_1 < SIZE(blk_size_d_1) .AND. bind_2 < SIZE(blk_size_d_2))
                                    IF (blk_size_d_1(bind_1 + 1) < blk_size_d_2(bind_2 + 1)) THEN
                                       bind_1 = bind_1 + 1
                                       bs = blk_size_d_1(bind_1)
                                       blk_size_d_2(bind_2 + 1) = blk_size_d_2(bind_2 + 1) - bs
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    ELSEIF (blk_size_d_1(bind_1 + 1) > blk_size_d_2(bind_2 + 1)) THEN
                                       bind_2 = bind_2 + 1
                                       bs = blk_size_d_2(bind_2)
                                       blk_size_d_1(bind_1 + 1) = blk_size_d_1(bind_1 + 1) - bs
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    ELSE
                                       bind_1 = bind_1 + 1
                                       bind_2 = bind_2 + 1
                                       bs = blk_size_d_1(bind_1)
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    END IF
                                 END DO

                                 IF (bind_1 < SIZE(blk_size_d_1)) THEN
                                    bind_1 = bind_1 + 1
                                    bs = blk_size_d_1(bind_1)
                                    size_sum = size_sum + bs
                                    isplit = isplit + 1
                                    blk_size_d_split(isplit) = bs
                                 END IF

                                 IF (bind_2 < SIZE(blk_size_d_2)) THEN
                                    bind_2 = bind_2 + 1
                                    bs = blk_size_d_2(bind_2)
                                    size_sum = size_sum + bs
                                    isplit = isplit + 1
                                    blk_size_d_split(isplit) = bs
                                 END IF

                                    IF (order_prv(idim) == 1) THEN
                                       ALLOCATE (blk_size_split_1_1, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (order_prv(idim) == 2) THEN
                                       ALLOCATE (blk_size_split_1_2, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (order_prv(idim) == 3) THEN
                                       ALLOCATE (blk_size_split_1_3, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (order_prv(idim) == 4) THEN
                                       ALLOCATE (blk_size_split_1_4, source=blk_size_d_split(:isplit))
                                    END IF

                                    IF (idim == 1) THEN
                                       ALLOCATE (blk_size_split_2_1, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (idim == 2) THEN
                                       ALLOCATE (blk_size_split_2_2, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (idim == 3) THEN
                                       ALLOCATE (blk_size_split_2_3, source=blk_size_d_split(:isplit))
                                    END IF
                                    IF (idim == 4) THEN
                                       ALLOCATE (blk_size_split_2_4, source=blk_size_d_split(:isplit))
                                    END IF

                                 DEALLOCATE (blk_size_d_split, blk_size_d_1, blk_size_d_2)
                              END DO

                                 IF (ndims_tensor(tensor1) == 2) THEN
                   CALL dbt_split_blocks_generic(tensor1, tensor1_split, blk_size_split_1_1, blk_size_split_1_2, nodata=nodata1)
                                    IF (move_prv .AND. .NOT. nodata1_prv) CALL dbt_clear(tensor1)
                   CALL dbt_split_blocks_generic(tensor2, tensor2_split, blk_size_split_2_1, blk_size_split_2_2, nodata=nodata2)
                                    IF (move_prv .AND. .NOT. nodata2_prv) CALL dbt_clear(tensor2)
                                 END IF
                                 IF (ndims_tensor(tensor1) == 3) THEN
                   CALL dbt_split_blocks_generic(tensor1, tensor1_split, blk_size_split_1_1, blk_size_split_1_2,&
                       & blk_size_split_1_3, nodata=nodata1)
                                    IF (move_prv .AND. .NOT. nodata1_prv) CALL dbt_clear(tensor1)
                   CALL dbt_split_blocks_generic(tensor2, tensor2_split, blk_size_split_2_1, blk_size_split_2_2,&
                       & blk_size_split_2_3, nodata=nodata2)
                                    IF (move_prv .AND. .NOT. nodata2_prv) CALL dbt_clear(tensor2)
                                 END IF
                                 IF (ndims_tensor(tensor1) == 4) THEN
                   CALL dbt_split_blocks_generic(tensor1, tensor1_split, blk_size_split_1_1, blk_size_split_1_2,&
                       & blk_size_split_1_3, blk_size_split_1_4, nodata=nodata1)
                                    IF (move_prv .AND. .NOT. nodata1_prv) CALL dbt_clear(tensor1)
                   CALL dbt_split_blocks_generic(tensor2, tensor2_split, blk_size_split_2_1, blk_size_split_2_2,&
                       & blk_size_split_2_3, blk_size_split_2_4, nodata=nodata2)
                                    IF (move_prv .AND. .NOT. nodata2_prv) CALL dbt_clear(tensor2)
                                 END IF

                           END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
                           SUBROUTINE dbt_crop(tensor_in, tensor_out, bounds, move_data)
                              TYPE(dbt_type), INTENT(INOUT) :: tensor_in
                              TYPE(dbt_type), INTENT(OUT) :: tensor_out
                              INTEGER, DIMENSION(2, ndims_tensor(tensor_in)), INTENT(IN) :: bounds
                              LOGICAL, INTENT(IN), OPTIONAL :: move_data

                              CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_crop'

                              INTEGER, DIMENSION(2, ndims_tensor(tensor_in)) :: blk_bounds
                              TYPE(dbt_iterator_type)                     :: iter
                              INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: blk_index, blk_size, blk_offset
                              LOGICAL :: found, move_data_prv
                              INTEGER :: handle, idim, iblk_out
                              INTEGER, DIMENSION(:, :), ALLOCATABLE :: blk_ind_out
                                 REAL(dp), DIMENSION(:,:), ALLOCATABLE :: block_2d, block_put_2d
                                 REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: block_3d, block_put_3d
                                 REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: block_4d, block_put_4d

                              CALL timeset(routineN, handle)

                              IF (PRESENT(move_data)) THEN
                                 move_data_prv = move_data
                              ELSE
                                 move_data_prv = .FALSE.
                              END IF

                              CALL dbt_create(tensor_in, tensor_out)

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor_in,tensor_out,bounds) &
!$OMP PRIVATE(iter,blk_ind_out,iblk_out,blk_index,blk_size,blk_offset,found,blk_bounds) &
!$OMP PRIVATE(block_2d,block_put_2d,block_3d,block_put_3d,block_4d,block_put_4d)

                              ! reserve blocks inside bounds
                              CALL dbt_iterator_start(iter, tensor_in)
                              ALLOCATE (blk_ind_out(dbt_iterator_num_blocks(iter), ndims_tensor(tensor_in)))
                              blk_ind_out(:, :) = 0
                              iblk_out = 0
                              blk_loop: DO WHILE (dbt_iterator_blocks_left(iter))
                                 CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)
                                 DO idim = 1, ndims_tensor(tensor_in)
                                    IF (bounds(1, idim) > blk_offset(idim) - 1 + blk_size(idim)) CYCLE blk_loop
                                    IF (bounds(2, idim) < blk_offset(idim)) CYCLE blk_loop
                                 END DO
                                 iblk_out = iblk_out + 1
                                 blk_ind_out(iblk_out, :) = blk_index
                              END DO blk_loop
                              CALL dbt_iterator_stop(iter)

                              CALL dbt_reserve_blocks(tensor_out, blk_ind_out(1:iblk_out, :))
                              DEALLOCATE (blk_ind_out)

                              ! copy blocks
                              CALL dbt_iterator_start(iter, tensor_out)
                              iter_loop: DO WHILE (dbt_iterator_blocks_left(iter))
                                 CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

                                 DO idim = 1, ndims_tensor(tensor_in)
                                    blk_bounds(1, idim) = MAX(bounds(1, idim) - blk_offset(idim) + 1, 1)
                                    blk_bounds(2, idim) = MIN(bounds(2, idim) - blk_offset(idim) + 1, blk_size(idim))
                                 END DO

                                    IF (ndims_tensor(tensor_in) == 2) THEN
                                       CALL dbt_get_block(tensor_in, blk_index, block_2d, found)

                             ALLOCATE (block_put_2d(blk_size(1), blk_size(2)))
                                       block_put_2d = 0.0_dp
   block_put_2d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2)) = &
           block_2d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2))
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_put_2d)
                                       DEALLOCATE (block_2d)
                                       DEALLOCATE (block_put_2d)
                                    END IF
                                    IF (ndims_tensor(tensor_in) == 3) THEN
                                       CALL dbt_get_block(tensor_in, blk_index, block_3d, found)

                             ALLOCATE (block_put_3d(blk_size(1), blk_size(2), blk_size(3)))
                                       block_put_3d = 0.0_dp
   block_put_3d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2), blk_bounds(1, 3):blk_bounds(2,3)) = &
           block_3d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2), blk_bounds(1, 3):blk_bounds(2,3))
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_put_3d)
                                       DEALLOCATE (block_3d)
                                       DEALLOCATE (block_put_3d)
                                    END IF
                                    IF (ndims_tensor(tensor_in) == 4) THEN
                                       CALL dbt_get_block(tensor_in, blk_index, block_4d, found)

                             ALLOCATE (block_put_4d(blk_size(1), blk_size(2), blk_size(3), blk_size(4)))
                                       block_put_4d = 0.0_dp
   block_put_4d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2), blk_bounds(1, 3):blk_bounds(2,3),&
       & blk_bounds(1, 4):blk_bounds(2,4)) = &
           block_4d(blk_bounds(1, 1):blk_bounds(2,1), blk_bounds(1, 2):blk_bounds(2,2), blk_bounds(1, 3):blk_bounds(2,3),&
               & blk_bounds(1, 4):blk_bounds(2,4))
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_put_4d)
                                       DEALLOCATE (block_4d)
                                       DEALLOCATE (block_put_4d)
                                    END IF
                              END DO iter_loop
                              CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
                              CALL dbt_finalize(tensor_out)

                              IF (move_data_prv) CALL dbt_clear(tensor_in)

                              ! transfer data for batched contraction
                              CALL dbt_copy_contraction_storage(tensor_in, tensor_out)

                              CALL timestop(handle)
                           END SUBROUTINE

                        END MODULE
