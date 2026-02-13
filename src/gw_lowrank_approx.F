MODULE gw_lowrank_approx
   USE gw_utils,                        ONLY: add_R,&
                                              analyt_conti_and_print,&
                                              de_init_bs_env,&
                                              get_V_tr_R,&
                                              is_cell_in_index_to_cell,&
                                              power,&
                                              time_to_freq,&
                                              write_reals_to_file,&
                                              write_complex_matrix_to_file
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              default_string_length
   USE kpoint_k_r_trafo_simple,         ONLY: add_kp_to_all_rs,&
                                              fm_add_kp_to_all_rs,&
                                              fm_rs_to_kp,&
                                              rs_to_kp
   USE machine,                         ONLY: m_walltime, m_mkdir
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mathlib,                         ONLY: gemm_square,&
                                              diag_complex
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_lowrank_approx'
   PUBLIC :: calculate_W_lowrank_projector, W_storage_pipeline, bs_env_do_LR

CONTAINS


! **************************************************************************************************
!> \brief Calculates the projector L such that M^a W M^a ⋍ L w L^T with w = L^T M^a W M^a L and a = low_rank_approx_M_exponent
!> \param training_sum The n_RI x n_RI training sum matrix (M^a W M^a) from which the low-rank projector is constructed
!>                     It is symmetrized in-place before the eigenvalue decomposition to ensure hermiticity.
!> \param L Complex projector matrix (n_RI x count) to be constructed
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_W_lowrank_projector(training_sum, L, bs_env, dbg_filenr)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)               :: training_sum
      TYPE(post_scf_bandstructure_type), POINTER                     :: bs_env
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT)    :: L
      INTEGER, OPTIONAL                                              :: dbg_filenr
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)                    :: eigenvectors

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)                       :: eigenvalues

      CHARACTER(LEN=*), PARAMETER                                    :: routineN = 'calculate_W_lowrank_projector'
      INTEGER                                                        :: handle, count, i, j, n_RI, ierr
      CHARACTER(LEN=default_string_length)                           :: filename

      REAL(KIND=dp) :: real_part_norm_trainingsum, imag_part_norm_trainingsum

      CALL timeset(routineN, handle)

      ! Check if rank or epsilon are given.
      IF(bs_env%low_rank_approx_W_rank > 0 .AND. bs_env%low_rank_approx_W_eps > 0.0_dp &
         .AND. bs_env%para_env%is_source()) THEN
         CPWARN("Both LOW_RANK_APPROX_W_RANK and ..._EPS are set. RANK used, EPS ignored.")
      END IF

      n_RI = SIZE(training_sum, 1)
      ALLOCATE(eigenvalues(n_RI))
      ALLOCATE(eigenvectors(n_RI, n_RI))

      ! IF (bs_env%unit_nr > 0) THEN
      !    PRINT *, "Max asymmetry in training sum (should be close to zero): ", &
      !       MAXVAL(ABS(training_sum - CONJG(TRANSPOSE(training_sum))))
      ! END IF

      ! Ensure training_sum is hermitian
      training_sum(:, :) = 0.5_dp * (training_sum(:,:) +  CONJG(TRANSPOSE(training_sum(:,:))))
         

      ! eigenvalue decomposition of the training sum
      CALL diag_complex(training_sum, eigenvectors, eigenvalues)

      real_part_norm_trainingsum = SUM(ABS(REAL(training_sum)))
      imag_part_norm_trainingsum = SUM(ABS(AIMAG(training_sum)))

      ! IF (bs_env%unit_nr > 0) THEN
      !    WRITE (bs_env%unit_nr, *) &
      !       "Imaginary part of training sum: ", imag_part_norm_trainingsum
      !    WRITE (bs_env%unit_nr, *) &
      !       "Real part of training sum: ", real_part_norm_trainingsum
      ! END IF

      IF (bs_env%low_rank_approx_W_rank > 0) THEN
         ! fixed rank mode
         count = MIN(bs_env%low_rank_approx_W_rank, n_RI)
      ELSE
         ! Epsilon threshold mode
         count = 0
         DO i = 1, n_RI
            IF (eigenvalues(i) > bs_env%low_rank_approx_W_eps) count = count + 1
         END DO
      END IF

      ! Construct projector L from eigenvectors with eigenvalues > threshold
      ALLOCATE (L(n_RI, count))

      DO j = 1, count
         i = n_RI - count + j ! Discard the n_RI - count smallest eigenvalues
         L(:, j) = eigenvectors(:, i)
      END DO

      IF (bs_env%para_env%is_source() .AND. PRESENT(dbg_filenr)) THEN
         CALL m_mkdir("W_kpoint_sum_evals", ierr)
         WRITE(filename, '(A,I4.4,A)') "W_kpoint_sum_evals/kp", dbg_filenr, ".txt"
         CALL write_reals_to_file(eigenvalues, TRIM(filename))

         CALL m_mkdir("W_kpoint_projectors", ierr)
         WRITE(filename, '(A,I4.4,A)') "W_kpoint_projectors/kp", dbg_filenr, ".txt"
         CALL write_complex_matrix_to_file(L, TRIM(filename))
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_W_lowrank_projector



! **************************************************************************************************
!> \brief Pipeline for storing/retrieving Real-Space Interaction W^R data
!> \param bs_env The bandstructure environment
!> \param freq_idx The index of the current frequency point (used for filename)
!> \param W_k The Input/Output large real-space matrix (n_RI x n_RI x n_images)
!> \param mode Action to perform: "write", "read", or "read_and_delete"
!>             - "write": Save W_k to disk (training phase)
!>             - "read": Load W_k from disk and broadcast
!>             - "read_and_delete": read W_k and delete the file
! **************************************************************************************************
   SUBROUTINE W_storage_pipeline(bs_env, freq_idx, W_k, mode)
      TYPE(post_scf_bandstructure_type), POINTER           :: bs_env
      INTEGER, INTENT(IN)                                  :: freq_idx
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)  :: W_k
      CHARACTER(LEN=*), INTENT(IN)                         :: mode

      CHARACTER(LEN=default_string_length)                 :: filename
      INTEGER                                              :: unit_nr, iostat
      LOGICAL                                              :: is_master, file_exists

      is_master = (bs_env%para_env%is_source())
      WRITE(filename, '(A,I3.3,A)') "W_realspace_freq_", freq_idx, ".dat"

      SELECT CASE (TRIM(mode))
      CASE ("write", "put", "send")
         ! WRITE MODE: Save data to disk (training phase)
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="WRITE", file_form="UNFORMATTED", &
                           file_status="REPLACE")
            WRITE(unit_nr) W_k
            CALL close_file(unit_nr)
         END IF
         CALL bs_env%para_env%sync()

      CASE ("read_and_delete", "get_and_delete", "receive_and_delete")
         ! READ AND DELETE MODE: Load, delete data from disk and broadcast
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="READ", file_form="UNFORMATTED", &
                           file_status="OLD")
            READ(unit_nr) W_k
            CALL close_file(unit_nr)
         END IF
         
         ! Broadcast the large W_R to all cores
         CALL bs_env%para_env%bcast(W_k)

         IF (is_master) THEN
            ! Delete the file immediately after writing
            INQUIRE(FILE=filename, EXIST=file_exists)
            IF (file_exists) THEN
               OPEN(NEWUNIT=unit_nr, FILE=filename, STATUS='OLD', IOSTAT=iostat)
               IF (iostat == 0) CLOSE(unit_nr, STATUS='DELETE')
            END IF
         END IF

         CALL bs_env%para_env%sync()

      CASE ("read", "get", "receive")
         ! READ MODE: Load data from disk and broadcast
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="READ", file_form="UNFORMATTED", &
                           file_status="OLD")
            READ(unit_nr) W_k
            CALL close_file(unit_nr)
         END IF
         
         ! Broadcast the large W_R to all cores
         CALL bs_env%para_env%bcast(W_k)

      CASE DEFAULT
         CPABORT("Invalid mode for W_storage_pipeline. Use 'write', 'read', or 'read_and_delete'.")
      END SELECT

   END SUBROUTINE W_storage_pipeline


   ! **************************************************************************************************
   !> \brief Checks if a low-rank approximation for W is requested.
   !> \param bs_env The bandstructure environment
   !> \return use_lr True if rank or epsilon threshold is set
   ! **************************************************************************************************
      PURE FUNCTION bs_env_do_LR(bs_env) RESULT(use_lr)
         TYPE(post_scf_bandstructure_type), INTENT(IN)      :: bs_env
         LOGICAL                                            :: use_lr

         use_lr = (bs_env%low_rank_approx_W_rank > 0) .OR. &
                  (bs_env%low_rank_approx_W_eps > 0.0_dp)

      END FUNCTION bs_env_do_LR

END MODULE gw_lowrank_approx