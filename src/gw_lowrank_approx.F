MODULE gw_lowrank_approx
   USE gw_utils,                        ONLY: add_R,&
                                              analyt_conti_and_print,&
                                              de_init_bs_env,&
                                              get_V_tr_R,&
                                              is_cell_in_index_to_cell,&
                                              power,&
                                              time_to_freq,&
                                              write_reals_to_file,&
                                              low_rank_filter
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              default_string_length
   USE kpoint_k_r_trafo_simple,         ONLY: add_kp_to_all_rs,&
                                              fm_add_kp_to_all_rs,&
                                              fm_rs_to_kp,&
                                              rs_to_kp
   USE machine,                         ONLY: m_walltime
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
   PUBLIC :: calculate_W_lowrank_projector, W_storage_pipeline

CONTAINS


! **************************************************************************************************
!> \brief Calculates the projector L such that M^a W M^a ⋍ L w L^T with w = L^T M^a W M^a L and a = low_rank_approx_M_exponent
!> \param W_R ...
!> \param M_R ...
!> \param L Currently complex projector (could be real for efficiency?)
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_W_lowrank_projector(training_sum, L, bs_env)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(IN)     :: training_sum
      TYPE(post_scf_bandstructure_type), POINTER                     :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT)    :: L

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)                       :: eigenvalues
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)                 :: eigenvectors

      CHARACTER(LEN=*), PARAMETER                                    :: routineN = 'calculate_W_lowrank_projector'
      INTEGER                                                        :: handle, count, i, j, n_RI

      REAL(KIND=dp) :: real_part_norm_trainingsum, imag_part_norm_trainingsum

      CALL timeset(routineN, handle)

      n_RI = SIZE(training_sum, 1)
      ALLOCATE(eigenvalues(n_RI), eigenvectors(n_RI, n_RI))
   
      ! eigenvalue decomposition of the training sum
      CALL diag_complex(training_sum, eigenvectors, eigenvalues)

      real_part_norm_trainingsum = SUM(ABS(REAL(training_sum)))

      imag_part_norm_trainingsum = SUM(ABS(AIMAG(training_sum)))

      PRINT *, "imaginary part of training sum (should be close to zero): ", imag_part_norm_trainingsum
      PRINT *, "real part of training sum: ", real_part_norm_trainingsum

      IF (bs_env%para_env%mepos == bs_env%para_env%source) THEN
         CALL write_reals_to_file(eigenvalues, "eigenvalues_W_lowrank_training_sum.txt")
      END IF

      ! Count eigenvalues > W_lowrank_approx_eps
      count = 0
      DO i = 1, n_RI
         IF (eigenvalues(i) > bs_env%low_rank_approx_W_eps) count = count + 1
      END DO

      ! Construct projector L from eigenvectors with eigenvalues > threshold
      ALLOCATE (L(n_RI, count))
      j = 0
      DO i = 1, n_RI
         IF (eigenvalues(i) > bs_env%low_rank_approx_W_eps) THEN
            j = j + 1
            L(:, j) = REAL(eigenvectors(:, i))
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_W_lowrank_projector



! **************************************************************************************************
!> \brief Pipeline for storing/retrieving Real-Space Interaction W^R data
!> \param bs_env The bandstructure environment
!> \param freq_idx The index of the current frequency point (used for filename)
!> \param MWM_R The Input/Output large real-space matrix (n_RI x n_RI x n_images)
!> \param mode Action to perform: "write", "read", or "read_and_delete"
!>             - "write": Save MWM_R to disk (training phase)
!>             - "read": Load MWM_R from disk and broadcast
!>             - "read_and_delete": read MWM_R and delete the file
! **************************************************************************************************
   SUBROUTINE W_storage_pipeline(bs_env, freq_idx, MWM_R, mode)
      TYPE(post_scf_bandstructure_type), POINTER           :: bs_env
      INTEGER, INTENT(IN)                                  :: freq_idx
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)  :: MWM_R
      CHARACTER(LEN=*), INTENT(IN)                         :: mode

      CHARACTER(LEN=default_string_length)                 :: filename
      INTEGER                                              :: unit_nr, iostat
      LOGICAL                                              :: is_master, file_exists

      is_master = (bs_env%para_env%mepos == bs_env%para_env%source)
      WRITE(filename, '(A,I3.3,A)') "W_realspace_freq_", freq_idx, ".dat"

      SELECT CASE (TRIM(mode))
      CASE ("write", "put", "send")
         ! WRITE MODE: Save data to disk (training phase)
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="WRITE", file_form="UNFORMATTED", &
                           file_status="REPLACE")
            WRITE(unit_nr) MWM_R
            CALL close_file(unit_nr)
         END IF
         CALL bs_env%para_env%sync()

      CASE ("read_and_delete", "get_and_delete", "receive_and_delete")
         ! READ AND DELETE MODE: Load, delete data from disk and broadcast
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="READ", file_form="UNFORMATTED", &
                           file_status="OLD")
            READ(unit_nr) MWM_R
            CALL close_file(unit_nr)
         END IF
         
         ! Broadcast the large W_R to all cores
         CALL bs_env%para_env%bcast(MWM_R)
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="WRITE", file_form="UNFORMATTED", &
                           file_status="REPLACE")
            WRITE(unit_nr) MWM_R
            CALL close_file(unit_nr)

            CALL bs_env%para_env%bcast(MWM_R)

            
            ! Delete the file immediately after writing
            INQUIRE(FILE=filename, EXIST=file_exists)
            IF (file_exists) THEN
               OPEN(NEWUNIT=unit_nr, FILE=filename, STATUS='OLD', IOSTAT=iostat)
               IF (iostat == 0) THEN
                  CLOSE(unit_nr, STATUS='DELETE')
               END IF
            END IF
         END IF

         CALL bs_env%para_env%sync()

      CASE ("read", "get", "receive")
         ! READ MODE: Load data from disk and broadcast
         IF (is_master) THEN
            CALL open_file(file_name=filename, unit_number=unit_nr, &
                           file_action="READ", file_form="UNFORMATTED", &
                           file_status="OLD")
            READ(unit_nr) MWM_R
            CALL close_file(unit_nr)
         END IF
         
         ! Broadcast the large W_R to all cores
         CALL bs_env%para_env%bcast(MWM_R)

      CASE DEFAULT
         CPABORT("Invalid mode for W_storage_pipeline. Use 'write', 'read', or 'read_and_delete'.")
      END SELECT

   END SUBROUTINE W_storage_pipeline



END MODULE gw_lowrank_approx