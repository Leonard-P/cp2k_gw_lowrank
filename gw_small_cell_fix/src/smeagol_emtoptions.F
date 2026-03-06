!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K+SMEAGOL interface.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE smeagol_emtoptions
   USE cell_types, ONLY: cell_type, &
                         scaled_to_real
   USE cp_dbcsr_api, ONLY: dbcsr_get_info, &
                           dbcsr_type
#if defined(__SMEAGOL)
   USE global_meshvar, ONLY: smeagolglobal_orbital_BS => orbital_BS
#endif
   USE input_constants, ONLY: smeagol_gridmethod_adaptive, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_chebyshev, &
                              smeagol_integraltype_gauss_legendre
   USE kinds, ONLY: default_string_length, &
                    dp
#if defined(__SMEAGOL)
   USE mbfield, ONLY: smeagolglobal_ZeemanBx => ZeemanBx, &
                      smeagolglobal_ZeemanBy => ZeemanBy, &
                      smeagolglobal_ZeemanBz => ZeemanBz
#endif
   USE message_passing, ONLY: mp_para_env_type
#if defined(__SMEAGOL)
   USE mselfenergies, ONLY: SetOptionsSelfEnergies
   USE negfcoop, ONLY: coopinfo
   USE negfmod, ONLY: &
      smeagolglobal_dsigmade => dsigmade, &
      smeagolglobal_calcmpsh => calcmpsh, &
      smeagolglobal_NeneT => NeneT, &
      smeagolglobal_rmoinit => rmoinit, &
      smeagolglobal_rmolast => rmolast, &
      smeagolglobal_rorbinit => rorbinit, &
      smeagolglobal_rorblast => rorblast, &
      smeagolglobal_maxkappa => maxkappa, &
      smeagolglobal_TEnergI => TEnergI, &
      smeagolglobal_TEnergF => TEnergF, &
      smeagolglobal_vgate => vgate, &
      smeagolglobal_zgater => zgater, &
      smeagolglobal_zgatel => zgatel, &
      smeagolglobal_alphag => alphag, &
      smeagolglobal_timereversal => timereversal, &
      smeagolglobal_TRCAddVCDFT => TRCAddVCDFT, &
      smeagolglobal_AddVgCompensatingIsolatedLocalCharges => AddVgCompensatingIsolatedLocalCharges, &
      smeagolglobal_Add_Rho_Gate => Add_Rho_Gate, &
      smeagolglobal_debug_gate => debug_gate, &
      smeagolglobal_empdoskso => empdoskso, &
      smeagolglobal_emtimings => emtimings, &
      smeagolglobal_ik_select => ik_select, &
      smeagolglobal_Q_gate => Q_gate, &
      smeagolglobal_RGlxmin => RGlxmin, &
      smeagolglobal_RGlxmax => RGlxmax, &
      smeagolglobal_RGlymin => RGlymin, &
      smeagolglobal_RGlymax => RGlymax, &
      smeagolglobal_RGlzmin => RGlzmin, &
      smeagolglobal_RGlzmax => RGlzmax, &
      smeagolglobal_pdosgs => pdosgs, &
      smeagolglobal_skipright => skipright, &
      smeagolglobal_sigmatodisk => sigmatodisk, &
      smeagolglobal_bs_add => bs_add, &
      smeagolglobal_bs_writetrc => bs_writetrc, &
      smeagolglobal_writetk => writetk, &
      smeagolglobal_setene => setene, &
      smeagolglobal_overwritehs => overwritehs, &
      smeagolglobal_transmembed => transmembed, &
      smeagolglobal_writeevsigma => writeevsigma, &
      smeagolglobal_skiptransm => skiptransm, &
      smeagolglobal_deauto => deauto, &
      smeagolglobal_geterrsigma => geterrsigma, &
      smeagolglobal_getsigma => getsigma, &
      smeagolglobal_trcef => trcef, &
      smeagolglobal_RhoSetZeroIfHZero => RhoSetZeroIfHZero, &
      smeagolglobal_m_usehinv => m_usehinv, &
      smeagolglobal_m_usevinv => m_usevinv, &
      smeagolglobal_m_callsvd => m_callsvd, &
      smeagolglobal_m_complexbands => m_complexbands, &
      smeagolglobal_m_dosleads => m_dosleads, &
      smeagolglobal_leadsdos => leadsdos, &
      smeagolglobal_leadspdos => leadspdos, &
      smeagolglobal_emldos => emldos, &
      smeagolglobal_emldos2 => emldos2, &
      smeagolglobal_curr_distKEne => curr_distKEne, &
      smeagolglobal_curr_distK => curr_distK, &
      smeagolglobal_curr_dist => curr_dist, &
      smeagolglobal_curr_distTRC => curr_distTRC, &
      smeagolglobal_em_NonSCF => em_NonSCF, &
      smeagolglobal_TransmissionRL => TransmissionRL, &
      smeagolglobal_emSTTk => emSTTk, &
      smeagolglobal_GetT0S => GetT0S, &
      smeagolglobal_TransmissionChannels => TransmissionChannels, &
      smeagolglobal_TransmissionChannelsWFS => TransmissionChannelsWFS, &
      smeagolglobal_TransmissionMatrix => TransmissionMatrix, &
      smeagolglobal_TransmissionMatrixWFS => TransmissionMatrixWFS, &
      smeagolglobal_TransmissionMatrixPDOS => TransmissionMatrixPDOS, &
      smeagolglobal_EM_EndCode1 => EM_EndCode1, &
      smeagolglobal_TransmOverk => TransmOverk, &
      smeagolglobal_emdos => emdos, &
      smeagolglobal_empdos => empdos, &
      smeagolglobal_empdosk => empdosk, &
      smeagolglobal_WriteIkTrcK => WriteIkTrcK, &
      smeagolglobal_cbspart => cbspart, &
      smeagolglobal_evec => evec, &
      smeagolglobal_negfon => negfon, &
      smeagolglobal_outinfo => outinfo, &
      smeagolglobal_WriteNk => WriteNk, &
      smeagolglobal_Set_RhoBoundary_Leads => Set_RhoBoundary_Leads, &
      smeagolglobal_Set_RhoBoundaryOverlap_Leads => Set_RhoBoundaryOverlap_Leads, &
      smeagolglobal_Set_HBoundary_Leads => Set_HBoundary_Leads, &
      smeagolglobal_Set_HLR_Zero => Set_HLR_Zero, &
      smeagolglobal_SkipNonEquilibriumRho => SkipNonEquilibriumRho, &
      smeagolglobal_SkipEquilibriumRho => SkipEquilibriumRho, &
      smeagolglobal_ComputeImpurityGfMatsubara => ComputeImpurityGfMatsubara, &
      smeagolglobal_PrintImpurityGfMatsubara => PrintImpurityGfMatsubara, &
      smeagolglobal_CallImpuritySolver => CallImpuritySolver, &
      smeagolglobal_maxdepth => maxdepth, &
      smeagolglobal_ndivisions => ndivisions, &
      smeagolglobal_storesigma => storesigma, &
      smeagolglobal_bs_method => bs_method, &
      smeagolglobal_bs_nmid => bs_nmid, &
      smeagolglobal_NRunSigmaMax => NRunSigmaMax, &
      smeagolglobal_bsrun => bsrun, &
      smeagolglobal_bsskip => bsskip, &
      smeagolglobal_bseskip => bseskip, &
      smeagolglobal_bssc => bssc, &
      smeagolglobal_inversion_solver => inversion_solver, &
      smeagolglobal_onekp => onekp, &
      smeagolglobal_nbss => nbss, &
      smeagolglobal_nSplit => nSplit, &
      smeagolglobal_TransmissionChannelsWFSSkipEne => TransmissionChannelsWFSSkipEne, &
      smeagolglobal_TransmissionChannelsWFSSkipKP => TransmissionChannelsWFSSkipKP, &
      smeagolglobal_TransmissionMatrixPDOSNWrite => TransmissionMatrixPDOSNWrite, &
      smeagolglobal_TransmissionMatrixiSetPhase => TransmissionMatrixiSetPhase, &
      smeagolglobal_MinChannelIndex => MinChannelIndex, &
      smeagolglobal_MaxChannelIndex => MaxChannelIndex, &
      smeagolglobal_N_Replace_L => N_Replace_L, &
      smeagolglobal_N_Replace_R => N_Replace_R, &
      smeagolglobal_curr_nb1 => curr_nb1, &
      smeagolglobal_curr_nb2 => curr_nb2, &
      smeagolglobal_NParallelK => NParallelK, &
      smeagolglobal_GetRhoSingleLead => GetRhoSingleLead, &
      smeagolglobal_bs_tol => bs_tol, &
      smeagolglobal_bs_min => bs_min, &
      smeagolglobal_HSTol => HSTol, &
      smeagolglobal_deltaimag => deltaimag, &
      smeagolglobal_deltaini => deltaini, &
      smeagolglobal_deltatode => deltatode, &
      smeagolglobal_critam => critam, &
      smeagolglobal_setemin => setemin, &
      smeagolglobal_setemax => setemax, &
      smeagolglobal_tolab => tolab, &
      smeagolglobal_weightrho => weightrho, &
      smeagolglobal_ThetaLeadsL => ThetaLeadsL, &
      smeagolglobal_PhiLeadsL => PhiLeadsL, &
      smeagolglobal_ThetaLeadsR => ThetaLeadsR, &
      smeagolglobal_PhiLeadsR => PhiLeadsR, &
      smeagolglobal_deltaimagtrc => deltaimagtrc, &
      smeagolglobal_m_tolki => m_tolki, &
      smeagolglobal_m_svdtolmax => m_svdtolmax, &
      smeagolglobal_m_svdtolmin => m_svdtolmin, &
      smeagolglobal_m_dsigmamax => m_dsigmamax, &
      smeagolglobal_m_rnoise => m_rnoise, &
      smeagolglobal_m_skipsvd => m_skipsvd, &
      smeagolglobal_m_svdtolzi => m_svdtolzi, &
      smeagolglobal_SigmaWideBand => SigmaWideBand, &
      smeagolglobal_emSTTkScale => emSTTkScale, &
      smeagolglobal_nlSetZeroRatio => nlSetZeroRatio, &
      smeagolglobal_nrSetZeroRatio => nrSetZeroRatio, &
      smeagolglobal_curr_fl_L => curr_fl_L, &
      smeagolglobal_curr_fr_L => curr_fr_L, &
      smeagolglobal_curr_fl_R => curr_fl_R, &
      smeagolglobal_curr_fr_R => curr_fr_R, &
      smeagolglobal_deltamin => deltamin, &
      smeagolglobal_tolsigma => tolsigma, &
      smeagolglobal_trcde => trcde, &
      smeagolglobal_deltabss_bs => deltabss_bs, &
      smeagolglobal_gamma_negf => gamma_negf, &
      smeagolglobal_emforces => emforces, &
      smeagolglobal_emSTT => emSTT, &
      smeagolglobal_emSTTLin => emSTTLin, &
      smeagolglobal_gridmethod => gridmethod, &
      smeagolglobal_integraltype => integraltype, &
      smeagolglobal_ndivxy => ndivxy, &
      smeagolglobal_ndivxyNL => ndivxyNL, &
      smeagolglobal_negf_base_comm => negf_base_comm, &
      smeagolglobal_nebss_bs => nebss_bs, &
      smeagolglobal_nprocs_hs => nprocs_hs
   USE sigma, ONLY: &
      smeagolglobal_sigma_method => sigma_method
#endif
      USE smeagol_control_types, ONLY: smeagol_control_type
      USE string_utilities, ONLY: integer_to_string
#include "./base/base_uses.f90"

      IMPLICIT NONE
      PRIVATE

      CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_emtoptions'

      PUBLIC :: reademtr, ReadOptionsNEGF_DFT, emtrans_options, emtrans_deallocate_global_arrays

   CONTAINS

      SUBROUTINE reademtr(smeagol_control, natoms, gamma_negf)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         INTEGER, INTENT(in)                                :: natoms
         LOGICAL, INTENT(in)                                :: gamma_negf

         CHARACTER(LEN=*), PARAMETER :: routineN = 'reademtr'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         smeagolglobal_gamma_negf = gamma_negf

         IF (smeagol_control%aux%AtmRVCte <= 0) smeagol_control%aux%AtmRVCte = natoms

            smeagolglobal_dsigmade = smeagol_control%aux%dsigmade
            smeagolglobal_calcmpsh = smeagol_control%aux%calcmpsh

            smeagolglobal_NeneT = smeagol_control%aux%NeneT
            smeagolglobal_rmoinit = smeagol_control%aux%rmoinit
            smeagolglobal_rmolast = smeagol_control%aux%rmolast
            smeagolglobal_rorbinit = smeagol_control%aux%rorbinit
            smeagolglobal_rorblast = smeagol_control%aux%rorblast

            smeagolglobal_maxkappa = smeagol_control%aux%maxkappa

            smeagolglobal_TEnergI = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%TEnergI
            smeagolglobal_TEnergF = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%TEnergF
            smeagolglobal_vgate = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%vgate
            smeagolglobal_zgater = smeagol_control%aux%zgater
            smeagolglobal_zgatel = smeagol_control%aux%zgatel
            smeagolglobal_alphag = smeagol_control%aux%alphag

         coopinfo%ccoop = smeagol_control%aux%COOPCalculate
         coopinfo%nbond = smeagol_control%aux%COOPNumberOfBonds
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(natoms)
         MARK_USED(gamma_negf)
#endif

         CALL timestop(handle)
      END SUBROUTINE reademtr

      SUBROUTINE ReadOptionsNEGF_DFT(smeagol_control, ucell, torqueflag, torquelin)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         TYPE(cell_type), POINTER                           :: ucell
         LOGICAL, INTENT(in)                                :: torqueflag, torquelin

         CHARACTER(LEN=*), PARAMETER :: routineN = 'ReadOptionsNEGF_DFT'

         INTEGER                                            :: handle
         REAL(kind=dp), DIMENSION(3)                        :: coord_real, coord_scaled

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         smeagolglobal_emSTT = torqueflag
         smeagolglobal_emSTTLin = torquelin

         ! In case of the original SIESTA+SMEAGOL, 'TimeReversal' keyword is enabled by default, therefore 'EM.TimeReversal' is also enabled.
         ! In case of this CP2K+SMEAGOL interface, the default value of 'timereversal' variable is .FALSE.
         IF (smeagol_control%aux%timereversal) THEN
            CALL cp_warn(__LOCATION__, &
                         "The SMEAGOL keyword 'EM.TimeReversal' is not supported yet.")
         END IF

            smeagolglobal_timereversal = smeagol_control%aux%timereversal
            smeagolglobal_TRCAddVCDFT = smeagol_control%aux%TRCAddVCDFT
            smeagolglobal_AddVgCompensatingIsolatedLocalCharges = smeagol_control%aux%AddVgCompensatingIsolatedLocalCharges
            smeagolglobal_Add_Rho_Gate = smeagol_control%aux%Add_Rho_Gate
            smeagolglobal_debug_gate = smeagol_control%aux%debug_gate
            smeagolglobal_empdoskso = smeagol_control%aux%empdoskso
            smeagolglobal_emtimings = smeagol_control%aux%emtimings

            smeagolglobal_ik_select = smeagol_control%aux%ik_select

            smeagolglobal_Q_gate = smeagol_control%aux%Q_gate

         ! readoptsnegf_negfmod_plist
         IF (.NOT. smeagol_control%aux%isexplicit_RGlxmin) smeagol_control%aux%RGlxmin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlymin) smeagol_control%aux%RGlymin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlzmin) smeagol_control%aux%RGlzmin = 0.0_dp
         IF (.NOT. smeagol_control%aux%isexplicit_RGlxmax) THEN
            coord_scaled(:) = [1.0_dp, 0.0_dp, 0.0_dp]
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlxmax = coord_real(1)
         END IF
         IF (.NOT. smeagol_control%aux%isexplicit_RGlymax) THEN
            coord_scaled(:) = [0.0_dp, 1.0_dp, 0.0_dp]
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlymax = coord_real(2)
         END IF
         IF (.NOT. smeagol_control%aux%isexplicit_RGlzmax) THEN
            coord_scaled(:) = [0.0_dp, 0.0_dp, 1.0_dp]
            CALL scaled_to_real(coord_real, coord_scaled, ucell)
            smeagol_control%aux%RGlzmax = coord_real(3)
         END IF
            smeagolglobal_RGlxmin = smeagol_control%aux%RGlxmin
            smeagolglobal_RGlxmax = smeagol_control%aux%RGlxmax
            smeagolglobal_RGlymin = smeagol_control%aux%RGlymin
            smeagolglobal_RGlymax = smeagol_control%aux%RGlymax
            smeagolglobal_RGlzmin = smeagol_control%aux%RGlzmin
            smeagolglobal_RGlzmax = smeagol_control%aux%RGlzmax

         ! options to add a Zeeman term to the Hamiltonian
            smeagolglobal_ZeemanBx = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%ZeemanBx
            smeagolglobal_ZeemanBy = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%ZeemanBy
            smeagolglobal_ZeemanBz = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%ZeemanBz
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(ucell)
         MARK_USED(torqueflag)
         MARK_USED(torquelin)
         ! local variables
         MARK_USED(coord_real)
         MARK_USED(coord_scaled)
#endif

         CALL timestop(handle)
      END SUBROUTINE ReadOptionsNEGF_DFT

      SUBROUTINE emtrans_options(smeagol_control, matrix_s, para_env, iter, istep, inicoor, iv, delta, nk)
         TYPE(smeagol_control_type), POINTER                :: smeagol_control
         TYPE(dbcsr_type), INTENT(in), POINTER              :: matrix_s
         TYPE(mp_para_env_type), POINTER                    :: para_env
         INTEGER, INTENT(in)                                :: iter, istep, inicoor, iv
         REAL(kind=dp), INTENT(in)                          :: delta
         INTEGER, INTENT(in)                                :: nk

         CHARACTER(LEN=*), PARAMETER :: routineN = 'emtrans_options'

         CHARACTER(len=default_string_length)               :: actual_val_str, expected_val_str
         INTEGER                                            :: GetRhoSingleLeadDefault, handle, i, iatom, n1, nblkcols_total, &
                                                               NParallelK, NParallelKbuf
         INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size

         CALL timeset(routineN, handle)

#if defined(__SMEAGOL)
         CPASSERT(ASSOCIATED(smeagol_control%aux))

         CALL dbcsr_get_info(matrix=matrix_s, nblkcols_total=nblkcols_total, &
                             col_blk_size=col_blk_size, col_blk_offset=col_blk_offset)
         ! number of atomic orbitals.
         ! The name n1 is meaningless, but it is used in the original SIESTA's version of emtrans_options() subroutine
         n1 = SUM(col_blk_size(1:nblkcols_total))

         IF ((iter == 1) .AND. (istep == inicoor) .AND. (iv == 0)) THEN

            IF (smeagol_control%aux%gridmethod == smeagol_gridmethod_traditional) THEN
               smeagolglobal_gridmethod = 'Traditional'
            ELSE IF (smeagol_control%aux%gridmethod == smeagol_gridmethod_adaptive) THEN
               smeagolglobal_gridmethod = 'Adaptivegrid'
            ELSE
               smeagolglobal_gridmethod = 'UNKNOWN'
            END IF

            IF (smeagol_control%aux%integraltype == smeagol_integraltype_gauss_legendre) THEN
               smeagolglobal_integraltype = 'gauss-legendre'
            ELSE IF (smeagol_control%aux%integraltype == smeagol_integraltype_gauss_chebyshev) THEN
               smeagolglobal_integraltype = 'gauss-chebyshev'
            ELSE
               smeagolglobal_integraltype = 'UNKNOWN'
            END IF

            smeagolglobal_negf_base_comm = para_env%get_handle()

            IF (MOD(smeagol_control%aux%ndivisions, 2) == 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "AM.NumberDivisions value must be odd.")
            END IF

            ! It seems that 'sigmatodisk' parameter is specific to the SIESTA interface.
            ! Consider setting smeagol_control%aux%sigmatodisk = .FALSE. and remove smeagol_control%aux%storesigma input keyword.
            IF (smeagol_control%aux%storesigma == 2) THEN
               smeagol_control%aux%sigmatodisk = .TRUE.
            ELSE IF (smeagol_control%aux%storesigma == 1) THEN
               smeagol_control%aux%sigmatodisk = .FALSE.
            END IF

            ! Bound states
            !     options:
            !       bs_add    : true => bound states are added
            !                   false=> normal smeagol
            !       bs_method : 0    => calculate bound states with effective Hamiltonian
            !                   1    => calculate bound states by adding a small
            !                           imaginary part to the selfenergies
            IF (smeagol_control%aux%bs_nmid == 0) smeagol_control%aux%bs_nmid = n1/2

            ! SC: Bad practice: real number comparison
            IF (smeagol_control%aux%SigmaWideBand /= 0.0_dp) smeagol_control%aux%m_svdtolzi = 0.0_dp

            IF (smeagol_control%aux%leadspdos) smeagol_control%aux%leadsdos = .TRUE.
            IF (smeagol_control%aux%curr_distKEne) smeagol_control%aux%curr_distK = .TRUE.
            IF (smeagol_control%aux%curr_distK) smeagol_control%aux%curr_dist = .TRUE.

            IF (smeagolglobal_emSTT .AND. smeagolglobal_emSTTLin .OR. smeagol_control%aux%curr_dist) THEN
               smeagol_control%aux%emldos2 = .TRUE.
               IF (smeagol_control%aux%curr_dist) THEN
                  GetRhoSingleLeadDefault = 3
               ELSE
                  GetRhoSingleLeadDefault = -3
               END IF
            ELSE
               GetRhoSingleLeadDefault = 0
            END IF

            ! current-induced forces
            ! The value of 'smeagol_control%emforces' is set in qs_energies().
            ! Calculation of forces is enabled automatically for certain run_types (energy_force, geo_opt, md) and disabled otherwise.
            IF (smeagol_control%aux%curr_dist) THEN
               smeagol_control%emforces = .TRUE.
            END IF

            IF (.NOT. smeagol_control%aux%isexplicit_nprocs_hs) smeagol_control%aux%nprocs_hs = smeagol_control%aux%nprocs_inverse
            smeagolglobal_nprocs_hs = smeagol_control%aux%nprocs_hs
           IF (.NOT. smeagol_control%aux%isexplicit_GetRhoSingleLead) smeagol_control%aux%GetRhoSingleLead = GetRhoSingleLeadDefault

            IF (smeagol_control%aux%MinChannelIndex < 1) smeagol_control%aux%MinChannelIndex = 1
            IF (smeagol_control%aux%MaxChannelIndex < 1) &
               smeagol_control%aux%MaxChannelIndex = smeagol_control%aux%MinChannelIndex + 4

            IF (smeagolglobal_emSTT .AND. smeagolglobal_emSTTLin .AND. smeagol_control%aux%GetRhoSingleLead /= -3) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.LDOSLeadsProjection should be set to -3. "// &
                            "If SpinTorque and STLinResp are T, otherwise the output "// &
                            "results for the spin transfer torque are incorrect.")
            END IF

            ! NParallelK
            NParallelK = smeagol_control%aux%NParallelK
            IF (MOD(para_env%num_pe, NParallelK) /= 0) then
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum must be a divisor of the total number of "// &
                            "MPI processes used in a run; resetting the value of NParallelK.")
               NParallelKbuf = NParallelK
               DO NParallelK = NParallelKbuf, 1, -1
                  IF (MOD(para_env%num_pe, NParallelK) == 0) EXIT
               END DO
            END IF

            IF (NParallelK > para_env%num_pe) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum can not be larger than the total number of "// &
                            "MPI processes used in a run; resetting the value of NParallelK.")
               NParallelK = -1
            END IF

            IF (NParallelK > nk) THEN
               CALL cp_warn(__LOCATION__, &
                            "EM.ParallelOverKNum can not be larger than the total number of "// &
                            "k-points used in a run; resetting the value of NParallelK.")
               NParallelK = -1
            END IF

            IF (NParallelK == -1) THEN
               DO NParallelK = nk, 1, -1
                  IF (MOD(para_env%num_pe, NParallelK) == 0) EXIT
               END DO
            END IF

            smeagol_control%aux%NParallelK = NParallelK

            IF (smeagol_control%aux%empdosk) smeagol_control%aux%empdos = .TRUE.
            IF (smeagol_control%aux%emldos2) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%TransmissionChannels) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%TransmissionMatrix) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%curr_dist) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%empdos) smeagol_control%aux%emdos = .TRUE.
            IF (smeagol_control%aux%m_skipsvd < 2.0_dp) smeagol_control%aux%m_skipsvd = 10.0_dp

            IF (smeagol_control%aux%CallImpuritySolver) smeagol_control%aux%ComputeImpurityGfMatsubara = .TRUE.

               smeagolglobal_pdosgs = smeagol_control%aux%pdosgs
               smeagolglobal_skipright = smeagol_control%aux%skipright
               smeagolglobal_sigmatodisk = smeagol_control%aux%sigmatodisk
               smeagolglobal_bs_add = smeagol_control%aux%bs_add
               smeagolglobal_bs_writetrc = smeagol_control%aux%bs_writetrc
               smeagolglobal_writetk = smeagol_control%aux%writetk
               smeagolglobal_setene = smeagol_control%aux%setene
               smeagolglobal_overwritehs = smeagol_control%aux%overwritehs
               smeagolglobal_transmembed = smeagol_control%aux%transmembed
               smeagolglobal_writeevsigma = smeagol_control%aux%writeevsigma
               smeagolglobal_skiptransm = smeagol_control%aux%skiptransm
               smeagolglobal_deauto = smeagol_control%aux%deauto
               smeagolglobal_geterrsigma = smeagol_control%aux%geterrsigma
               smeagolglobal_getsigma = smeagol_control%aux%getsigma
               smeagolglobal_trcef = smeagol_control%aux%trcef
               smeagolglobal_RhoSetZeroIfHZero = smeagol_control%aux%RhoSetZeroIfHZero
               smeagolglobal_m_usehinv = smeagol_control%aux%m_usehinv
               smeagolglobal_m_usevinv = smeagol_control%aux%m_usevinv
               smeagolglobal_m_callsvd = smeagol_control%aux%m_callsvd
               smeagolglobal_m_complexbands = smeagol_control%aux%m_complexbands
               smeagolglobal_m_dosleads = smeagol_control%aux%m_dosleads
               smeagolglobal_leadsdos = smeagol_control%aux%leadsdos
               smeagolglobal_leadspdos = smeagol_control%aux%leadspdos
               smeagolglobal_emldos = smeagol_control%aux%emldos
               smeagolglobal_emldos2 = smeagol_control%aux%emldos2
               smeagolglobal_curr_distKEne = smeagol_control%aux%curr_distKEne
               smeagolglobal_curr_distK = smeagol_control%aux%curr_distK
               smeagolglobal_curr_dist = smeagol_control%aux%curr_dist
               smeagolglobal_curr_distTRC = smeagol_control%aux%curr_distTRC
               smeagolglobal_em_NonSCF = smeagol_control%aux%em_NonSCF
               smeagolglobal_TransmissionRL = smeagol_control%aux%TransmissionRL
               smeagolglobal_emSTTk = smeagol_control%aux%emSTTk
               smeagolglobal_GetT0S = smeagol_control%aux%GetT0S
               smeagolglobal_TransmissionChannels = smeagol_control%aux%TransmissionChannels
               smeagolglobal_TransmissionChannelsWFS = smeagol_control%aux%TransmissionChannelsWFS
               smeagolglobal_TransmissionMatrix = smeagol_control%aux%TransmissionMatrix
               smeagolglobal_TransmissionMatrixWFS = smeagol_control%aux%TransmissionMatrixWFS
               smeagolglobal_TransmissionMatrixPDOS = smeagol_control%aux%TransmissionMatrixPDOS
               smeagolglobal_EM_EndCode1 = smeagol_control%aux%EM_EndCode1
               smeagolglobal_TransmOverk = smeagol_control%aux%TransmOverk
               smeagolglobal_emdos = smeagol_control%aux%emdos
               smeagolglobal_empdos = smeagol_control%aux%empdos
               smeagolglobal_empdosk = smeagol_control%aux%empdosk
               smeagolglobal_WriteIkTrcK = smeagol_control%aux%WriteIkTrcK
               smeagolglobal_cbspart = smeagol_control%aux%cbspart
               smeagolglobal_evec = smeagol_control%aux%evec
               smeagolglobal_negfon = smeagol_control%aux%negfon
               smeagolglobal_outinfo = smeagol_control%aux%outinfo
               smeagolglobal_WriteNk = smeagol_control%aux%WriteNk
               smeagolglobal_Set_RhoBoundary_Leads = smeagol_control%aux%Set_RhoBoundary_Leads
               smeagolglobal_Set_RhoBoundaryOverlap_Leads = smeagol_control%aux%Set_RhoBoundaryOverlap_Leads
               smeagolglobal_Set_HBoundary_Leads = smeagol_control%aux%Set_HBoundary_Leads
               smeagolglobal_Set_HLR_Zero = smeagol_control%aux%Set_HLR_Zero
               smeagolglobal_SkipNonEquilibriumRho = smeagol_control%aux%SkipNonEquilibriumRho
               smeagolglobal_SkipEquilibriumRho = smeagol_control%aux%SkipEquilibriumRho
               smeagolglobal_ComputeImpurityGfMatsubara = smeagol_control%aux%ComputeImpurityGfMatsubara
               smeagolglobal_PrintImpurityGfMatsubara = smeagol_control%aux%PrintImpurityGfMatsubara
               smeagolglobal_CallImpuritySolver = smeagol_control%aux%CallImpuritySolver
            smeagolglobal_emforces = smeagol_control%emforces

            ! emtoptions_local_ilist, emtoptions_local_explicit_ilist
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NxLeft) smeagol_control%aux%Sigma_NxLeft = smeagol_control%aux%Sigma_Nx
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NyLeft) smeagol_control%aux%Sigma_NyLeft = smeagol_control%aux%Sigma_Ny
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NxRight) smeagol_control%aux%Sigma_NxRight = smeagol_control%aux%Sigma_Nx
            IF (.NOT. smeagol_control%aux%isexplicit_Sigma_NyRight) smeagol_control%aux%Sigma_NyRight = smeagol_control%aux%Sigma_Ny

            smeagolglobal_ndivxy(1) = smeagol_control%aux%Sigma_Nx
            smeagolglobal_ndivxy(2) = smeagol_control%aux%Sigma_Ny
            smeagolglobal_ndivxyNL(1, 1) = smeagol_control%aux%Sigma_NxLeft
            smeagolglobal_ndivxyNL(1, 2) = smeagol_control%aux%Sigma_NyLeft
            smeagolglobal_ndivxyNL(2, 1) = smeagol_control%aux%Sigma_NxRight
            smeagolglobal_ndivxyNL(2, 2) = smeagol_control%aux%Sigma_NyRight

            CALL SetOptionsSelfEnergies(smeagolglobal_ndivxyNL, 2)

               smeagolglobal_maxdepth = smeagol_control%aux%maxdepth
               smeagolglobal_ndivisions = smeagol_control%aux%ndivisions
               smeagolglobal_storesigma = smeagol_control%aux%storesigma
               smeagolglobal_bs_method = smeagol_control%aux%bs_method
               smeagolglobal_bs_nmid = smeagol_control%aux%bs_nmid
               smeagolglobal_NRunSigmaMax = smeagol_control%aux%NRunSigmaMax
               smeagolglobal_bsrun = smeagol_control%aux%bsrun
               smeagolglobal_bsskip = smeagol_control%aux%bsskip
               smeagolglobal_bseskip = smeagol_control%aux%bseskip
               smeagolglobal_bssc = smeagol_control%aux%bssc
               smeagolglobal_inversion_solver = smeagol_control%aux%inversion_solver
               smeagolglobal_onekp = smeagol_control%aux%onekp
               smeagolglobal_nbss = smeagol_control%aux%nbss
               smeagolglobal_nSplit = smeagol_control%aux%nSplit
               smeagolglobal_TransmissionChannelsWFSSkipEne = smeagol_control%aux%TransmissionChannelsWFSSkipEne
               smeagolglobal_TransmissionChannelsWFSSkipKP = smeagol_control%aux%TransmissionChannelsWFSSkipKP
               smeagolglobal_TransmissionMatrixPDOSNWrite = smeagol_control%aux%TransmissionMatrixPDOSNWrite
               smeagolglobal_TransmissionMatrixiSetPhase = smeagol_control%aux%TransmissionMatrixiSetPhase
               smeagolglobal_MinChannelIndex = smeagol_control%aux%MinChannelIndex
               smeagolglobal_MaxChannelIndex = smeagol_control%aux%MaxChannelIndex
               smeagolglobal_N_Replace_L = smeagol_control%aux%N_Replace_L
               smeagolglobal_N_Replace_R = smeagol_control%aux%N_Replace_R
               smeagolglobal_curr_nb1 = smeagol_control%aux%curr_nb1
               smeagolglobal_curr_nb2 = smeagol_control%aux%curr_nb2
               smeagolglobal_NParallelK = smeagol_control%aux%NParallelK

               smeagolglobal_GetRhoSingleLead = smeagol_control%aux%GetRhoSingleLead

               smeagolglobal_sigma_method = smeagol_control%aux%sigma_method

               smeagolglobal_bs_tol = smeagol_control%aux%bs_tol
               smeagolglobal_bs_min = smeagol_control%aux%bs_min
               smeagolglobal_HSTol = smeagol_control%aux%HSTol
               smeagolglobal_deltaimag = smeagol_control%aux%deltaimag
               smeagolglobal_deltaini = smeagol_control%aux%deltaini
               smeagolglobal_deltatode = smeagol_control%aux%deltatode
               smeagolglobal_critam = smeagol_control%aux%critam
               smeagolglobal_setemin = smeagol_control%aux%setemin
               smeagolglobal_setemax = smeagol_control%aux%setemax
               smeagolglobal_tolab = smeagol_control%aux%tolab
               smeagolglobal_weightrho = smeagol_control%aux%weightrho
               smeagolglobal_ThetaLeadsL = smeagol_control%aux%ThetaLeadsL
               smeagolglobal_PhiLeadsL = smeagol_control%aux%PhiLeadsL
               smeagolglobal_ThetaLeadsR = smeagol_control%aux%ThetaLeadsR
               smeagolglobal_PhiLeadsR = smeagol_control%aux%PhiLeadsR
               smeagolglobal_deltaimagtrc = smeagol_control%aux%deltaimagtrc
               smeagolglobal_m_tolki = smeagol_control%aux%m_tolki
               smeagolglobal_m_svdtolmax = smeagol_control%aux%m_svdtolmax
               smeagolglobal_m_svdtolmin = smeagol_control%aux%m_svdtolmin
               smeagolglobal_m_dsigmamax = smeagol_control%aux%m_dsigmamax
               smeagolglobal_m_rnoise = smeagol_control%aux%m_rnoise
               smeagolglobal_m_skipsvd = smeagol_control%aux%m_skipsvd
               smeagolglobal_m_svdtolzi = smeagol_control%aux%m_svdtolzi
               smeagolglobal_SigmaWideBand = smeagol_control%aux%SigmaWideBand
               smeagolglobal_emSTTkScale = smeagol_control%aux%emSTTkScale
               smeagolglobal_nlSetZeroRatio = smeagol_control%aux%nlSetZeroRatio
               smeagolglobal_nrSetZeroRatio = smeagol_control%aux%nrSetZeroRatio
               smeagolglobal_curr_fl_L = smeagol_control%aux%curr_fl_L
               smeagolglobal_curr_fr_L = smeagol_control%aux%curr_fr_L
               smeagolglobal_curr_fl_R = smeagol_control%aux%curr_fl_R
               smeagolglobal_curr_fr_R = smeagol_control%aux%curr_fr_R

            ! +++ bound states
            IF (.NOT. smeagol_control%aux%isexplicit_deltamin) smeagol_control%aux%deltamin = delta

               smeagolglobal_deltamin = smeagol_control%aux%deltamin

               smeagolglobal_tolsigma = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%tolsigma
               smeagolglobal_trcde = smeagol_control%to_smeagol_energy_units*smeagol_control%aux%trcde

            ! At present read_options_ImpuritySolver() is a dummy libsmeagol.a subroutine,
            ! so we do not call it
            !IF (smeagol_control%ComputeImpurityGfMatsubara) THEN
            !   CALL read_options_ImpuritySolver(CallImpuritySolver,n1)
            !END IF
         END IF

         ! *** bound states (re-implemented ReadSpeciesBS() and ReadBSSubSystemBoundaries() subroutines)

         ! bound- state-related global allocatable arrays
         ALLOCATE (smeagolglobal_deltabss_bs(smeagol_control%aux%nbss))
         ALLOCATE (smeagolglobal_nebss_bs(smeagol_control%aux%nbss, 2))
         ALLOCATE (smeagolglobal_orbital_BS(n1))

         ! BS.SubSystemsDelta
         IF (ALLOCATED(smeagol_control%aux%deltabss_bs)) THEN
            smeagolglobal_deltabss_bs(:) = smeagol_control%aux%deltabss_bs(:)
         ELSE
            smeagolglobal_deltabss_bs(:) = smeagol_control%aux%deltamin
         END IF

         ! BS.SubSystemsBoundaries
         IF (ALLOCATED(smeagol_control%aux%nebss_bs)) THEN
            IF (MAXVAL(smeagol_control%aux%nebss_bs) > n1) THEN
               CALL integer_to_string(MAXVAL(smeagol_control%aux%nebss_bs), actual_val_str)
               CALL integer_to_string(n1, expected_val_str)
               CALL cp_abort(__LOCATION__, &
                             "The largest index in BS.SubSystemsBoundaries section ("//TRIM(actual_val_str)// &
                             ") exceeds the number of atomic orbitals ("//TRIM(expected_val_str)//").")
            END IF
            smeagolglobal_nebss_bs(:, :) = smeagol_control%aux%nebss_bs(:, :)
         ELSE
            smeagolglobal_nebss_bs(1, 1) = 1
            smeagolglobal_nebss_bs(1, 2) = n1/smeagol_control%aux%nbss
            DO i = 2, smeagol_control%aux%nbss
               smeagolglobal_nebss_bs(i, 1) = smeagolglobal_nebss_bs(i - 1, 2) + 1
               smeagolglobal_nebss_bs(i, 2) = i*n1/smeagol_control%aux%nbss
            END DO
            smeagolglobal_nebss_bs(smeagol_control%aux%nbss, 2) = n1
         END IF

         ! AM.AtomListBS
         IF (ALLOCATED(smeagol_control%aux%atomlist_bs)) THEN
            IF (MAXVAL(smeagol_control%aux%atomlist_bs) > nblkcols_total) THEN
               CALL integer_to_string(MAXVAL(smeagol_control%aux%atomlist_bs), actual_val_str)
               CALL integer_to_string(nblkcols_total, expected_val_str)
               CALL cp_abort(__LOCATION__, &
                             "The largest atomic index in AM.AtomListBS keyword ("//TRIM(actual_val_str)// &
                             ") exceeds the number of atoms ("//TRIM(expected_val_str)//") in the extended molecule.")
            END IF

            smeagolglobal_orbital_BS(:) = .FALSE.
            DO i = 1, SIZE(smeagol_control%aux%atomlist_bs)
               iatom = smeagol_control%aux%atomlist_bs(i)
               smeagolglobal_orbital_BS(col_blk_offset(iatom):col_blk_offset(iatom) + col_blk_size(iatom) - 1) = .TRUE.
            END DO
         ELSE
            smeagolglobal_orbital_BS(:) = .TRUE.
         END IF
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
         MARK_USED(smeagol_control)
         MARK_USED(matrix_s)
         MARK_USED(para_env)
         MARK_USED(iter)
         MARK_USED(istep)
         MARK_USED(inicoor)
         MARK_USED(iv)
         MARK_USED(delta)
         MARK_USED(nk)
         ! local variables
         MARK_USED(actual_val_str)
         MARK_USED(col_blk_offset)
         MARK_USED(col_blk_size)
         MARK_USED(expected_val_str)
         MARK_USED(GetRhoSingleLeadDefault)
         MARK_USED(i)
         MARK_USED(iatom)
         MARK_USED(n1)
         MARK_USED(nblkcols_total)
         MARK_USED(NParallelK)
         MARK_USED(NParallelKbuf)
#endif

         CALL timestop(handle)
      END SUBROUTINE emtrans_options

      SUBROUTINE emtrans_deallocate_global_arrays()

#if defined(__SMEAGOL)
         IF (ALLOCATED(smeagolglobal_deltabss_bs)) DEALLOCATE (smeagolglobal_deltabss_bs)
         IF (ALLOCATED(smeagolglobal_nebss_bs)) DEALLOCATE (smeagolglobal_nebss_bs)
         IF (ALLOCATED(smeagolglobal_orbital_BS)) DEALLOCATE (smeagolglobal_orbital_BS)
#endif

      END SUBROUTINE emtrans_deallocate_global_arrays

   END MODULE smeagol_emtoptions

