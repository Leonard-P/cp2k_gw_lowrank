!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input control types for NEGF/SMEAGOL transport calculations.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************

MODULE smeagol_control_types
   USE cp_units, ONLY: cp_unit_from_cp2k
   USE input_constants, ONLY: smeagol_bulklead_leftright, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_legendre, &
                              smeagol_runtype_bulktransport, &
                              smeagol_runtype_emtransport
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: default_string_length, &
                    dp
   USE physcon, ONLY: kelvin
   USE string_utilities, ONLY: integer_to_string
   USE util, ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_control_types'

   PUBLIC :: smeagol_control_type
   PUBLIC :: smeagol_control_create, smeagol_control_release, read_smeagol_control

! **************************************************************************************************
!> \brief SMEAGOL-related input parameters
! **************************************************************************************************
   TYPE smeagol_control_type
      LOGICAL                                            :: smeagol_enabled = .FALSE.

      !> type of calculation
      INTEGER                                            :: run_type = smeagol_runtype_bulktransport

      !> regression test mode
      LOGICAL                                            :: do_regtest = .FALSE.

      !> current-induced forces. It is set automatically based on GLOBAL/run_type
      LOGICAL                                            :: emforces = .FALSE.

      !> scale factor to convert from CP2K (Hartree) to SMEAGOL (Rydberg) default energy unit
      REAL(kind=dp)                                      :: to_smeagol_energy_units = 2.0_dp

      !> number of cell images along i and j cell vectors
      INTEGER, DIMENSION(2)                              :: n_cell_images = (/1, 1/)

      !> what lead (bulk transport calculation)
      INTEGER                                            :: lead_label = smeagol_bulklead_leftright

      !> The length of the SMEAGOL project name is limited by 20 characters (hardcoded in SMEAGOL)
      CHARACTER(len=20)                                  :: project_name = "PROJECT_NAME"

      TYPE(smeagol_aux_control_type), POINTER            :: aux => NULL()

   END TYPE smeagol_control_type

! **************************************************************************************************
!> \brief SMEAGOL-related auxiliary input parameters. They remain unallocated when
!>        SMEAGOL support is disabled.
! **************************************************************************************************
   TYPE smeagol_aux_control_type
      !> MD or GEO_OPT iteration. In contrast with other components of smeagol_control_type that are read from the input file,
      !> these variables are initialised at the first NEGF iteration.
      !>
      !> index of MD or GEO_OPT iteration level.
      !> 0 if there is neither MD nor GEO_OPT iteration level (e.g. single-point energy calculation).
      !> -1 if it is not initialised by run_smeagol_emtrans() subroutine
      INTEGER                                            :: md_iter_level = -1
      !> The starting step value for MD / GEO_OPT iterations. The default value 0 can be overrited via STEP_START_VAL input keyword.
      INTEGER                                            :: md_first_step = 0

      !> BS.SubSystemsDelta(1:BS.Subsystems)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: deltabss_bs
      !> BS.SubSystemsBoundaries(1:BS.Subsystems, 1:2)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nebss_bs
      !>
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_bs

      REAL(kind=dp)                                      :: temperature = 300.0_dp/kelvin

      ! reademtr()
         LOGICAL                                            :: TrCoeff = .false.
         LOGICAL                                            :: skiplast = .false.
         LOGICAL                                            :: CalcIETS = .false.

         LOGICAL                                            :: isexplicit_maxR = .FALSE.

         LOGICAL                                            :: dsigmade = .false.
         LOGICAL                                            :: calcmpsh = .false.

         LOGICAL                                            :: COOPCalculate = .false.

         INTEGER                                            :: NEnergR = 288
         INTEGER                                            :: NEnergIC = 32
         INTEGER                                            :: NEnergIL = 16
         INTEGER                                            :: NPoles = 16
         INTEGER                                            :: SpinCL = 0
         INTEGER                                            :: AtmLVCte = 1
         INTEGER                                            :: AtmRVCte = -1
         INTEGER                                            :: NSlices = 1
         INTEGER                                            :: nlinmix = 1
         INTEGER                                            :: tmdskip = 0
         INTEGER                                            :: tmdsampling = 1

         INTEGER                                            :: NeneT = 100
         INTEGER                                            :: rmoinit = 1
         INTEGER                                            :: rmolast = 1
         INTEGER                                            :: rorbinit = 1
         INTEGER                                            :: rorblast = 1

         INTEGER                                            :: COOPNumberOfBonds = 0

         REAL(kind=dp)                                      :: Delta = 1.0e-10_dp

         REAL(kind=dp)                                      :: maxkappa = 0.0_dp

         REAL(kind=dp)                                      :: VBias = 0.0_dp*0.5_dp
         REAL(kind=dp)                                      :: EnergLB = -7.0_dp*0.5_dp
         REAL(kind=dp)                                      :: minL = 0.0_dp

         REAL(kind=dp)                                      :: maxR = 0.0_dp

         REAL(kind=dp)                                      :: TEnergI = -0.3675_dp*0.5_dp
         REAL(kind=dp)                                      :: TEnergF = -0.3675_dp*0.5_dp
         REAL(kind=dp)                                      :: vgate = 0.0_dp*0.5_dp

         REAL(kind=dp)                                      :: zgater = 0.0_dp
         REAL(kind=dp)                                      :: zgatel = 0.0_dp
         REAL(kind=dp)                                      :: alphag = 1.0_dp

      ! ReadOptionsNEGF_DFT()
         LOGICAL                                            :: timereversal = .false.
         LOGICAL                                            :: TRCAddVCDFT = .false.
         LOGICAL                                            :: AddVgCompensatingIsolatedLocalCharges = .true.
         LOGICAL                                            :: Add_Rho_Gate = .false.
         LOGICAL                                            :: debug_gate = .false.
         LOGICAL                                            :: empdoskso = .false.
         LOGICAL                                            :: emtimings = .false.

         LOGICAL                                            :: isexplicit_RGlxmin = .FALSE.
         LOGICAL                                            :: isexplicit_RGlxmax = .FALSE.
         LOGICAL                                            :: isexplicit_RGlymin = .FALSE.
         LOGICAL                                            :: isexplicit_RGlymax = .FALSE.
         LOGICAL                                            :: isexplicit_RGlzmin = .FALSE.
         LOGICAL                                            :: isexplicit_RGlzmax = .FALSE.

         INTEGER                                            :: ik_select = 0

         REAL(kind=dp)                                      :: Q_gate = 0.0_dp

         REAL(kind=dp)                                      :: RGlxmin = 0.0_dp
         REAL(kind=dp)                                      :: RGlxmax = 0.0_dp
         REAL(kind=dp)                                      :: RGlymin = 0.0_dp
         REAL(kind=dp)                                      :: RGlymax = 0.0_dp
         REAL(kind=dp)                                      :: RGlzmin = 0.0_dp
         REAL(kind=dp)                                      :: RGlzmax = 0.0_dp

         REAL(kind=dp)                                      :: ZeemanBx = 0.0_dp*0.5_dp
         REAL(kind=dp)                                      :: ZeemanBy = 0.0_dp*0.5_dp
         REAL(kind=dp)                                      :: ZeemanBz = 0.0_dp*0.5_dp

      ! emtrans_options()
      INTEGER                                            :: gridmethod = smeagol_gridmethod_traditional
      INTEGER                                            :: integraltype = smeagol_integraltype_gauss_legendre

         LOGICAL                                            :: pdosgs = .false.
         LOGICAL                                            :: skipright = .false.
         LOGICAL                                            :: sigmatodisk = .false.
         LOGICAL                                            :: bs_add = .false.
         LOGICAL                                            :: bs_writetrc = .false.
         LOGICAL                                            :: writetk = .false.
         LOGICAL                                            :: setene = .false.
         LOGICAL                                            :: overwritehs = .true.
         LOGICAL                                            :: transmembed = .true.
         LOGICAL                                            :: writeevsigma = .false.
         LOGICAL                                            :: skiptransm = .false.
         LOGICAL                                            :: deauto = .false.
         LOGICAL                                            :: geterrsigma = .false.
         LOGICAL                                            :: getsigma = .false.
         LOGICAL                                            :: trcef = .true.
         LOGICAL                                            :: RhoSetZeroIfHZero = .false.
         LOGICAL                                            :: m_usehinv = .false.
         LOGICAL                                            :: m_usevinv = .false.
         LOGICAL                                            :: m_callsvd = .false.
         LOGICAL                                            :: m_complexbands = .false.
         LOGICAL                                            :: m_dosleads = .false.
         LOGICAL                                            :: leadsdos = .false.
         LOGICAL                                            :: leadspdos = .false.
         LOGICAL                                            :: emldos = .false.
         LOGICAL                                            :: emldos2 = .false.
         LOGICAL                                            :: curr_distKEne = .false.
         LOGICAL                                            :: curr_distK = .false.
         LOGICAL                                            :: curr_dist = .false.
         LOGICAL                                            :: curr_distTRC = .false.
         LOGICAL                                            :: em_NonSCF = .false.
         LOGICAL                                            :: TransmissionRL = .false.
         LOGICAL                                            :: emSTTk = .false.
         LOGICAL                                            :: GetT0S = .false.
         LOGICAL                                            :: TransmissionChannels = .false.
         LOGICAL                                            :: TransmissionChannelsWFS = .false.
         LOGICAL                                            :: TransmissionMatrix = .false.
         LOGICAL                                            :: TransmissionMatrixWFS = .false.
         LOGICAL                                            :: TransmissionMatrixPDOS = .false.
         LOGICAL                                            :: EM_EndCode1 = .false.
         LOGICAL                                            :: TransmOverk = .false.
         LOGICAL                                            :: emdos = .false.
         LOGICAL                                            :: empdos = .false.
         LOGICAL                                            :: empdosk = .false.
         LOGICAL                                            :: WriteIkTrcK = .false.
         LOGICAL                                            :: cbspart = .false.
         LOGICAL                                            :: evec = .false.
         LOGICAL                                            :: negfon = .false.
         LOGICAL                                            :: outinfo = .false.
         LOGICAL                                            :: WriteNk = .false.
         LOGICAL                                            :: Set_RhoBoundary_Leads = .true.
         LOGICAL                                            :: Set_RhoBoundaryOverlap_Leads = .true.
         LOGICAL                                            :: Set_HBoundary_Leads = .true.
         LOGICAL                                            :: Set_HLR_Zero = .true.
         LOGICAL                                            :: SkipNonEquilibriumRho = .false.
         LOGICAL                                            :: SkipEquilibriumRho = .false.
         LOGICAL                                            :: ComputeImpurityGfMatsubara = .false.
         LOGICAL                                            :: PrintImpurityGfMatsubara = .false.
         LOGICAL                                            :: CallImpuritySolver = .false.

         LOGICAL                                            :: isexplicit_nprocs_hs = .FALSE.
         LOGICAL                                            :: isexplicit_Sigma_NxLeft = .FALSE.
         LOGICAL                                            :: isexplicit_Sigma_NyLeft = .FALSE.
         LOGICAL                                            :: isexplicit_Sigma_NxRight = .FALSE.
         LOGICAL                                            :: isexplicit_Sigma_NyRight = .FALSE.

         LOGICAL                                            :: isexplicit_GetRhoSingleLead = .FALSE.

         LOGICAL                                            :: isexplicit_deltamin = .FALSE.

         INTEGER                                            :: nprocs_inverse = 1
         INTEGER                                            :: Sigma_Nx = 1
         INTEGER                                            :: Sigma_Ny = 1

         INTEGER                                            :: nprocs_hs = 1
         INTEGER                                            :: Sigma_NxLeft = 1
         INTEGER                                            :: Sigma_NyLeft = 1
         INTEGER                                            :: Sigma_NxRight = 1
         INTEGER                                            :: Sigma_NyRight = 1

         INTEGER                                            :: maxdepth = 15
         INTEGER                                            :: ndivisions = 9
         INTEGER                                            :: storesigma = 0
         INTEGER                                            :: bs_method = 0
         INTEGER                                            :: bs_nmid = 0
         INTEGER                                            :: NRunSigmaMax = 3
         INTEGER                                            :: bsrun = 0
         INTEGER                                            :: bsskip = 1
         INTEGER                                            :: bseskip = 1
         INTEGER                                            :: bssc = 1
         INTEGER                                            :: inversion_solver = 1
         INTEGER                                            :: onekp = 0
         INTEGER                                            :: nbss = 2
         INTEGER                                            :: nSplit = 0
         INTEGER                                            :: TransmissionChannelsWFSSkipEne = 100
         INTEGER                                            :: TransmissionChannelsWFSSkipKP = 100
         INTEGER                                            :: TransmissionMatrixPDOSNWrite = 20
         INTEGER                                            :: TransmissionMatrixiSetPhase = -1
         INTEGER                                            :: MinChannelIndex = 1
         INTEGER                                            :: MaxChannelIndex = 0
         INTEGER                                            :: N_Replace_L = 1
         INTEGER                                            :: N_Replace_R = 1
         INTEGER                                            :: curr_nb1 = -1
         INTEGER                                            :: curr_nb2 = -1
         INTEGER                                            :: NParallelK = 1

         INTEGER                                            :: GetRhoSingleLead = 0

         INTEGER                                            :: sigma_method = 1

         REAL(kind=dp)                                      :: bs_tol = 1.0e-05_dp
         REAL(kind=dp)                                      :: bs_min = 0.0_dp
         REAL(kind=dp)                                      :: HSTol = 1.0e-10_dp
         REAL(kind=dp)                                      :: deltaimag = 0.0_dp
         REAL(kind=dp)                                      :: deltaini = 1.0e-03_dp
         REAL(kind=dp)                                      :: deltatode = 1.5_dp
         REAL(kind=dp)                                      :: critam = 2.0_dp
         REAL(kind=dp)                                      :: setemin = -1.0_dp
         REAL(kind=dp)                                      :: setemax = 1.0_dp
         REAL(kind=dp)                                      :: tolab = 1.0e-06_dp
         REAL(kind=dp)                                      :: weightrho = 0.5_dp
         REAL(kind=dp)                                      :: ThetaLeadsL = 0.0_dp
         REAL(kind=dp)                                      :: PhiLeadsL = 0.0_dp
         REAL(kind=dp)                                      :: ThetaLeadsR = 0.0_dp
         REAL(kind=dp)                                      :: PhiLeadsR = 0.0_dp
         REAL(kind=dp)                                      :: deltaimagtrc = 0.0_dp
         REAL(kind=dp)                                      :: m_tolki = 1.0e-06_dp
         REAL(kind=dp)                                      :: m_svdtolmax = 5.0e-12_dp
         REAL(kind=dp)                                      :: m_svdtolmin = 1.0e-15_dp
         REAL(kind=dp)                                      :: m_dsigmamax = 1.0e-05_dp
         REAL(kind=dp)                                      :: m_rnoise = 1.0_dp
         REAL(kind=dp)                                      :: m_skipsvd = 1.0e+03_dp
         REAL(kind=dp)                                      :: m_svdtolzi = 1.0e-08_dp
         REAL(kind=dp)                                      :: SigmaWideBand = 0.0_dp
         REAL(kind=dp)                                      :: emSTTkScale = 1.0_dp
         REAL(kind=dp)                                      :: nlSetZeroRatio = 1.0_dp
         REAL(kind=dp)                                      :: nrSetZeroRatio = 1.0_dp
         REAL(kind=dp)                                      :: curr_fl_L = 1.0_dp
         REAL(kind=dp)                                      :: curr_fr_L = 0.0_dp
         REAL(kind=dp)                                      :: curr_fl_R = 0.0_dp
         REAL(kind=dp)                                      :: curr_fr_R = 1.0_dp

         REAL(kind=dp)                                      :: deltamin = 0.0_dp

         REAL(kind=dp)                                      :: tolsigma = 1.0e-05_dp*0.5_dp
         REAL(kind=dp)                                      :: trcde = 0.0_dp*0.5_dp

      ! interface options
         REAL(kind=dp)                                      :: HartreeLeadsLeft = 0.0_dp
         REAL(kind=dp)                                      :: HartreeLeadsRight = 0.0_dp
         REAL(kind=dp)                                      :: HartreeLeadsBottom = 0.0_dp*0.5_dp
         LOGICAL                                            :: isexplicit_HartreeLeadsBottom = .FALSE.

   END TYPE smeagol_aux_control_type
CONTAINS

! **************************************************************************************************
!> \brief allocate control options for SMEAGOL calculation
!> \param smeagol_control an object to create
! **************************************************************************************************
   SUBROUTINE smeagol_control_create(smeagol_control)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      CHARACTER(len=*), PARAMETER :: routineN = 'smeagol_control_create'

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(smeagol_control))
      CALL timeset(routineN, handle)

      ALLOCATE (smeagol_control)
#if defined(__SMEAGOL)
      ALLOCATE (smeagol_control%aux)
#else
      NULLIFY (smeagol_control%aux)
#endif

      CALL timestop(handle)
   END SUBROUTINE smeagol_control_create

! **************************************************************************************************
!> \brief release SMEAGOL control object
!> \param smeagol_control an object to release
! **************************************************************************************************
   SUBROUTINE smeagol_control_release(smeagol_control)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control

      CHARACTER(len=*), PARAMETER :: routineN = 'smeagol_control_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(smeagol_control)) THEN
         IF (ASSOCIATED(smeagol_control%aux)) THEN
            IF (ALLOCATED(smeagol_control%aux%nebss_bs)) DEALLOCATE (smeagol_control%aux%nebss_bs)
            IF (ALLOCATED(smeagol_control%aux%deltabss_bs)) DEALLOCATE (smeagol_control%aux%deltabss_bs)
            IF (ALLOCATED(smeagol_control%aux%atomlist_bs)) DEALLOCATE (smeagol_control%aux%atomlist_bs)

            DEALLOCATE (smeagol_control%aux)
         END IF
         DEALLOCATE (smeagol_control)
      END IF

      CALL timestop(handle)
   END SUBROUTINE smeagol_control_release

! **************************************************************************************************
!> \brief Read SMEAGOL-related input parameters.
!> \param smeagol_control SMEAGOL control parameters
!> \param smeagol_section SMEAGOL input section
! **************************************************************************************************
   SUBROUTINE read_smeagol_control(smeagol_control, smeagol_section)
      TYPE(smeagol_control_type), POINTER                :: smeagol_control
      TYPE(section_vals_type), POINTER                   :: smeagol_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_smeagol_control'

      CHARACTER(len=default_string_length)               :: project_name
      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: n_cell_images_ptr

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(smeagol_section, "_SECTION_PARAMETERS_", l_val=smeagol_control%smeagol_enabled)

      ! SMEAGOL project name
      CALL section_vals_val_get(smeagol_section, "PROJECT_NAME", c_val=project_name)
      smeagol_control%project_name = project_name(:20)
      IF (smeagol_control%smeagol_enabled .AND. smeagol_control%project_name /= project_name) THEN
         CALL cp_warn(__LOCATION__, &
                      "SMEAGOL limits the length of the project name to 20 characters. "// &
                      "The project name is truncated to '"//TRIM(smeagol_control%project_name)//"'.")
      END IF

      ! reademtr() enum keywords
      CALL section_vals_val_get(smeagol_section, "RUN_TYPE", i_val=smeagol_control%run_type)

#if !defined(__SMEAGOL)
      IF (smeagol_control%run_type == smeagol_runtype_emtransport) THEN
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support. SMEAGOL RUN_TYPE EMTransport is not available")
      END IF
#endif

      CALL section_vals_val_get(smeagol_section, "REGRESSION_TEST", l_val=smeagol_control%do_regtest)

      CALL section_vals_val_get(smeagol_section, "BulkLead", i_val=smeagol_control%lead_label)

      ! NOTE: keyword NIMAGES_IJ seems to be identical to ignored SMEAGOL keywords BulkTransvCellSizeX and BulkTransvCellSizeY
      NULLIFY (n_cell_images_ptr)
      CALL section_vals_val_get(smeagol_section, "NIMAGES_IJ", i_vals=n_cell_images_ptr)
      smeagol_control%n_cell_images(1:2) = n_cell_images_ptr(1:2)

      ! Hartree > Rydberg scaling factor
      smeagol_control%to_smeagol_energy_units = cp_unit_from_cp2k(1.0_dp, "RY")

      IF (ASSOCIATED(smeagol_control%aux)) CALL read_smeagol_aux_control(smeagol_control%aux, smeagol_section)

      CALL timestop(handle)
   END SUBROUTINE read_smeagol_control

! **************************************************************************************************
!> \brief Read SMEAGOL-related auxiliary input parameters.
!> \param smeagol_control SMEAGOL control parameters
!> \param smeagol_section SMEAGOL input section
! **************************************************************************************************
   SUBROUTINE read_smeagol_aux_control(smeagol_control, smeagol_section)
      TYPE(smeagol_aux_control_type), POINTER            :: smeagol_control
      TYPE(section_vals_type), POINTER                   :: smeagol_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_smeagol_aux_control'

      CHARACTER(len=default_string_length)               :: nvals_str
      INTEGER                                            :: handle, i, n, nrep, n_unique
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ivec, indices
      INTEGER, DIMENSION(:), POINTER                     :: iptr
      LOGICAL                                            :: is_explicit, do_abort
      REAL(kind=dp), DIMENSION(:), POINTER               :: rptr
      TYPE(section_vals_type), POINTER                   :: subsection

      CALL timeset(routineN, handle)

      smeagol_control%md_iter_level = -1
      smeagol_control%md_first_step = 0

      CALL section_vals_val_get(smeagol_section, "TEMPERATURE", r_val=smeagol_control%temperature)

      ! reademtr() logical keywords
      ! The following keywords (variables) that are part of the reademtr_local_llist list are read but unused
      ! CB.WriteComplexBands (WriteComplexBands = .FALSE.) unimplemented.
      ! Presumably, some of these keywords should be removed.
         CALL section_vals_val_get(smeagol_section, "TrCoefficients", l_val=smeagol_control%TrCoeff)
         CALL section_vals_val_get(smeagol_section, "SkipLastIter", l_val=smeagol_control%skiplast)
         CALL section_vals_val_get(smeagol_section, "CalcIETS", l_val=smeagol_control%CalcIETS)
         CALL section_vals_val_get(smeagol_section, "Sigma.DSigmaDE", l_val=smeagol_control%dsigmade)
         CALL section_vals_val_get(smeagol_section, "CalcMPSH", l_val=smeagol_control%calcmpsh)
         CALL section_vals_val_get(smeagol_section, "EM.COOPCalculate", l_val=smeagol_control%COOPCalculate)

      ! reademtr() integer keywords
         CALL section_vals_val_get(smeagol_section, "NEnergReal", i_val=smeagol_control%NEnergR)
         CALL section_vals_val_get(smeagol_section, "NenergImCircle", i_val=smeagol_control%NEnergIC)
         CALL section_vals_val_get(smeagol_section, "NenergImLine", i_val=smeagol_control%NEnergIL)
         CALL section_vals_val_get(smeagol_section, "NPoles", i_val=smeagol_control%NPoles)
         CALL section_vals_val_get(smeagol_section, "SpinConfLeads", i_val=smeagol_control%SpinCL)
         CALL section_vals_val_get(smeagol_section, "AtomLeftVcte", i_val=smeagol_control%AtmLVCte)
         CALL section_vals_val_get(smeagol_section, "AtomRightVcte", i_val=smeagol_control%AtmRVCte)
         CALL section_vals_val_get(smeagol_section, "NSlices", i_val=smeagol_control%NSlices)
         CALL section_vals_val_get(smeagol_section, "NumberLinearMix", i_val=smeagol_control%nlinmix)
         CALL section_vals_val_get(smeagol_section, "MD.TRCSkip", i_val=smeagol_control%tmdskip)
         CALL section_vals_val_get(smeagol_section, "MD.TRCSampling", i_val=smeagol_control%tmdsampling)
         CALL section_vals_val_get(smeagol_section, "NTransmPoints", i_val=smeagol_control%NeneT)
         CALL section_vals_val_get(smeagol_section, "MPSHAtomFirst", i_val=smeagol_control%rmoinit)
         CALL section_vals_val_get(smeagol_section, "MPSHAtomLast", i_val=smeagol_control%rmolast)
         CALL section_vals_val_get(smeagol_section, "MPSHOrbFirst", i_val=smeagol_control%rorbinit)
         CALL section_vals_val_get(smeagol_section, "MPSHOrbLast", i_val=smeagol_control%rorblast)
         CALL section_vals_val_get(smeagol_section, "EM.COOPNumberOfBonds", i_val=smeagol_control%COOPNumberOfBonds)

      ! reademtr() real-valued keywords
         CALL section_vals_val_get(smeagol_section, "Delta", r_val=smeagol_control%Delta)
         CALL section_vals_val_get(smeagol_section, "CB.MaxKappa", r_val=smeagol_control%maxkappa)

      ! reademtr() physical-valued keywords
         CALL section_vals_val_get(smeagol_section, "VBias", r_val=smeagol_control%VBias)
         CALL section_vals_val_get(smeagol_section, "EnergLowestBound", r_val=smeagol_control%EnergLB)
         CALL section_vals_val_get(smeagol_section, "ZLeftVcte", r_val=smeagol_control%minL)
         CALL section_vals_val_get(smeagol_section, "ZRightVcte", &
                                   r_val=smeagol_control%maxR, &
                                   explicit=smeagol_control%isexplicit_maxR)
         CALL section_vals_val_get(smeagol_section, "InitTransmRange", r_val=smeagol_control%TEnergI)
         CALL section_vals_val_get(smeagol_section, "FinalTransmRange", r_val=smeagol_control%TEnergF)
         CALL section_vals_val_get(smeagol_section, "VGate", r_val=smeagol_control%vgate)
         CALL section_vals_val_get(smeagol_section, "ZVGateR", r_val=smeagol_control%zgater)
         CALL section_vals_val_get(smeagol_section, "ZVGateL", r_val=smeagol_control%zgatel)
         CALL section_vals_val_get(smeagol_section, "ZBroadeningG", r_val=smeagol_control%alphag)

      ! ReadOptionsNEGF_DFT() logical keywords
         CALL section_vals_val_get(smeagol_section, "EM.TimeReversal", l_val=smeagol_control%timereversal)
         CALL section_vals_val_get(smeagol_section, "EM.TRCAddVCDFT", l_val=smeagol_control%TRCAddVCDFT)
         CALL section_vals_val_get(smeagol_section, "EM.AddVgIsolatedLocalCharges", l_val=smeagol_control%AddVgCompensatingIsolated&
             &LocalCharges)
         CALL section_vals_val_get(smeagol_section, "EM.AddRhoGate", l_val=smeagol_control%Add_Rho_Gate)
         CALL section_vals_val_get(smeagol_section, "EM.DebugRhoGate", l_val=smeagol_control%debug_gate)
         CALL section_vals_val_get(smeagol_section, "EMPDOSKSO", l_val=smeagol_control%empdoskso)
         CALL section_vals_val_get(smeagol_section, "EM.Timings", l_val=smeagol_control%emtimings)

      ! ReadOptionsNEGF_DFT() integer keywords
         CALL section_vals_val_get(smeagol_section, "Ik_Select", i_val=smeagol_control%ik_select)

      ! ReadOptionsNEGF_DFT() real-valued keywords
         CALL section_vals_val_get(smeagol_section, "EM.NetRhoGateCharge", r_val=smeagol_control%Q_gate)

      ! ReadOptionsNEGF_DFT() physical-valued keywords
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLxMin", &
                                   r_val=smeagol_control%RGlxmin, &
                                   explicit=smeagol_control%isexplicit_RGlxmin)
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLxMax", &
                                   r_val=smeagol_control%RGlxmax, &
                                   explicit=smeagol_control%isexplicit_RGlxmax)
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLyMin", &
                                   r_val=smeagol_control%RGlymin, &
                                   explicit=smeagol_control%isexplicit_RGlymin)
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLyMax", &
                                   r_val=smeagol_control%RGlymax, &
                                   explicit=smeagol_control%isexplicit_RGlymax)
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLzMin", &
                                   r_val=smeagol_control%RGlzmin, &
                                   explicit=smeagol_control%isexplicit_RGlzmin)
         CALL section_vals_val_get(smeagol_section, "EM.RhoGateLzMax", &
                                   r_val=smeagol_control%RGlzmax, &
                                   explicit=smeagol_control%isexplicit_RGlzmax)

         CALL section_vals_val_get(smeagol_section, "ZeemanTermBx", r_val=smeagol_control%ZeemanBx)
         CALL section_vals_val_get(smeagol_section, "ZeemanTermBy", r_val=smeagol_control%ZeemanBy)
         CALL section_vals_val_get(smeagol_section, "ZeemanTermBz", r_val=smeagol_control%ZeemanBz)

      ! emtrans_options() enum keywords
      CALL section_vals_val_get(smeagol_section, "EnergyGridType", i_val=smeagol_control%gridmethod)
      CALL section_vals_val_get(smeagol_section, "TypeOfIntegral", i_val=smeagol_control%integraltype)

      ! emtrans_options() logical keywords
         CALL section_vals_val_get(smeagol_section, "Sigma.PDOS", l_val=smeagol_control%pdosgs)
         CALL section_vals_val_get(smeagol_section, "Sigma.SkipRight", l_val=smeagol_control%skipright)
         CALL section_vals_val_get(smeagol_section, "Sigma.WriteToDisk", l_val=smeagol_control%sigmatodisk)
         CALL section_vals_val_get(smeagol_section, "BS.Add", l_val=smeagol_control%bs_add)
         CALL section_vals_val_get(smeagol_section, "BS.WriteTransmission", l_val=smeagol_control%bs_writetrc)
         CALL section_vals_val_get(smeagol_section, "Current.WriteTransmission", l_val=smeagol_control%writetk)
         CALL section_vals_val_get(smeagol_section, "SetEnergyRange", l_val=smeagol_control%setene)
         CALL section_vals_val_get(smeagol_section, "Sigma.OverwriteHS", l_val=smeagol_control%overwritehs)
         CALL section_vals_val_get(smeagol_section, "TransmissionEmbed", l_val=smeagol_control%transmembed)
         CALL section_vals_val_get(smeagol_section, "Sigma.WriteEV", l_val=smeagol_control%writeevsigma)
         CALL section_vals_val_get(smeagol_section, "Sigma.SkipTransmission", l_val=smeagol_control%skiptransm)
         CALL section_vals_val_get(smeagol_section, "TRCDEAuto", l_val=smeagol_control%deauto)
         CALL section_vals_val_get(smeagol_section, "Sigma.CheckAccuracy", l_val=smeagol_control%geterrsigma)
         CALL section_vals_val_get(smeagol_section, "Sigma.Write", l_val=smeagol_control%getsigma)
         CALL section_vals_val_get(smeagol_section, "TRCScaleEf", l_val=smeagol_control%trcef)
         CALL section_vals_val_get(smeagol_section, "EM.SetRhoZeroIfHZero", l_val=smeagol_control%RhoSetZeroIfHZero)
         CALL section_vals_val_get(smeagol_section, "Sigma.InvertH1", l_val=smeagol_control%m_usehinv)
         CALL section_vals_val_get(smeagol_section, "Sigma.InvertV", l_val=smeagol_control%m_usevinv)
         CALL section_vals_val_get(smeagol_section, "Sigma.CSVD", l_val=smeagol_control%m_callsvd)
         CALL section_vals_val_get(smeagol_section, "Sigma.ComplexBands", l_val=smeagol_control%m_complexbands)
         CALL section_vals_val_get(smeagol_section, "Sigma.DOSVV", l_val=smeagol_control%m_dosleads)
         CALL section_vals_val_get(smeagol_section, "TRC.LeadsDOS", l_val=smeagol_control%leadsdos)
         CALL section_vals_val_get(smeagol_section, "TRC.LeadsPDOS", l_val=smeagol_control%leadspdos)
         CALL section_vals_val_get(smeagol_section, "EM.LDOS3", l_val=smeagol_control%emldos)
         CALL section_vals_val_get(smeagol_section, "EM.LDOS", l_val=smeagol_control%emldos2)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentDensityK_Ene", l_val=smeagol_control%curr_distKEne)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentDensityK", l_val=smeagol_control%curr_distK)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentDensity", l_val=smeagol_control%curr_dist)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentDensityTRC", l_val=smeagol_control%curr_distTRC)
         CALL section_vals_val_get(smeagol_section, "EM.NonSelfConsistentRun", l_val=smeagol_control%em_NonSCF)
         CALL section_vals_val_get(smeagol_section, "EM.TransmissionRL", l_val=smeagol_control%TransmissionRL)
         CALL section_vals_val_get(smeagol_section, "EM.SpinTorqueK", l_val=smeagol_control%emSTTk)
         CALL section_vals_val_get(smeagol_section, "EM.TRCSpinComponents", l_val=smeagol_control%GetT0S)
         CALL section_vals_val_get(smeagol_section, "EM.TRCChannels", l_val=smeagol_control%TransmissionChannels)
         CALL section_vals_val_get(smeagol_section, "EM.TRCChannelsWFS", l_val=smeagol_control%TransmissionChannelsWFS)
         CALL section_vals_val_get(smeagol_section, "EM.TMatrix", l_val=smeagol_control%TransmissionMatrix)
         CALL section_vals_val_get(smeagol_section, "EM.TMatrixWFS", l_val=smeagol_control%TransmissionMatrixWFS)
         CALL section_vals_val_get(smeagol_section, "EM.TMatrixPDOS", l_val=smeagol_control%TransmissionMatrixPDOS)
         CALL section_vals_val_get(smeagol_section, "EM.EndCode1", l_val=smeagol_control%EM_EndCode1)
         CALL section_vals_val_get(smeagol_section, "TransmissionOverK", l_val=smeagol_control%TransmOverk)
         CALL section_vals_val_get(smeagol_section, "TRC.EMDOS", l_val=smeagol_control%emdos)
         CALL section_vals_val_get(smeagol_section, "TRC.EMPDOS", l_val=smeagol_control%empdos)
         CALL section_vals_val_get(smeagol_section, "TRC.EMPDOSK", l_val=smeagol_control%empdosk)
         CALL section_vals_val_get(smeagol_section, "TRC.TrcOverKWriteIk", l_val=smeagol_control%WriteIkTrcK)
         CALL section_vals_val_get(smeagol_section, "CBS.SubSystem", l_val=smeagol_control%cbspart)
         CALL section_vals_val_get(smeagol_section, "CBS.PrintWaveFunction", l_val=smeagol_control%evec)
         CALL section_vals_val_get(smeagol_section, "EM.OrderN", l_val=smeagol_control%negfon)
         CALL section_vals_val_get(smeagol_section, "EM.OutputInfo", l_val=smeagol_control%outinfo)
         CALL section_vals_val_get(smeagol_section, "EM.WriteNk", l_val=smeagol_control%WriteNk)
         CALL section_vals_val_get(smeagol_section, "EM.SetRhoBoundaryFromLeads", l_val=smeagol_control%Set_RhoBoundary_Leads)
         CALL section_vals_val_get(smeagol_section, "EM.SetRhoBoundaryOverlapFromLeads",&
             & l_val=smeagol_control%Set_RhoBoundaryOverlap_Leads)
         CALL section_vals_val_get(smeagol_section, "EM.SetHBoundaryFromLeads", l_val=smeagol_control%Set_HBoundary_Leads)
         CALL section_vals_val_get(smeagol_section, "EM.SetHLRZero", l_val=smeagol_control%Set_HLR_Zero)
         CALL section_vals_val_get(smeagol_section, "EM.SkipNonEquilibriumRho", l_val=smeagol_control%SkipNonEquilibriumRho)
         CALL section_vals_val_get(smeagol_section, "EM.SkipEquilibriumRho", l_val=smeagol_control%SkipEquilibriumRho)
         CALL section_vals_val_get(smeagol_section, "EM.ComputeImpurityGfMatsubara", l_val=smeagol_control%ComputeImpurityGfMatsuba&
             &ra)
         CALL section_vals_val_get(smeagol_section, "EM.PrintImpurityGfMatsubara", l_val=smeagol_control%PrintImpurityGfMatsubara)
         CALL section_vals_val_get(smeagol_section, "EM.CallImpuritySolver", l_val=smeagol_control%CallImpuritySolver)

      ! emtrans_options() integer keywords
         CALL section_vals_val_get(smeagol_section, "EM.NProcessorsInverse", i_val=smeagol_control%nprocs_inverse)
         CALL section_vals_val_get(smeagol_section, "Sigma.Nx", i_val=smeagol_control%Sigma_Nx)
         CALL section_vals_val_get(smeagol_section, "Sigma.Ny", i_val=smeagol_control%Sigma_Ny)

         CALL section_vals_val_get(smeagol_section, "EM.NProcessorsDistributeMatrices", &
                                   i_val=smeagol_control%nprocs_hs, &
                                   explicit=smeagol_control%isexplicit_nprocs_hs)
         CALL section_vals_val_get(smeagol_section, "Sigma.NxLeft", &
                                   i_val=smeagol_control%Sigma_NxLeft, &
                                   explicit=smeagol_control%isexplicit_Sigma_NxLeft)
         CALL section_vals_val_get(smeagol_section, "Sigma.NyLeft", &
                                   i_val=smeagol_control%Sigma_NyLeft, &
                                   explicit=smeagol_control%isexplicit_Sigma_NyLeft)
         CALL section_vals_val_get(smeagol_section, "Sigma.NxRight", &
                                   i_val=smeagol_control%Sigma_NxRight, &
                                   explicit=smeagol_control%isexplicit_Sigma_NxRight)
         CALL section_vals_val_get(smeagol_section, "Sigma.NyRight", &
                                   i_val=smeagol_control%Sigma_NyRight, &
                                   explicit=smeagol_control%isexplicit_Sigma_NyRight)

         CALL section_vals_val_get(smeagol_section, "AM.MaxDepth", i_val=smeagol_control%maxdepth)
         CALL section_vals_val_get(smeagol_section, "AM.NumberDivisions", i_val=smeagol_control%ndivisions)
         CALL section_vals_val_get(smeagol_section, "SIGMA.Save", i_val=smeagol_control%storesigma)
         CALL section_vals_val_get(smeagol_section, "BS.Method", i_val=smeagol_control%bs_method)
         CALL section_vals_val_get(smeagol_section, "BS.MiddleOrbital", i_val=smeagol_control%bs_nmid)
         CALL section_vals_val_get(smeagol_section, "Sigma.NRerunSigmaMax", i_val=smeagol_control%NRunSigmaMax)
         CALL section_vals_val_get(smeagol_section, "BS.TypeOfRun", i_val=smeagol_control%bsrun)
         CALL section_vals_val_get(smeagol_section, "BS.Skip", i_val=smeagol_control%bsskip)
         CALL section_vals_val_get(smeagol_section, "BS.ESkip", i_val=smeagol_control%bseskip)
         CALL section_vals_val_get(smeagol_section, "BS.SetOccupation", i_val=smeagol_control%bssc)
         CALL section_vals_val_get(smeagol_section, "EM.InverseSolver", i_val=smeagol_control%inversion_solver)
         CALL section_vals_val_get(smeagol_section, "EM.OneKP", i_val=smeagol_control%onekp)
         CALL section_vals_val_get(smeagol_section, "BS.Subsystems", i_val=smeagol_control%nbss)
         CALL section_vals_val_get(smeagol_section, "EM.TRCSpinComponentsOribalSplit", i_val=smeagol_control%nSplit)
         CALL section_vals_val_get(smeagol_section, "EM.TRCChannelsWFS_SkipEne", i_val=smeagol_control%TransmissionChannelsWFSSkipE&
             &ne)
         CALL section_vals_val_get(smeagol_section, "EM.TRCChannelsWFS_SkipKP", i_val=smeagol_control%TransmissionChannelsWFSSkipKP)
         CALL section_vals_val_get(smeagol_section, "EM.TMatrixPDOSNWrite", i_val=smeagol_control%TransmissionMatrixPDOSNWrite)
         CALL section_vals_val_get(smeagol_section, "EM.TMatrixSetPhaseIndex", i_val=smeagol_control%TransmissionMatrixiSetPhase)
         CALL section_vals_val_get(smeagol_section, "EM.TRCMinChannelIndex", i_val=smeagol_control%MinChannelIndex)
         CALL section_vals_val_get(smeagol_section, "EM.TRCMaxChannelIndex", i_val=smeagol_control%MaxChannelIndex)
         CALL section_vals_val_get(smeagol_section, "EM.NSetHBoundaryFromLeadsL", i_val=smeagol_control%N_Replace_L)
         CALL section_vals_val_get(smeagol_section, "EM.NSetHBoundaryFromLeadsR", i_val=smeagol_control%N_Replace_R)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxOrbital1", i_val=smeagol_control%curr_nb1)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxOrbital2", i_val=smeagol_control%curr_nb2)
         CALL section_vals_val_get(smeagol_section, "EM.ParallelOverKNum", i_val=smeagol_control%NParallelK)

         CALL section_vals_val_get(smeagol_section, "EM.LDOSLeadsProjection", &
                                   i_val=smeagol_control%GetRhoSingleLead, &
                                   explicit=smeagol_control%isexplicit_GetRhoSingleLead)

         CALL section_vals_val_get(smeagol_section, "Sigma.Method", i_val=smeagol_control%sigma_method)

      ! emtrans_options() real-valued keywords
         CALL section_vals_val_get(smeagol_section, "BS.Tolerance", r_val=smeagol_control%bs_tol)
         CALL section_vals_val_get(smeagol_section, "BS.Minimum", r_val=smeagol_control%bs_min)
         CALL section_vals_val_get(smeagol_section, "EM.HSTol", r_val=smeagol_control%HSTol)
         CALL section_vals_val_get(smeagol_section, "Sigma.EImag", r_val=smeagol_control%deltaimag)
         CALL section_vals_val_get(smeagol_section, "AM.DeltaInit", r_val=smeagol_control%deltaini)
         CALL section_vals_val_get(smeagol_section, "AM.DeltaToDe", r_val=smeagol_control%deltatode)
         CALL section_vals_val_get(smeagol_section, "AM.Criterium", r_val=smeagol_control%critam)
         CALL section_vals_val_get(smeagol_section, "SetEmin", r_val=smeagol_control%setemin)
         CALL section_vals_val_get(smeagol_section, "SetEmax", r_val=smeagol_control%setemax)
         CALL section_vals_val_get(smeagol_section, "Sigma.Tolab", r_val=smeagol_control%tolab)
         CALL section_vals_val_get(smeagol_section, "EM.WeightRho", r_val=smeagol_control%weightrho)
         CALL section_vals_val_get(smeagol_section, "EM.RotateSpinLeadsLeft.Theta", r_val=smeagol_control%ThetaLeadsL)
         CALL section_vals_val_get(smeagol_section, "EM.RotateSpinLeadsLeft.Phi", r_val=smeagol_control%PhiLeadsL)
         CALL section_vals_val_get(smeagol_section, "EM.RotateSpinLeadsRight.Theta", r_val=smeagol_control%ThetaLeadsR)
         CALL section_vals_val_get(smeagol_section, "EM.RotateSpinLeadsRight.Phi", r_val=smeagol_control%PhiLeadsR)
         CALL section_vals_val_get(smeagol_section, "DeltaTransmission", r_val=smeagol_control%deltaimagtrc)
         CALL section_vals_val_get(smeagol_section, "Sigma.Dkimag", r_val=smeagol_control%m_tolki)
         CALL section_vals_val_get(smeagol_section, "Sigma.DSVDMax", r_val=smeagol_control%m_svdtolmax)
         CALL section_vals_val_get(smeagol_section, "Sigma.DSVDMin", r_val=smeagol_control%m_svdtolmin)
         CALL section_vals_val_get(smeagol_section, "Sigma.DMax", r_val=smeagol_control%m_dsigmamax)
         CALL section_vals_val_get(smeagol_section, "Sigma.RNoise", r_val=smeagol_control%m_rnoise)
         CALL section_vals_val_get(smeagol_section, "Sigma.SkipSVD", r_val=smeagol_control%m_skipsvd)
         CALL section_vals_val_get(smeagol_section, "Sigma.SVDTolZero", r_val=smeagol_control%m_svdtolzi)
         CALL section_vals_val_get(smeagol_section, "Sigma.WideBandImag", r_val=smeagol_control%SigmaWideBand)
         CALL section_vals_val_get(smeagol_section, "EM.SpinTorqueKScale", r_val=smeagol_control%emSTTkScale)
         CALL section_vals_val_get(smeagol_section, "EM.SetHLRZeroRatioLeft", r_val=smeagol_control%nlSetZeroRatio)
         CALL section_vals_val_get(smeagol_section, "EM.SetHLRZeroRatioRight", r_val=smeagol_control%nrSetZeroRatio)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxFL_L", r_val=smeagol_control%curr_fl_L)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxFR_L", r_val=smeagol_control%curr_fr_L)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxFL_R", r_val=smeagol_control%curr_fl_R)
         CALL section_vals_val_get(smeagol_section, "EM.CurrentFluxFR_R", r_val=smeagol_control%curr_fr_R)

         CALL section_vals_val_get(smeagol_section, "BS.Delta", &
                                   r_val=smeagol_control%deltamin, &
                                   explicit=smeagol_control%isexplicit_deltamin)

      ! emtrans_options() physical-valued keywords
         CALL section_vals_val_get(smeagol_section, "Sigma.WarnInOutput", r_val=smeagol_control%tolsigma)
         CALL section_vals_val_get(smeagol_section, "TRCDE", r_val=smeagol_control%trcde)

      ! interface physical-valued keywords
         CALL section_vals_val_get(smeagol_section, "HartreeLeadsLeft", r_val=smeagol_control%HartreeLeadsLeft)
         CALL section_vals_val_get(smeagol_section, "HartreeLeadsRight", r_val=smeagol_control%HartreeLeadsRight)
         CALL section_vals_val_get(smeagol_section, "HartreeLeadsBottom", &
                                   r_val=smeagol_control%HartreeLeadsBottom, &
                                   explicit=smeagol_control%isexplicit_HartreeLeadsBottom)

      ! Bound states
      IF (smeagol_control%nbss <= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "BS.Subsystems keyword should have a positive integer value.")
      END IF

      subsection => section_vals_get_subs_vals(smeagol_section, "BS.SubSystemsBoundaries")
      CALL section_vals_get(subsection, explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(subsection, "_DEFAULT_KEYWORD_", n_rep_val=nrep, explicit=is_explicit)
         IF (.NOT. is_explicit) nrep = 0

         do_abort = .FALSE.
         IF (smeagol_control%nbss == nrep) THEN
            IF (is_explicit) THEN
               ALLOCATE (smeagol_control%nebss_bs(nrep, 2))
               DO i = 1, nrep
                  CALL section_vals_val_get(subsection, "_DEFAULT_KEYWORD_", i_rep_val=i, i_vals=iptr)
                  IF (SIZE(iptr) == 2) THEN
                     IF (MINVAL(iptr) > 0 .AND. iptr(1) <= iptr(2)) THEN
                        smeagol_control%nebss_bs(i, 1:2) = iptr(1:2)
                     ELSE
                        do_abort = .TRUE.
                        EXIT
                     END IF
                  ELSE
                     do_abort = .TRUE.
                     EXIT
                  END IF
               END DO
            END IF
         ELSE
            do_abort = .TRUE.
         END IF

         IF (do_abort) THEN
            CALL integer_to_string(smeagol_control%nbss, nvals_str)
            CALL cp_abort(__LOCATION__, &
                          "BS.SubSystemsBoundaries section is expected to have BS.Subsystems ("//TRIM(nvals_str)// &
                          ") lines with two positive integer numbers on each line.")
         END IF

      END IF

      CALL section_vals_val_get(smeagol_section, "BS.SubSystemsDelta", explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(smeagol_section, "BS.SubSystemsDelta", r_vals=rptr)
         IF (smeagol_control%nbss /= SIZE(rptr)) THEN ! do negative deltas make sense ? .OR. MINVAL(rptr) < 0
            CALL integer_to_string(smeagol_control%nbss, nvals_str)
            CALL cp_abort(__LOCATION__, &
                          "BS.SubSystemsDelta keyword is expected to have BS.Subsystems ("//TRIM(nvals_str)//") real numbers.")
         END IF

         n = SIZE(rptr)
         ALLOCATE (smeagol_control%deltabss_bs(n))
         smeagol_control%deltabss_bs(1:n) = rptr(1:n)
         !ELSE allocate and set deltabss_bs() to deltamin
      END IF

      CALL section_vals_val_get(smeagol_section, "AM.AtomListBS", explicit=is_explicit)
      IF (is_explicit) THEN
         CALL section_vals_val_get(smeagol_section, "AM.AtomListBS", i_vals=iptr)

         ! subsys is not available at this stage, so we cannot check that MINVAL(iptr) <= natoms so far
         IF (SIZE(iptr) == 0 .OR. MINVAL(iptr) <= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "All atomic indices in AM.AtomListBS should be positive integer numbers.")
         END IF

         n = SIZE(iptr)
         ALLOCATE (ivec(n), indices(n))
         ivec(1:n) = iptr(1:n)
         CALL sort(ivec, n, indices)

         n_unique = 1
         DO i = 2, n
            IF (ivec(i) > ivec(i - 1)) n_unique = n_unique + 1
         END DO

         ALLOCATE (smeagol_control%atomlist_bs(n_unique))
         n_unique = 1
         smeagol_control%atomlist_bs(1) = ivec(1)

         DO i = 2, n
            IF (ivec(i) > ivec(i - 1)) THEN
               n_unique = n_unique + 1
               smeagol_control%atomlist_bs(n_unique) = ivec(i)
            END IF
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE read_smeagol_aux_control
END MODULE smeagol_control_types
