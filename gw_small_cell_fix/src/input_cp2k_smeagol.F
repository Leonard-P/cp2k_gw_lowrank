!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Functions that build SMEAGOL input section.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE input_cp2k_smeagol
   USE bibliography, ONLY: Ahart2024
   USE cp_output_handling, ONLY: add_last_numeric, &
                                 cp_print_key_section_create, &
                                 debug_print_level, &
                                 high_print_level, &
                                 low_print_level, &
                                 medium_print_level, &
                                 silent_print_level
   USE cp_units, ONLY: cp_unit_to_cp2k
   USE input_constants, ONLY: smeagol_bulklead_left, &
                              smeagol_bulklead_leftright, &
                              smeagol_bulklead_right, &
                              smeagol_gridmethod_adaptive, &
                              smeagol_gridmethod_traditional, &
                              smeagol_integraltype_gauss_chebyshev, &
                              smeagol_integraltype_gauss_legendre, &
                              smeagol_runtype_bulktransport, &
                              smeagol_runtype_emtransport
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type
   USE input_val_types, ONLY: char_t, &
                              integer_t, &
                              logical_t, &
                              real_t
   USE kinds, ONLY: dp
   USE physcon, ONLY: kelvin
   USE string_utilities, ONLY: s2a
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_smeagol'

   PUBLIC :: create_dft_smeagol_section

CONTAINS

! **************************************************************************************************
!> \brief Create SMEAGOL input section.
!> \param section input section
! **************************************************************************************************
   SUBROUTINE create_dft_smeagol_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      ! +11 -- explicitly defined keywords
      !       _SECTION_PARAMETERS_, BS.SubSystemsDelta, AM.AtomListBS, TEMPERATURE, PROJECT_NAME, RUN_TYPE, REGRESSION_TEST,
      !        BulkLead, NIMAGES_IJ, EnergyGridType, TypeOfIntegral
      CALL section_create(section, __LOCATION__, name="SMEAGOL", &
                          description="Parameters for NEGF/SMEAGOL simulation.", &
                          n_subsections=1, repeats=.FALSE., &
                          citations=[Ahart2024])

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls activation of SMEAGOL interface", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! +++ reademtr()

      ! BS.SubSystemsBoundaries subsection
      NULLIFY (subsection)
      CALL section_create(subsection, __LOCATION__, name="BS.SubSystemsBoundaries", &
                          description="Bound state subsystems' boundaries", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL keyword_create(keyword, __LOCATION__, name="_DEFAULT_KEYWORD_", &
                          description="BS.Subsystems lines containing two integer numbers", &
                          usage="lower_bound upper_bound", &
                          repeats=.TRUE., n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      ! BS.SubSystemsDelta
      CALL keyword_create(keyword, __LOCATION__, name="BS.SubSystemsDelta", &
                          description="Specifies a list of atoms to include.", &
                          usage="BS.SubSystemsDelta {real} {real} .. {real}", repeats=.FALSE., &
                          n_var=-1, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! use a list of atomic indices (AM.SpeciesBS) instead of list of atomic kind labels (AM.SpeciesBS)
      CALL keyword_create(keyword, __LOCATION__, name="AM.AtomListBS", &
                          description="Specifies a list of atoms to include.", &
                          usage="AM.AtomListBS {integer} {integer} .. {integer}", repeats=.FALSE., &
                          n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Electronic temperature of the extended molecule.
      CALL keyword_create(keyword, __LOCATION__, name="TEMPERATURE", &
                          description="Electronic temperature of the extended molecule."// &
                          " Temperature of semi-infinite leads are taken from SMEAR / ELECTRONIC_TEMPERATURE"// &
                          " keyword of relevant bulk-transport calculations.", &
                          n_var=1, type_of_var=real_t, unit_str="K", &
                          default_r_val=300.0_dp/kelvin)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROJECT_NAME", &
                          description="Name of this SMEAGOL project. SMEAGOL writes certain quantities into files named as "// &
                          "smeagol_project.EXTENSION. The length of the project name should not exceed 20 characters.", &
                          n_var=1, type_of_var=char_t, default_c_val="smeagol_project")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! reademtr() enum keywords

      CALL keyword_create(keyword, __LOCATION__, name="RUN_TYPE", &
                          description="Type of SMEAGOL calculation", &
                          default_i_val=smeagol_runtype_bulktransport, &
                          enum_c_vals=s2a("BulkTransport", "EMTransport"), &
                          enum_desc=s2a( &
                          "Save Kohn-Sham and density matrices of semi-infinite leads.", &
                          "Run transport calculation for the extended molecule."), &
                          enum_i_vals=[smeagol_runtype_bulktransport, smeagol_runtype_emtransport])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="REGRESSION_TEST", &
                          description="Enables regression test mode. "// &
                          "Kohn-Sham and density matrices of semi-infinite leads are not written down.", &
                          n_var=1, type_of_var=logical_t, default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="BulkLead", &
                          description="Lead label for bulk transport calculation", &
                          default_i_val=smeagol_bulklead_leftright, &
                          enum_c_vals=s2a("LR", "L", "R"), &
                          enum_desc=s2a( &
                          "Identical left and right leads.", &
                          "Left bulk lead.", &
                          "Right bulk lead."), &
                          enum_i_vals=[smeagol_bulklead_leftright, smeagol_bulklead_left, smeagol_bulklead_right])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! force the number of cell images along i and j cell vectors
      CALL keyword_create(keyword, __LOCATION__, name="NIMAGES_IJ", &
                          description="Number of cell images along i and j cell vectors. "// &
                          "If missing, determine the number of images based on sparsity pattern of the overlap matrix.", &
                          usage="NIMAGES_IJ 3 3", &
                          n_var=2, type_of_var=integer_t, default_i_vals=[0, 0])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! reademtr() logical keywords
         CALL keyword_create(keyword, __LOCATION__, name="TrCoefficients", &
                             description="TrCoefficients", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SkipLastIter", &
                             description="SkipLastIter", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="CalcIETS", &
                             description="CalcIETS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="Sigma.DSigmaDE", &
                             description="Sigma.DSigmaDE", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="CalcMPSH", &
                             description="CalcMPSH", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="EM.COOPCalculate", &
                             description="EM.COOPCalculate", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! reademtr() integer keywords
         CALL keyword_create(keyword, __LOCATION__, name="NEnergReal", &
                             description="NEnergReal", &
                             n_var=1, type_of_var=integer_t, default_i_val=288)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="NenergImCircle", &
                             description="NenergImCircle", &
                             n_var=1, type_of_var=integer_t, default_i_val=32)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="NenergImLine", &
                             description="NenergImLine", &
                             n_var=1, type_of_var=integer_t, default_i_val=16)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="NPoles", &
                             description="NPoles", &
                             n_var=1, type_of_var=integer_t, default_i_val=16)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SpinConfLeads", &
                             description="SpinConfLeads", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AtomLeftVcte", &
                             description="AtomLeftVcte", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AtomRightVcte", &
                             description="AtomRightVcte", &
                             n_var=1, type_of_var=integer_t, default_i_val=-1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="NSlices", &
                             description="NSlices", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="NumberLinearMix", &
                             description="NumberLinearMix", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MD.TRCSkip", &
                             description="MD.TRCSkip", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MD.TRCSampling", &
                             description="MD.TRCSampling", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="NTransmPoints", &
                             description="NTransmPoints", &
                             n_var=1, type_of_var=integer_t, default_i_val=100)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MPSHAtomFirst", &
                             description="MPSHAtomFirst", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MPSHAtomLast", &
                             description="MPSHAtomLast", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MPSHOrbFirst", &
                             description="MPSHOrbFirst", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="MPSHOrbLast", &
                             description="MPSHOrbLast", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="EM.COOPNumberOfBonds", &
                             description="EM.COOPNumberOfBonds", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! reademtr() real-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="Delta", &
                             description="Delta", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-10_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="CB.MaxKappa", &
                             description="CB.MaxKappa", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! reademtr() physical-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="VBias", &
                             description="VBias", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EnergLowestBound", &
                             description="EnergLowestBound", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=-7.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="ZLeftVcte", &
                             description="ZLeftVcte", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="ZRightVcte", &
                             description="ZRightVcte", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="InitTransmRange", &
                             description="InitTransmRange", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=-0.3675_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="FinalTransmRange", &
                             description="FinalTransmRange", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=-0.3675_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="VGate", &
                             description="VGate", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="ZVGateR", &
                             description="ZVGateR", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="ZVGateL", &
                             description="ZVGateL", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="ZBroadeningG", &
                             description="ZBroadeningG", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! +++ ReadOptionsNEGF_DFT()

      ! ReadOptionsNEGF_DFT() logical keywords
         CALL keyword_create(keyword, __LOCATION__, name="EM.TimeReversal", &
                             description="EM.TimeReversal", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCAddVCDFT", &
                             description="EM.TRCAddVCDFT", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.AddVgIsolatedLocalCharges", &
                             description="EM.AddVgIsolatedLocalCharges", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.AddRhoGate", &
                             description="EM.AddRhoGate", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.DebugRhoGate", &
                             description="EM.DebugRhoGate", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EMPDOSKSO", &
                             description="EMPDOSKSO", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.Timings", &
                             description="EM.Timings", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! ReadOptionsNEGF_DFT() integer keywords
         CALL keyword_create(keyword, __LOCATION__, name="Ik_Select", &
                             description="Ik_Select", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! ReadOptionsNEGF_DFT() real-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="EM.NetRhoGateCharge", &
                             description="EM.NetRhoGateCharge", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! ReadOptionsNEGF_DFT() physical-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLxMin", &
                             description="EM.RhoGateLxMin", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLxMax", &
                             description="EM.RhoGateLxMax", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLyMin", &
                             description="EM.RhoGateLyMin", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLyMax", &
                             description="EM.RhoGateLyMax", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLzMin", &
                             description="EM.RhoGateLzMin", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RhoGateLzMax", &
                             description="EM.RhoGateLzMax", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="ZeemanTermBx", &
                             description="ZeemanTermBx", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="ZeemanTermBy", &
                             description="ZeemanTermBy", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="ZeemanTermBz", &
                             description="ZeemanTermBz", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! +++ emtrans_options()

      ! emtrans_options() enum keywords
      CALL keyword_create(keyword, __LOCATION__, name="EnergyGridType", &
                          description="Type of energy grid", &
                          default_i_val=smeagol_gridmethod_traditional, &
                          enum_c_vals=s2a("TRADITIONAL", "ADAPTIVEGRID"), &
                          enum_desc=s2a( &
                          "Traditional grid.", &
                          "Adaptive grid."), &
                          enum_i_vals=[smeagol_gridmethod_traditional, smeagol_gridmethod_adaptive])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TypeOfIntegral", &
                          description="Method to integrate Green's functions", &
                          default_i_val=smeagol_integraltype_gauss_legendre, &
                          enum_c_vals=s2a("GAUSS-LEGENDRE", "GAUSS-CHEBYSHEV"), &
                          enum_desc=s2a( &
                          "Gauss-Legendre quadrature method.", &
                          "Gauss-Chebyshev quadrature method."), &
                          enum_i_vals=[smeagol_integraltype_gauss_legendre, smeagol_integraltype_gauss_chebyshev])
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! emtrans_options() logical keywords
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.PDOS", &
                             description="Sigma.PDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.SkipRight", &
                             description="Sigma.SkipRight", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.WriteToDisk", &
                             description="Sigma.WriteToDisk", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.Add", &
                             description="BS.Add", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.WriteTransmission", &
                             description="BS.WriteTransmission", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Current.WriteTransmission", &
                             description="Current.WriteTransmission", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SetEnergyRange", &
                             description="SetEnergyRange", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.OverwriteHS", &
                             description="Sigma.OverwriteHS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TransmissionEmbed", &
                             description="TransmissionEmbed", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.WriteEV", &
                             description="Sigma.WriteEV", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.SkipTransmission", &
                             description="Sigma.SkipTransmission", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRCDEAuto", &
                             description="TRCDEAuto", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.CheckAccuracy", &
                             description="Sigma.CheckAccuracy", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Write", &
                             description="Sigma.Write", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRCScaleEf", &
                             description="TRCScaleEf", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetRhoZeroIfHZero", &
                             description="EM.SetRhoZeroIfHZero", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.InvertH1", &
                             description="Sigma.InvertH1", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.InvertV", &
                             description="Sigma.InvertV", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.CSVD", &
                             description="Sigma.CSVD", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.ComplexBands", &
                             description="Sigma.ComplexBands", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.DOSVV", &
                             description="Sigma.DOSVV", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.LeadsDOS", &
                             description="TRC.LeadsDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.LeadsPDOS", &
                             description="TRC.LeadsPDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.LDOS3", &
                             description="EM.LDOS3", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.LDOS", &
                             description="EM.LDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentDensityK_Ene", &
                             description="EM.CurrentDensityK_Ene", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentDensityK", &
                             description="EM.CurrentDensityK", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentDensity", &
                             description="EM.CurrentDensity", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentDensityTRC", &
                             description="EM.CurrentDensityTRC", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.NonSelfConsistentRun", &
                             description="EM.NonSelfConsistentRun", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TransmissionRL", &
                             description="EM.TransmissionRL", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SpinTorqueK", &
                             description="EM.SpinTorqueK", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCSpinComponents", &
                             description="EM.TRCSpinComponents", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCChannels", &
                             description="EM.TRCChannels", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCChannelsWFS", &
                             description="EM.TRCChannelsWFS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TMatrix", &
                             description="EM.TMatrix", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TMatrixWFS", &
                             description="EM.TMatrixWFS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TMatrixPDOS", &
                             description="EM.TMatrixPDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.EndCode1", &
                             description="EM.EndCode1", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TransmissionOverK", &
                             description="TransmissionOverK", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.EMDOS", &
                             description="TRC.EMDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.EMPDOS", &
                             description="TRC.EMPDOS", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.EMPDOSK", &
                             description="TRC.EMPDOSK", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRC.TrcOverKWriteIk", &
                             description="TRC.TrcOverKWriteIk", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="CBS.SubSystem", &
                             description="CBS.SubSystem", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="CBS.PrintWaveFunction", &
                             description="CBS.PrintWaveFunction", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.OrderN", &
                             description="EM.OrderN", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.OutputInfo", &
                             description="EM.OutputInfo", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.WriteNk", &
                             description="EM.WriteNk", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetRhoBoundaryFromLeads", &
                             description="EM.SetRhoBoundaryFromLeads", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetRhoBoundaryOverlapFromLeads", &
                             description="EM.SetRhoBoundaryOverlapFromLeads", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetHBoundaryFromLeads", &
                             description="EM.SetHBoundaryFromLeads", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetHLRZero", &
                             description="EM.SetHLRZero", &
                             n_var=1, type_of_var=logical_t, default_l_val=.true., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SkipNonEquilibriumRho", &
                             description="EM.SkipNonEquilibriumRho", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SkipEquilibriumRho", &
                             description="EM.SkipEquilibriumRho", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.ComputeImpurityGfMatsubara", &
                             description="EM.ComputeImpurityGfMatsubara", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.PrintImpurityGfMatsubara", &
                             description="EM.PrintImpurityGfMatsubara", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CallImpuritySolver", &
                             description="EM.CallImpuritySolver", &
                             n_var=1, type_of_var=logical_t, default_l_val=.false., &
                             lone_keyword_l_val=.TRUE.)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="EM.NProcessorsInverse", &
                             description="EM.NProcessorsInverse", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Nx", &
                             description="Sigma.Nx", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Ny", &
                             description="Sigma.Ny", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="EM.NProcessorsDistributeMatrices", &
                             description="EM.NProcessorsDistributeMatrices", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.NxLeft", &
                             description="Sigma.NxLeft", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.NyLeft", &
                             description="Sigma.NyLeft", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.NxRight", &
                             description="Sigma.NxRight", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.NyRight", &
                             description="Sigma.NyRight", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="AM.MaxDepth", &
                             description="AM.MaxDepth", &
                             n_var=1, type_of_var=integer_t, default_i_val=15)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AM.NumberDivisions", &
                             description="AM.NumberDivisions", &
                             n_var=1, type_of_var=integer_t, default_i_val=9)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SIGMA.Save", &
                             description="SIGMA.Save", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.Method", &
                             description="BS.Method", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.MiddleOrbital", &
                             description="BS.MiddleOrbital", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.NRerunSigmaMax", &
                             description="Sigma.NRerunSigmaMax", &
                             n_var=1, type_of_var=integer_t, default_i_val=3)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.TypeOfRun", &
                             description="BS.TypeOfRun", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.Skip", &
                             description="BS.Skip", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.ESkip", &
                             description="BS.ESkip", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.SetOccupation", &
                             description="BS.SetOccupation", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.InverseSolver", &
                             description="EM.InverseSolver", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.OneKP", &
                             description="EM.OneKP", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.Subsystems", &
                             description="BS.Subsystems", &
                             n_var=1, type_of_var=integer_t, default_i_val=2)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCSpinComponentsOribalSplit", &
                             description="EM.TRCSpinComponentsOribalSplit", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCChannelsWFS_SkipEne", &
                             description="EM.TRCChannelsWFS_SkipEne", &
                             n_var=1, type_of_var=integer_t, default_i_val=100)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCChannelsWFS_SkipKP", &
                             description="EM.TRCChannelsWFS_SkipKP", &
                             n_var=1, type_of_var=integer_t, default_i_val=100)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TMatrixPDOSNWrite", &
                             description="EM.TMatrixPDOSNWrite", &
                             n_var=1, type_of_var=integer_t, default_i_val=20)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TMatrixSetPhaseIndex", &
                             description="EM.TMatrixSetPhaseIndex", &
                             n_var=1, type_of_var=integer_t, default_i_val=-1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCMinChannelIndex", &
                             description="EM.TRCMinChannelIndex", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.TRCMaxChannelIndex", &
                             description="EM.TRCMaxChannelIndex", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.NSetHBoundaryFromLeadsL", &
                             description="EM.NSetHBoundaryFromLeadsL", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.NSetHBoundaryFromLeadsR", &
                             description="EM.NSetHBoundaryFromLeadsR", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxOrbital1", &
                             description="EM.CurrentFluxOrbital1", &
                             n_var=1, type_of_var=integer_t, default_i_val=-1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxOrbital2", &
                             description="EM.CurrentFluxOrbital2", &
                             n_var=1, type_of_var=integer_t, default_i_val=-1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.ParallelOverKNum", &
                             description="EM.ParallelOverKNum", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="EM.LDOSLeadsProjection", &
                             description="EM.LDOSLeadsProjection", &
                             n_var=1, type_of_var=integer_t, default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Method", &
                             description="Sigma.Method", &
                             n_var=1, type_of_var=integer_t, default_i_val=1)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! emtrans_options() real-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="BS.Tolerance", &
                             description="BS.Tolerance", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-05_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="BS.Minimum", &
                             description="BS.Minimum", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.HSTol", &
                             description="EM.HSTol", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-10_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.EImag", &
                             description="Sigma.EImag", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AM.DeltaInit", &
                             description="AM.DeltaInit", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-03_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AM.DeltaToDe", &
                             description="AM.DeltaToDe", &
                             n_var=1, type_of_var=real_t, default_r_val=1.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="AM.Criterium", &
                             description="AM.Criterium", &
                             n_var=1, type_of_var=real_t, default_r_val=2.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SetEmin", &
                             description="SetEmin", &
                             n_var=1, type_of_var=real_t, default_r_val=-1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="SetEmax", &
                             description="SetEmax", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Tolab", &
                             description="Sigma.Tolab", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-06_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.WeightRho", &
                             description="EM.WeightRho", &
                             n_var=1, type_of_var=real_t, default_r_val=0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RotateSpinLeadsLeft.Theta", &
                             description="EM.RotateSpinLeadsLeft.Theta", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RotateSpinLeadsLeft.Phi", &
                             description="EM.RotateSpinLeadsLeft.Phi", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RotateSpinLeadsRight.Theta", &
                             description="EM.RotateSpinLeadsRight.Theta", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.RotateSpinLeadsRight.Phi", &
                             description="EM.RotateSpinLeadsRight.Phi", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="DeltaTransmission", &
                             description="DeltaTransmission", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.Dkimag", &
                             description="Sigma.Dkimag", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-06_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.DSVDMax", &
                             description="Sigma.DSVDMax", &
                             n_var=1, type_of_var=real_t, default_r_val=5.0e-12_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.DSVDMin", &
                             description="Sigma.DSVDMin", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-15_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.DMax", &
                             description="Sigma.DMax", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-05_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.RNoise", &
                             description="Sigma.RNoise", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.SkipSVD", &
                             description="Sigma.SkipSVD", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e+03_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.SVDTolZero", &
                             description="Sigma.SVDTolZero", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0e-08_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.WideBandImag", &
                             description="Sigma.WideBandImag", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SpinTorqueKScale", &
                             description="EM.SpinTorqueKScale", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetHLRZeroRatioLeft", &
                             description="EM.SetHLRZeroRatioLeft", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.SetHLRZeroRatioRight", &
                             description="EM.SetHLRZeroRatioRight", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxFL_L", &
                             description="EM.CurrentFluxFL_L", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxFR_L", &
                             description="EM.CurrentFluxFR_L", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxFL_R", &
                             description="EM.CurrentFluxFL_R", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="EM.CurrentFluxFR_R", &
                             description="EM.CurrentFluxFR_R", &
                             n_var=1, type_of_var=real_t, default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="BS.Delta", &
                             description="BS.Delta", &
                             n_var=1, type_of_var=real_t, default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! emtrans_options() physical-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="Sigma.WarnInOutput", &
                             description="Sigma.WarnInOutput", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=1.0e-05_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="TRCDE", &
                             description="TRCDE", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ! interface physical-valued keywords
         CALL keyword_create(keyword, __LOCATION__, name="HartreeLeadsLeft", &
                             description="HartreeLeadsLeft", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
         CALL keyword_create(keyword, __LOCATION__, name="HartreeLeadsRight", &
                             description="HartreeLeadsRight", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Bohr", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="HartreeLeadsBottom", &
                             description="HartreeLeadsBottom", &
                             n_var=1, type_of_var=real_t, &
                             unit_str="Ry", default_r_val=0.0_dp*0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

   END SUBROUTINE create_dft_smeagol_section

END MODULE input_cp2k_smeagol
