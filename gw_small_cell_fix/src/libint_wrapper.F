!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE libint_wrapper

#if(__LIBINT)
#include <libint2/libint2_params.h>
#include <libint2/config.h>
#endif

! maximum angular momentum to be supported in CP2K-LIBINT interface

   USE ISO_C_BINDING, ONLY: C_F_POINTER, &
                            C_F_PROCPOINTER, &
                            C_NULL_PTR, &
                            C_FUNPTR
   USE kinds, ONLY: dp
#if(__LIBINT)
   USE libint_f, ONLY: &
      libint2_build, libint2_build_eri, libint2_build_eri1, libint2_cleanup_eri, &
      libint2_cleanup_eri1, libint2_init_eri, libint2_init_eri1, libint2_static_cleanup, &
      libint2_static_init, libint_t, libint2_max_am_eri, libint2_init_3eri, libint2_cleanup_3eri, &
      libint2_init_2eri, libint2_cleanup_2eri, &
      libint2_build_2eri, libint2_build_3eri, libint2_build_3eri1, libint2_cleanup_3eri1, libint2_init_3eri1, &
      libint2_build_2eri1, libint2_cleanup_2eri1, libint2_init_2eri1
#endif
   USE orbital_pointers, ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cp_libint_t, prim_data_f_size, build_eri_size, build_deriv1_eri_size, &
             libint_max_am, libderiv_max_am1, cp_libint_get_eris, cp_libint_get_derivs, &
             cp_libint_init_eri, cp_libint_init_eri1, cp_libint_cleanup_eri, &
             cp_libint_cleanup_eri1, cp_libint_static_init, cp_libint_static_cleanup, &
             get_ssss_f_val, cp_libint_set_contrdepth, cp_libint_set_params_eri_screen, &
             cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, &
             cp_libint_init_3eri, cp_libint_cleanup_3eri, cp_libint_get_3eris, &
             cp_libint_init_2eri, cp_libint_cleanup_2eri, cp_libint_get_2eris, &
             cp_libint_get_3eri_derivs, cp_libint_init_3eri1, cp_libint_cleanup_3eri1, &
             cp_libint_get_2eri_derivs, cp_libint_init_2eri1, cp_libint_cleanup_2eri1

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_wrapper'

#if(__LIBINT)
   INTEGER, PARAMETER :: libint_max_am = libint2_max_am_eri
#else
   INTEGER, PARAMETER :: libint_max_am = 0
#endif

   INTEGER, PARAMETER :: libderiv_max_am1 = libint_max_am
   INTEGER, PARAMETER :: prim_data_f_size = 4*(libint_max_am) + 1
   INTEGER, PARAMETER :: libint_vrr_classes_size = 2*(libint_max_am) + 1
   INTEGER, PARAMETER :: libint_dvrr_classes_size = 2*(libderiv_max_am1) + 1
   INTEGER, PARAMETER :: build_eri_size = libint_max_am
   INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1

   TYPE :: cp_libint_t
      PRIVATE
#if(__LIBINT)
      TYPE(libint_t), DIMENSION(1) :: prv
#else
      INTEGER :: unused = -1
#endif
   END TYPE cp_libint_t

CONTAINS

   SUBROUTINE cp_libint_set_params_eri_screen(libint, A, B, C, D, P, Q, W, ZetaInv, EtaInv, ZetapEtaInv, Rho, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      INTEGER, INTENT(IN)                     :: m_max
      REAL(KIND=dp), DIMENSION(:)             :: F

#if(__LIBINT)
      libint%prv(1)%AB_x(1) = A(1) - B(1)
      libint%prv(1)%AB_y(1) = A(2) - B(2)
      libint%prv(1)%AB_z(1) = A(3) - B(3)

      libint%prv(1)%CD_x(1) = C(1) - D(1)
      libint%prv(1)%CD_y(1) = C(2) - D(2)
      libint%prv(1)%CD_z(1) = C(3) - D(3)

      libint%prv(1)%PA_x(1) = P(1) - A(1)
      libint%prv(1)%PA_y(1) = P(2) - A(2)
      libint%prv(1)%PA_z(1) = P(3) - A(3)

      libint%prv(1)%QC_x(1) = Q(1) - C(1)
      libint%prv(1)%QC_y(1) = Q(2) - C(2)
      libint%prv(1)%QC_z(1) = Q(3) - C(3)

      libint%prv(1)%WP_x(1) = W(1) - P(1)
      libint%prv(1)%WP_y(1) = W(2) - P(2)
      libint%prv(1)%WP_z(1) = W(3) - P(3)

      libint%prv(1)%WQ_x(1) = W(1) - Q(1)
      libint%prv(1)%WQ_y(1) = W(2) - Q(2)
      libint%prv(1)%WQ_z(1) = W(3) - Q(3)

      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#if 4*LIBINT2_MAX_AM_eri > 0 - 1
         IF (0 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0 (1) &
            = F(0+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 1 - 1
         IF (1 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1 (1) &
            = F(1+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 2 - 1
         IF (2 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2 (1) &
            = F(2+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 3 - 1
         IF (3 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3 (1) &
            = F(3+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 4 - 1
         IF (4 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4 (1) &
            = F(4+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 5 - 1
         IF (5 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5 (1) &
            = F(5+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 6 - 1
         IF (6 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6 (1) &
            = F(6+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 7 - 1
         IF (7 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7 (1) &
            = F(7+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 8 - 1
         IF (8 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8 (1) &
            = F(8+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 9 - 1
         IF (9 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9 (1) &
            = F(9+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 10 - 1
         IF (10 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10 (1) &
            = F(10+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 11 - 1
         IF (11 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11 (1) &
            = F(11+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 12 - 1
         IF (12 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12 (1) &
            = F(12+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 13 - 1
         IF (13 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13 (1) &
            = F(13+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 14 - 1
         IF (14 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14 (1) &
            = F(14+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 15 - 1
         IF (15 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15 (1) &
            = F(15+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 16 - 1
         IF (16 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16 (1) &
            = F(16+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 17 - 1
         IF (17 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17 (1) &
            = F(17+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 18 - 1
         IF (18 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18 (1) &
            = F(18+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 19 - 1
         IF (19 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19 (1) &
            = F(19+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 20 - 1
         IF (20 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20 (1) &
            = F(20+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 21 - 1
         IF (21 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_21 (1) &
            = F(21+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 22 - 1
         IF (22 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_22 (1) &
            = F(22+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 23 - 1
         IF (23 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_23 (1) &
            = F(23+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 24 - 1
         IF (24 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_24 (1) &
            = F(24+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 25 - 1
         IF (25 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_25 (1) &
            = F(25+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 26 - 1
         IF (26 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_26 (1) &
            = F(26+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 27 - 1
         IF (27 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_27 (1) &
            = F(27+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 28 - 1
         IF (28 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_28 (1) &
            = F(28+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 29 - 1
         IF (29 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_29 (1) &
            = F(29+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 30 - 1
         IF (30 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_30 (1) &
            = F(30+1)

#endif
#if 4*LIBINT2_MAX_AM_eri > 31 - 1
         IF (31 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_31 (1) &
            = F(31+1)

#endif

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)

      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_set_params_eri_screen

   SUBROUTINE cp_libint_set_params_eri_deriv(libint, A, B, C, D, P, Q, W, zeta_A, zeta_B, zeta_C, zeta_D, &
                                             ZetaInv, EtaInv, ZetapEtaInv, Rho, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: zeta_A, zeta_B, zeta_C, zeta_D, ZetaInv, EtaInv, ZetapEtaInv, Rho

      INTEGER, INTENT(IN)                     :: m_max
      REAL(KIND=dp), DIMENSION(:)             :: F

#if(__LIBINT)
      REAL(KIND=dp)                           :: gammap, gammaq, gammapq, rhop, rhoq
      libint%prv(1)%AB_x(1) = A(1) - B(1)
      libint%prv(1)%AB_y(1) = A(2) - B(2)
      libint%prv(1)%AB_z(1) = A(3) - B(3)

      libint%prv(1)%CD_x(1) = C(1) - D(1)
      libint%prv(1)%CD_y(1) = C(2) - D(2)
      libint%prv(1)%CD_z(1) = C(3) - D(3)

      libint%prv(1)%PA_x(1) = P(1) - A(1)
      libint%prv(1)%PA_y(1) = P(2) - A(2)
      libint%prv(1)%PA_z(1) = P(3) - A(3)

      libint%prv(1)%PB_x(1) = P(1) - B(1)
      libint%prv(1)%PB_y(1) = P(2) - B(2)
      libint%prv(1)%PB_z(1) = P(3) - B(3)

      libint%prv(1)%QC_x(1) = Q(1) - C(1)
      libint%prv(1)%QC_y(1) = Q(2) - C(2)
      libint%prv(1)%QC_z(1) = Q(3) - C(3)

      libint%prv(1)%WP_x(1) = W(1) - P(1)
      libint%prv(1)%WP_y(1) = W(2) - P(2)
      libint%prv(1)%WP_z(1) = W(3) - P(3)

      libint%prv(1)%WQ_x(1) = W(1) - Q(1)
      libint%prv(1)%WQ_y(1) = W(2) - Q(2)
      libint%prv(1)%WQ_z(1) = W(3) - Q(3)

      libint%prv(1)%two_alpha0_bra(1) = 2.0_dp*Zeta_A
      libint%prv(1)%two_alpha0_ket(1) = 2.0_dp*Zeta_B
      libint%prv(1)%two_alpha1_ket(1) = 2.0_dp*Zeta_D

      gammap = Zeta_A + Zeta_B
      gammaq = Zeta_C + Zeta_D
      gammapq = gammap*gammaq/(gammap + gammaq)
      libint%prv(1)%alpha1_rho_over_zeta2(1) = Zeta_A*gammapq/(gammap*gammap)
      libint%prv(1)%alpha2_rho_over_zeta2(1) = Zeta_B*gammapq/(gammap*gammap)
      libint%prv(1)%alpha4_rho_over_eta2(1) = Zeta_D*gammapq/(gammaq*gammaq)
      libint%prv(1)%alpha1_over_zetapluseta(1) = Zeta_A/(gammap + gammaq)
      libint%prv(1)%alpha2_over_zetapluseta(1) = Zeta_B/(gammap + gammaq)
      libint%prv(1)%alpha4_over_zetapluseta(1) = Zeta_D/(gammap + gammaq)

      rhop = Zeta_A*Zeta_B/gammap
      rhoq = Zeta_C*Zeta_D/gammaq
      libint%prv(1)%rho12_over_alpha1(1) = rhop/Zeta_A

      libint%prv(1)%rho34_over_alpha3(1) = rhoq/Zeta_C

      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#if 4*LIBINT2_MAX_AM_eri > 0 - 1
         IF (0 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(0+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 1 - 1
         IF (1 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(1+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 2 - 1
         IF (2 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(2+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 3 - 1
         IF (3 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(3+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 4 - 1
         IF (4 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(4+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 5 - 1
         IF (5 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(5+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 6 - 1
         IF (6 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(6+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 7 - 1
         IF (7 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(7+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 8 - 1
         IF (8 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(8+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 9 - 1
         IF (9 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(9+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 10 - 1
         IF (10 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(10+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 11 - 1
         IF (11 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(11+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 12 - 1
         IF (12 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(12+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 13 - 1
         IF (13 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(13+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 14 - 1
         IF (14 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(14+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 15 - 1
         IF (15 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(15+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 16 - 1
         IF (16 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(16+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 17 - 1
         IF (17 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(17+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 18 - 1
         IF (18 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(18+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 19 - 1
         IF (19 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(19+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 20 - 1
         IF (20 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(20+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 21 - 1
         IF (21 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_21 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(21+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 22 - 1
         IF (22 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_22 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(22+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 23 - 1
         IF (23 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_23 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(23+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 24 - 1
         IF (24 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_24 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(24+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 25 - 1
         IF (25 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_25 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(25+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 26 - 1
         IF (26 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_26 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(26+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 27 - 1
         IF (27 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_27 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(27+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 28 - 1
         IF (28 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_28 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(28+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 29 - 1
         IF (29 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_29 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(29+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 30 - 1
         IF (30 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_30 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(30+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 31 - 1
         IF (31 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_31 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(31+1)
#endif

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(zeta_A)
      MARK_USED(zeta_B)
      MARK_USED(zeta_C)
      MARK_USED(zeta_D)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_set_params_eri_deriv

   SUBROUTINE cp_libint_set_params_eri(libint, A, B, C, D, ZetaInv, EtaInv, ZetapEtaInv, Rho, P, Q, W, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      REAL(KIND=dp), DIMENSION(:)             :: F

      INTEGER, INTENT(IN)                     :: m_max

#if(__LIBINT)
      libint%prv(1)%AB_x(1) = A(1) - B(1)
      libint%prv(1)%AB_y(1) = A(2) - B(2)
      libint%prv(1)%AB_z(1) = A(3) - B(3)

      libint%prv(1)%CD_x(1) = C(1) - D(1)
      libint%prv(1)%CD_y(1) = C(2) - D(2)
      libint%prv(1)%CD_z(1) = C(3) - D(3)

      libint%prv(1)%PA_x(1) = P(1) - A(1)
      libint%prv(1)%PA_y(1) = P(2) - A(2)
      libint%prv(1)%PA_z(1) = P(3) - A(3)

      libint%prv(1)%QC_x(1) = Q(1) - C(1)
      libint%prv(1)%QC_y(1) = Q(2) - C(2)
      libint%prv(1)%QC_z(1) = Q(3) - C(3)

      libint%prv(1)%WP_x(1) = W(1) - P(1)
      libint%prv(1)%WP_y(1) = W(2) - P(2)
      libint%prv(1)%WP_z(1) = W(3) - P(3)

      libint%prv(1)%WQ_x(1) = W(1) - Q(1)
      libint%prv(1)%WQ_y(1) = W(2) - Q(2)
      libint%prv(1)%WQ_z(1) = W(3) - Q(3)

      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#if 4*LIBINT2_MAX_AM_eri > 0 - 1

         IF (0 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(0+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 1 - 1

         IF (1 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(1+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 2 - 1

         IF (2 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(2+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 3 - 1

         IF (3 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(3+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 4 - 1

         IF (4 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_4 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(4+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 5 - 1

         IF (5 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_5 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(5+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 6 - 1

         IF (6 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_6 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(6+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 7 - 1

         IF (7 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_7 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(7+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 8 - 1

         IF (8 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_8 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(8+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 9 - 1

         IF (9 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_9 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(9+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 10 - 1

         IF (10 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_10 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(10+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 11 - 1

         IF (11 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_11 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(11+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 12 - 1

         IF (12 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_12 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(12+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 13 - 1

         IF (13 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_13 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(13+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 14 - 1

         IF (14 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_14 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(14+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 15 - 1

         IF (15 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_15 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(15+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 16 - 1

         IF (16 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_16 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(16+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 17 - 1

         IF (17 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_17 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(17+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 18 - 1

         IF (18 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(18+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 19 - 1

         IF (19 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_19 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(19+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 20 - 1

         IF (20 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_20 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(20+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 21 - 1

         IF (21 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_21 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(21+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 22 - 1

         IF (22 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_22 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(22+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 23 - 1

         IF (23 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_23 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(23+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 24 - 1

         IF (24 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_24 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(24+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 25 - 1

         IF (25 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_25 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(25+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 26 - 1

         IF (26 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_26 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(26+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 27 - 1

         IF (27 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_27 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(27+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 28 - 1

         IF (28 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_28 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(28+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 29 - 1

         IF (29 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_29 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(29+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 30 - 1

         IF (30 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_30 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(30+1)
#endif
#if 4*LIBINT2_MAX_AM_eri > 31 - 1

         IF (31 <= m_max) &
            libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_31 (1) & ! ERROR: __LIBINT_MAX_AM is too large
            = F(31+1)
#endif

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_set_params_eri
! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_eris(n_d, n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_eri(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_eris

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eris(n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_3eri(n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eris

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eri_derivs(n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      INTEGER                                            :: a_mysize(1)
      REAL(dp), DIMENSION(:, :), POINTER                    :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER :: i

      CALL C_F_PROCPOINTER(libint2_build_3eri1(n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      ALLOCATE (p_work(a_mysize(1), 9))

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 9
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=a_mysize)
         p_work(:, i) = p_work_tmp
      END DO
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eri_derivs

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eri_derivs(n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      INTEGER                                            :: a_mysize(1)
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER :: i

      CALL C_F_PROCPOINTER(libint2_build_2eri1(n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      ALLOCATE (p_work(a_mysize(1), 6))

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 6
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=a_mysize)
         p_work(:, i) = p_work_tmp
      END DO
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_2eri_derivs

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eris(n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_2eri(n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_2eris

! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param work_forces ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_derivs(n_d, n_c, n_b, n_a, lib, work_forces, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco &
                          (n_d), 12)                                      :: work_forces
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER                                            :: i, k
#endif

#if(__LIBINT)
      CALL C_F_PROCPOINTER(libint2_build_eri1(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      DO k = 1, 12
         IF (k == 4 .OR. k == 5 .OR. k == 6) CYCLE
         CALL C_F_POINTER(lib%prv(1)%targets(k), p_work, SHAPE=a_mysize)
         DO i = 1, a_mysize(1)
            work_forces(i, k) = p_work(i)
         END DO
      END DO
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(work_forces)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_derivs

   FUNCTION get_ssss_f_val(lib)
      TYPE(cp_libint_t) :: lib
      REAL(KIND=dp) :: get_ssss_f_val

      get_ssss_f_val = 0

#if(__LIBINT)
      get_ssss_f_val = lib%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0(1)
#else
      MARK_USED(lib)
      get_ssss_f_val = 0.0_dp
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END FUNCTION get_ssss_f_val

   SUBROUTINE cp_libint_init_eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_eri

   SUBROUTINE cp_libint_set_contrdepth(lib, contrdepth)
      TYPE(cp_libint_t) :: lib
      INTEGER           :: contrdepth
#if(__LIBINT)
#if LIBINT_CONTRACTED_INTS
      lib%prv(1)%contrdepth = contrdepth
#else
      MARK_USED(lib)
      MARK_USED(contrdepth)
#endif
#else
      MARK_USED(lib)
      MARK_USED(contrdepth)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_set_contrdepth

   SUBROUTINE cp_libint_init_eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_eri1

   SUBROUTINE cp_libint_init_3eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_3eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_3eri

   SUBROUTINE cp_libint_init_3eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_3eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_3eri1

   SUBROUTINE cp_libint_init_2eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_2eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_2eri1

   SUBROUTINE cp_libint_init_2eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_2eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_2eri

   SUBROUTINE cp_libint_cleanup_eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_eri

   SUBROUTINE cp_libint_cleanup_eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_eri1

   SUBROUTINE cp_libint_cleanup_3eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_3eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_3eri

   SUBROUTINE cp_libint_cleanup_3eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_3eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_3eri1

   SUBROUTINE cp_libint_cleanup_2eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_2eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_2eri1

   SUBROUTINE cp_libint_cleanup_2eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_2eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_2eri

   SUBROUTINE cp_libint_static_init()
#if(__LIBINT)
      CALL libint2_static_init()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_static_init

   SUBROUTINE cp_libint_static_cleanup()
#if(__LIBINT)
      CALL libint2_static_cleanup()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_static_cleanup

END MODULE libint_wrapper
