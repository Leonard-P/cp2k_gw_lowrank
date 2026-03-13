/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/
#ifndef DBM_MULTIPLY_H
#define DBM_MULTIPLY_H

#include "dbm_matrix.h"

#include <stdbool.h>
#include <stdint.h>

/*******************************************************************************
 * \brief Performs a multiplication of two dbm_matrix_t matrices,
          as  C := alpha * op( A ) * op( B ) + beta * C.

          The filter_eps parameter is used for on-the-fly pre-screening:
          individual multiplications are skipped if the product of the frobenius
          norms of the left- and right-matrix blocks are less than the specified
          epsilon divided by the maximum number of possible multiplies in each
          row.

          The filter_eps_post parameter controls the final post-multiply
          filtering of C (blocks with frobenius norm below this threshold are
          removed). When set to a negative value, filter_eps is used instead
          (legacy behavior).

          The nblocks_per_row parameter, when non-NULL, provides externally
          computed global block-counts per row of op(A). This overrides the
          internal per-row counting + MPI reduction, allowing the caller to
          supply counts that span multiple sub-groups (e.g. TAS CASE 2).
          When NULL, the counts are computed internally as before.
 * \author Ole Schuett
 ******************************************************************************/
void dbm_multiply(const bool transa, const bool transb, const double alpha,
                  const dbm_matrix_t *matrix_a, const dbm_matrix_t *matrix_b,
                  const double beta, dbm_matrix_t *matrix_c,
                  const bool retain_sparsity, const double filter_eps,
                  const double filter_eps_post, const int *nblocks_per_row,
                  int64_t *flop);

#endif

// EOF
