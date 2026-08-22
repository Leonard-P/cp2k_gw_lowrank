# CP2K Fork: Low-Rank G0W0

This is a fork of [CP2K](https://github.com/cp2k/cp2k). It was created for a Bachelor's thesis to assess the accuracy of a low-rank approximation (LRA) to the screened Coulomb interaction $W$ within the small-cell G0W0 code (`gw_small_cell_full_kp.F`).

## Branches

### `feature/gw_lowrank`
Applies a low-rank approximation to the matrices $W(\mathbf{k}, i\omega)$ independently for each k-point. 

The map $L$ is computed by diagonalizing the training sum $\sum_\omega (M^a W(i\omega) M^a)^2$ with a metric exponent $a$. Its remaining eigenvectors form the columns of the map $L$, and the compressed matrix is $w = L^\dagger W L$. To evaluate accuracy without modifying the self-energy code, this implementation immediately maps back: 
$W \approx M^{-1-a} L L^\dagger M^a W M^a L L^\dagger M^{-1-a}$

This code thus measures accuracy only, no computational speedups are realized. The findings suggest speedups are possible if $\Sigma$ is computed in k-space (as $a=0$ was found to be necessary for a good LRA, so that $M(k)^{-1}$ needs to be multiplied in $k$-space after the LRA).

Respective input parameters for the `&GW` section:
* `LOW_RANK_APPROX_W_EPS`: Threshold $\epsilon$ below which eigenvalues of the training sum are discarded, alternatively, use
* `LOW_RANK_APPROX_W_RANK`: Sets a fixed rank (discards the $N_\mathrm{RI} - \mathrm{rank}$ eigenvectors with the smallest values).
* `LOW_RANK_APPROX_M_EXPONENT`: Adjusts the metric exponent $a$ (e.g., `-0.5` or `0.0`).

### `feature/gw_lowrank_kindependent`
Analogous to the previous branch, but computes a single map $L$ using a training sum over all frequency points and *all* k-points. That is, one universal low-rank subspace for all matrices $W$.