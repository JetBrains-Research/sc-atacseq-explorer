#!/usr/bin/env python
#
# This is a copy of https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/analysis/diffexp.py
# Adapted for Python 3.5 by JetBrains Research
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp
import scipy.stats
from sklearn.utils import sparsefuncs
import sys

SSEQ_ZETA_QUANTILE = 0.995


def estimate_size_factors(x):
    """ Estimate size factors (related to cell RNA content and GEM-to-GEM technical variance)
    Args:
      x - Sparse matrix (csc) of counts (feature x cell)
    Returns:
      Array of floats, one per cell.
    """
    counts_per_cell = np.squeeze(np.asarray(x.sum(axis=0)))
    size_factors = counts_per_cell.astype(np.float64) / np.median(counts_per_cell)
    return size_factors


def compute_sseq_params(x, zeta_quantile=SSEQ_ZETA_QUANTILE):
    """ Compute global parameters for the sSeq differential expression method.
    The key parameters are the shrunken feature-wise dispersions.

    This method was published in:
    Yu D, et al. (2013) Shrinkage estimation of dispersion in Negative Binomial models for RNA-seq experiments with small sample size.
    Bioinformatics. 29: 1275-1282. doi: 10.1093/bioinformatics/btt143
    Args:
      x - Sparse matrix (csc) of counts (feature x cell)
      zeta_quantile (float) - Quantile of method-of-moments dispersion estimates to
                              use as the shrinkage target zeta.
    Returns:
      A dictionary containing the sSeq parameters and some diagnostic info.
    """
    # Number of cells
    N = x.shape[1]

    # Number of features
    G = x.shape[0]

    # Estimate size factors and normalize the matrix for quick mean/var calcs
    size_factors = estimate_size_factors(x)
    # Cast to float to prevent truncation of 1 -> 0 for size factors < 1
    x_norm = scipy.sparse.csc_matrix(x, dtype=np.float64, copy=True)
    sparsefuncs.inplace_column_scale(x_norm, 1.0 / size_factors)

    # Estimate featurewise mean, variance, and dispersion by the method of moments
    # assuming that each feature follows a negative-binomial distribution.
    mean_g = np.squeeze(np.asarray(x_norm.mean(axis=1, dtype=np.float64)))
    # V[X] = E[X^2] - E[X]^2
    mean_sq_g = np.squeeze(np.asarray(x_norm.multiply(x_norm).mean(axis=1, dtype=np.float64)))
    var_g = mean_sq_g - np.square(mean_g)

    # Method of moments estimate of feature-wise dispersion (phi)
    # Only use features with non-zero variance in the following estimation
    use_g = var_g > 0
    phi_mm_g = np.zeros(G)
    phi_mm_g[use_g] = np.maximum(0, (float(N) * var_g[use_g] - mean_g[use_g] * np.sum(1.0 / size_factors)) /
                                 (np.square(mean_g[use_g]) * np.sum(1.0 / size_factors)))

    # Estimate the optimal global target dispersion (zeta_hat).
    # The true optimal zeta is that which minimizes the MSE vs the true dispersions.
    # The featurewise dispersions will be "shrunk" towards our estimate of zeta.

    # Use a high quantile of the MoM dispersion as our shrinkage target
    # per the rule of thumb in Yu, et al.
    zeta_hat = np.nanpercentile(phi_mm_g[use_g], 100.0 * zeta_quantile)

    # Compute delta, the optimal shrinkage towards zeta_hat
    # This defines a linear function that shrinks the MoM dispersion estimates
    mean_phi_mm_g = np.mean(phi_mm_g[use_g])
    delta = (np.sum(np.square(phi_mm_g[use_g] - mean_phi_mm_g)) / float(G - 1)) / \
            (np.sum(np.square(phi_mm_g[use_g] - zeta_hat)) / float(G - 2))

    # Compute the shrunken dispersion estimates
    # Interpolate between the MoM estimates and zeta_hat by delta
    phi_g = np.full(G, np.nan)
    if np.any(phi_mm_g[use_g] > 0):
        phi_g[use_g] = (1 - delta) * phi_mm_g[use_g] + delta * zeta_hat
    else:
        phi_g[use_g] = 0.0

    return {
        'N': N,
        'G': G,
        'size_factors': size_factors,
        'mean_g': mean_g,
        'var_g': var_g,
        'use_g': use_g,
        'phi_mm_g': phi_mm_g,
        'eval_zeta': None,
        'eval_asd': None,
        'asd_slope': None,
        'zeta_hat': zeta_hat,
        'delta': delta,
        'phi_g': phi_g,
    }


def neg_bin_log_pmf(k, mu, phi):
    """ Log(PMF) of negative binomial distribution with mean mu and dispersion phi,
    conveniently parameterized.
    Args:
      k (int) - NB random variable
      mu (float) - mean
      phi (float) - dispersion
    Returns:
      The log of the pmf at k. """
    r = 1.0 / phi
    return gammaln(r + k) - (gammaln(r) + gammaln(k + 1)) + \
           k * np.log(mu / (r + mu)) + \
           r * np.log(r / (r + mu))


def nb_exact_test(x_a, x_b, size_factor_a, size_factor_b, mu, phi):
    """ Compute p-value for a pairwise exact test using the negative binomial.
    Args:
      x_a (int) - Total count for a single feature in group A
      x_b (int) - Total count for a single feature in group B
      size_factor_a (float) - Sum of size factors for group A
      size_factor_b (float) - Sum of size factors for group B
      mu (float) - Common mean count for this feature
      phi (float) - Common dispersion for this feature
    Returns:
      p-value (float); the probability that a random pair of counts under the null hypothesis is more extreme
      than the observed pair of counts. """
    size_factor_a = float(size_factor_a)
    size_factor_b = float(size_factor_b)
    mu = float(mu)
    phi = float(phi)

    if (x_a + x_b) == 0:
        return 1.0
    if phi == 0:
        return 1.0
    if size_factor_a == 0 or size_factor_b == 0:
        return 1.0

    all_x_a = np.arange(0, 1 + x_a + x_b, 1)
    all_x_b = np.arange(x_a + x_b, -1, -1)

    def log_prob(x, size_factor):
        return neg_bin_log_pmf(x, size_factor * mu, phi / size_factor)

    log_p_obs = log_prob(x_a, size_factor_a) + log_prob(x_b, size_factor_b)
    log_p_all = log_prob(all_x_a, size_factor_a) + log_prob(all_x_b, size_factor_b)

    more_extreme = log_p_all <= log_p_obs
    if np.sum(more_extreme) == 0:
        return 0.0

    return np.exp(logsumexp(log_p_all[log_p_all <= log_p_obs]) - logsumexp(log_p_all))


def nb_asymptotic_test(x_a, x_b, size_factor_a, size_factor_b, mu, phi):
    """ Compute p-value for a pairwise exact test using a fast beta approximation
    to the conditional joint distribution of (x_a, x_b).
    Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, 321-332
    "It is based a method-of-moments gamma approximation to the negative binomial distribution."
      - Personal communication w/ author
    Adapted from implementation in the "edgeR" package:
      https://github.com/Bioconductor-mirror/edgeR/blob/1ab290c9585335cf99bb41f50cfce2ce4d40f907/R/exactTestBetaApprox.R
    This function is vectorized. It always returns a vector even if the inputs are scalar.
    Args:
      x_a (int/np.array) - Total count for a single feature in group A
      x_b (int/np.array) - Total count for a single feature in group B
      size_factor_a (float) - Sum of size factors for group A
      size_factor_b (float) - Sum of size factors for group B
      mu (float/np.array) - Common mean count for this feature
      phi (float/np.array) - Common dispersion for this feature
    Returns:
      p-value (np.array); the probability that a random pair of counts under the null hypothesis is more extreme
      than the observed pair of counts. """
    x_a = np.array(x_a, ndmin=1, copy=False)
    x_b = np.array(x_b, ndmin=1, copy=False)
    mu = np.array(mu, ndmin=1, copy=False)
    phi = np.array(phi, ndmin=1, copy=False)

    alpha = size_factor_a * mu / (1 + phi * mu)
    beta = (size_factor_b / size_factor_a) * alpha

    total = x_a + x_b

    median = scipy.stats.beta.median(alpha, beta)
    left = ((x_a + 0.5) / total) < median
    right = np.logical_not(left)

    p = np.empty(len(x_a))
    p[left] = 2 * scipy.stats.beta.cdf((x_a[left] + 0.5) / total[left],
                                       alpha[left], beta[left])
    # If X ~ Beta(a, b) then 1 - X ~ Beta(b, a)
    # This avoids the asymmetry in beta.cdf
    p[right] = 2 * (scipy.stats.beta.cdf((x_b[right] + 0.5) / total[right],
                                         beta[right], alpha[right]))
    return p


def adjust_pvalue_bh(p):
    """ Multiple testing correction of p-values using the Benjamini-Hochberg procedure """
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return q[np.argsort(descending)]


def numpy_logical_and_list(list_of_logicals):
    assert (len(list_of_logicals) >= 2)
    output = list_of_logicals[0]
    for i in range(1, len(list_of_logicals)):
        output = np.logical_and(output, list_of_logicals[i])
    return output


def sseq_differential_expression(x, cond_a, cond_b, sseq_params, big_count=900):
    """ Run sSeq pairwise differential expression test.
      Args:
        x - Sparse matrix (csc) of counts (feature x cell)
        cond_a (np.array(int)): Indices of cells in group A
        cond_b (np.array(int)): Indices of cells in group B
        sseq_params (dict): Precomputed global parameters
        big_count (int): Use asymptotic approximation if both counts > this
      Returns:
        A pd.DataFrame with DE results for group A relative to group B """
    x_a = x[:, cond_a]
    x_b = x[:, cond_b]

    # Number of features
    G = x.shape[0]

    # Size factors
    size_factor_a = np.sum(sseq_params['size_factors'][cond_a])
    size_factor_b = np.sum(sseq_params['size_factors'][cond_b])

    # Compute p-value for each feature
    p_values = np.ones(G)

    feature_sums_a = np.squeeze(np.asarray(x_a.sum(axis=1)))
    feature_sums_b = np.squeeze(np.asarray(x_b.sum(axis=1)))

    big = numpy_logical_and_list([sseq_params['use_g'], feature_sums_a > big_count, feature_sums_b > big_count])
    small = np.logical_and(sseq_params['use_g'], np.logical_not(big))

    sys.stderr.write("Computing %d exact tests and %d asymptotic tests.\n" %
                     (np.sum(small), np.sum(big)))

    # Compute exact test for small-count features
    for i in np.flatnonzero(small):
        p_values[i] = nb_exact_test(feature_sums_a[i], feature_sums_b[i],
                                    size_factor_a, size_factor_b,
                                    sseq_params['mean_g'][i],
                                    sseq_params['phi_g'][i])
    # Compute asymptotic approximation for big-count features
    p_values[big] = nb_asymptotic_test(feature_sums_a[big],
                                       feature_sums_b[big],
                                       size_factor_a, size_factor_b,
                                       sseq_params['mean_g'][big],
                                       sseq_params['phi_g'][big])
    # Adjust p-values for multiple testing correction
    # Only adjust the features that were actually tested
    adj_p_values = p_values.copy()
    adj_p_values[sseq_params['use_g']] = adjust_pvalue_bh(p_values[sseq_params['use_g']])

    de_result = pd.DataFrame({
        'tested': sseq_params['use_g'],
        'sum_a': feature_sums_a,
        'sum_b': feature_sums_b,
        'common_mean': sseq_params['mean_g'],
        'common_dispersion': sseq_params['phi_g'],
        'norm_mean_a': feature_sums_a / size_factor_a,
        'norm_mean_b': feature_sums_b / size_factor_b,
        'p_value': p_values,
        'adjusted_p_value': adj_p_values,
        # Introduce a pseudocount into log2(fold_change)
        'log2_fold_change': np.log2((1 + feature_sums_a) / (1 + size_factor_a)) - \
                            np.log2((1 + feature_sums_b) / (1 + size_factor_b))
    })

    return de_result

