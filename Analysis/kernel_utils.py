import numpy as np
import pandas as pd
import os.path
import scipy.linalg as la
import scipy.stats as stats
import warnings
from scipy.spatial.distance import pdist, squareform

pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)


# K - kernel matrix (2d numpy array)
# out_file - name of output file
# num_snps - #SNPs used to estimate kinship (not sure what gcta uses it for). Not relevant if bed_nold_file is specified.
# K_index - the names of the individuals in the matrix K. Only relevant if we want to use variable #SNPs for every entry
# bed_nold_file - name of Plink file with the SNPs that were used to estimate kinship. If not specified, the parameter num_snps will be used
def write_grm(K, out_file, num_snps=500000, K_index=None, bed_nold_file=None):
    # fill GRM columns 0,1,3
    n = K.shape[0]

    grm = np.zeros((n * (n - 1) / 2 + n, 4))
    tril_indices = np.tril_indices(n)
    grm[:, 0] = tril_indices[0] + 1
    grm[:, 1] = tril_indices[1] + 1
    grm[:, 2] = num_snps
    grm[:, 3] = K[tril_indices]

    # fill #non-missing SNPs columns
    if (bed_nold_file is not None):
        assert K_index is not None, 'K_index not provided'
        from pysnptools.snpreader.bed import Bed
        bed_nold = Bed(bed_nold_file, count_A1=True).read()
        notNan = (~np.isnan(bed_nold.val)).astype(np.float)
        notNan_K = notNan.dot(notNan.T)
        id2ind = dict([])
        for ind_i, ind in enumerate(bed_nold.iid[:, 1].astype(np.int)): id2ind[ind] = ind_i
        tril_indices_nold = [None, None]
        tril_indices_nold[0] = [id2ind[K_index[ind]] for ind in tril_indices[0]]
        tril_indices_nold[1] = [id2ind[K_index[ind]] for ind in tril_indices[1]]
        grm[:, 2] = notNan_K[tril_indices_nold]

    pd_grm = pd.DataFrame(grm, columns=['i', 'j', 'num_SNPs', 'K'])
    pd_grm['i'] = pd_grm['i'].astype(np.int)
    pd_grm['j'] = pd_grm['j'].astype(np.int)
    pd_grm.to_csv(out_file, compression='gzip', header=False, index=False, sep='\t', float_format='%0.6e')


def D2K(D):
    n = D.shape[0]
    centerM = np.eye(n) - 1.0 / n
    K = -0.5 * centerM.dot(D * D).dot(centerM)
    s, U = la.eigh(K)
    K = U.dot(U.T * np.abs(s[:, np.newaxis]))
    return K


def BrayCurtis(X_orig, is_log_abundance=True, zero_min_value=True):
    if is_log_abundance:
        X = 10 ** X_orig
    else:
        X = X_orig.copy()
    if zero_min_value: X[X_orig == np.min(X_orig)] = 0

    ####X = np.array([[1,2,3,0], [3,2,1,0], [0,0,0,1]])
    D = squareform(pdist(X, metric='braycurtis'))
    return D


def e_matrix(distance_matrix):
    """Compute E matrix from a distance matrix.
    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998."""
    return distance_matrix * distance_matrix / -2.0


def f_matrix(E_matrix):
    """Compute F matrix from E matrix.
    Centring step: for each element, the mean of the corresponding
    row and column are substracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
    row_means = E_matrix.mean(axis=1, keepdims=True)
    col_means = E_matrix.mean(axis=0, keepdims=True)
    matrix_mean = E_matrix.mean()
    return E_matrix - row_means - col_means + matrix_mean


def PCoA(D, suppress_warning=False, return_explained_variance=False):
    E_matrix = e_matrix(D)
    F_matrix = f_matrix(E_matrix)
    eigvals, eigvecs = la.eigh(F_matrix)
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0
    if (np.any(eigvals < 0) and not suppress_warning):
        warnings.warn(
            "The result contains negative eigenvalues."
            " Please compare their magnitude with the magnitude of some"
            " of the largest positive eigenvalues. If the negative ones"
            " are smaller, it's probably safe to ignore them, but if they"
            " are large in magnitude, the results won't be useful. See the"
            " Notes section for more details. The smallest eigenvalue is"
            " {0} and the largest is {1}.".format(eigvals.min(),
                                                  eigvals.max()),
            RuntimeWarning
        )

    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]
    num_positive = (eigvals >= 0).sum()
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)
    coordinates = eigvecs * np.sqrt(eigvals)

    if return_explained_variance: return coordinates, eigvals / eigvals.sum()
    return coordinates


def spearman_mat_vec(D, y):
    Dy = squareform(pdist(np.row_stack(y)))
    return spearman_mat_mat(D, Dy)


def spearman_mat_mat(D1, D2):
    try:
        D1 = D1.values
    except:
        pass
    try:
        D2 = D2.values
    except:
        pass
    return stats.spearmanr(D1[np.triu_indices(D1.shape[0], 1)], D2[np.triu_indices(D2.shape[0], 1)])[0]


def pearson_mat_vec(D, y):
    Dy = squareform(pdist(np.row_stack(y)))
    return pearson_mat_mat(D, Dy)


def pearson_mat_mat(D1, D2):
    try:
        D1 = D1.values
    except:
        pass
    try:
        D2 = D2.values
    except:
        pass
    return stats.pearsonr(D1[np.triu_indices(D1.shape[0], 1)], D2[np.triu_indices(D2.shape[0], 1)])[0]


