import numpy as np
import pandas as pd


pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)


# K - kernel matrix (2d numpy array)
# out_file - name of output file
# num_snps - #SNPs used to estimate kinship(not sure what gcta uses it for). Not relevant if bed_nold_file is specified.
# K_index - the names of the individuals in the matrix K. Only relevant if we want to use variable #SNPs for every entry
# bed_nold_file - name of Plink file with the SNPs that were used to estimate kinship. If not specified,
# the parameter num_snps will be used
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
    if bed_nold_file is not None:
        assert K_index is not None, 'K_index not provided'
        from pysnptools.snpreader.bed import Bed
        bed_nold = Bed(bed_nold_file, count_A1=True).read()
        notNan = (~np.isnan(bed_nold.val)).astype(np.float)
        notNan_K = notNan.dot(notNan.T)
        id2ind = dict([])
        for ind_i, ind in enumerate(bed_nold.iid[:, 1].astype(np.int)):
            id2ind[ind] = ind_i
        tril_indices_nold = [None, None]
        tril_indices_nold[0] = [id2ind[K_index[ind]] for ind in tril_indices[0]]
        tril_indices_nold[1] = [id2ind[K_index[ind]] for ind in tril_indices[1]]
        grm[:, 2] = notNan_K[tril_indices_nold]

    pd_grm = pd.DataFrame(grm, columns=['i', 'j', 'num_SNPs', 'K'])
    pd_grm['i'] = pd_grm['i'].astype(np.int)
    pd_grm['j'] = pd_grm['j'].astype(np.int)
    pd_grm.to_csv(out_file, compression='gzip', header=False, index=False, sep='\t', float_format='%0.6e')
