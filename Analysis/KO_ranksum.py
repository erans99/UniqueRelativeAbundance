import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection


# rel_p - pathways/modules lists
# rel_kos_df  - ko's X pathways/modules matrix
#                        (1- ko is in the pathway/module; 0 - ko is not in pathway/module)
# temp_K0_res - a matrix that shows ko's rank of p-value associated with a phenotype
#               (i.e association between phenotype and a ko)
def test_enrichment(temp_K0_res, rel_kos_df, rel_p):
    ko_enrichment_res = pd.DataFrame(index=rel_p, columns=['pval', 'natural_pval', 'fdr', 'pathway',
                                                           'in_fraction_of_KOs'])

    for pathway in rel_p:
        _, p0 = mannwhitneyu(temp_K0_res.loc[rel_kos_df[rel_kos_df[pathway] == 1].index]['rank'],
                             temp_K0_res.loc[rel_kos_df[rel_kos_df[pathway] == 0].index]['rank'],
                             alternative='less')
        _, p1 = mannwhitneyu(temp_K0_res.loc[rel_kos_df[rel_kos_df[pathway] == 1].index]['rank'],
                             temp_K0_res.loc[rel_kos_df[rel_kos_df[pathway] == 0].index]['rank'],
                             alternative='greater')
        if p0 < p1:

            ko_enrichment_res.loc[pathway, 'natural_pval'] = p0
            ko_enrichment_res.loc[pathway, 'pval'] = -(np.log10(p0))
        else:
            ko_enrichment_res.loc[pathway, 'natural_pval'] = p1
            ko_enrichment_res.loc[pathway, 'pval'] = np.log10(p1)

        ko_enrichment_res.loc[pathway, 'pathway'] = [pathway]
        ko_enrichment_res.loc[pathway, 'in_fraction_of_KOs'] = rel_kos_df[pathway].sum() / rel_kos_df.shape[0]
        ko_enrichment_res['fdr'] = fdrcorrection(ko_enrichment_res['natural_pval'])[1]
    ko_enrichment_res = ko_enrichment_res.sort_values('fdr')
    return ko_enrichment_res
