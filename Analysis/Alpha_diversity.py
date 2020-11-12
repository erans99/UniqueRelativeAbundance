from skbio.diversity.alpha import shannon
import pandas as pd
import os
import numpy as np

basepath = '/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/Alpha-Diversity'
phenotypes_thresholds = os.path.join(basepath,'Figures - phenotype_ranges.csv')


def alpha_div(x, method):
    return method((x.values*10**7).astype(int))


def compute_alpha_div(M_df, method):
    alpha = M_df.apply(alpha_div, 1, args=(method, ))
    return alpha


def get_outliers_by_phenotype(phenotype_name, phenotypes_df):
    phenotype_threshold = pd.read_csv(phenotypes_thresholds, index_col=0)
    if phenotype_name in phenotype_threshold.index:
        if phenotype_threshold.loc[phenotype_name,'Selection'] == 'Log':
            log = True
            lower = 10**phenotype_threshold.loc[phenotype_name, 'Log lower']
            upper = 10**phenotype_threshold.loc[phenotype_name, 'Log upper']
        else:
            log = False
            lower = phenotype_threshold.loc[phenotype_name, 'Linear lower']
            upper = phenotype_threshold.loc[phenotype_name, 'Linear upper']
        is_outlier = (phenotypes_df[phenotype_name] < lower).values | (phenotypes_df[phenotype_name] > upper).values
    else:
        log = False
        is_outlier = [False]*len(phenotypes_df[phenotype_name])
    return is_outlier, log


def report_alpha_diverity():
    phenotypes_df = pd.read_csv(os.path.join(basepath, 'all_phenotypes_with_age_gender__il.csv')).set_index('client_id')
    MB_df = pd.read_csv(os.path.join(basepath, 'sgb__il.csv')).set_index('client_id')
    assert (MB_df.index == phenotypes_df.index).all()
    for m in [shannon]:
        alpha = compute_alpha_div(MB_df, m)
        print("Using %s" % m.func_name)
        alpha.to_csv(os.path.join(basepath, 'Alpha-Shannon__il.csv'))
        pheno_df = pd.DataFrame(index=phenotypes_df.index, columns=['alpha'] + phenotypes_df.columns.values)
        pheno_df['alpha'] = alpha[pheno_df.index]
        pheno_df = pheno_df.sort_values('alpha')
        for phenotype in phenotypes_df.columns:
            phenotype_data = phenotypes_df[phenotype].dropna()
            is_outlier, log = get_outliers_by_phenotype(phenotype, phenotypes_df)
            if np.sum(is_outlier) > 0:
                phenotypes_df.loc[is_outlier, phenotype] = np.nan
            phenotype_data = phenotype_data.dropna()
            if log:
                phenotype_data = np.log10(phenotype_data)
            pheno_df.loc[phenotype_data.index,phenotype] = phenotype_data
        pheno_df = pheno_df.dropna(axis=1, how='all')
        pheno_df = pheno_df.loc[:, pheno_df.std(axis=0) > 0]
        pheno_df.to_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__il.csv'))
    return pheno_df


if __name__ == "__main__":
    report_alpha_diverity()
