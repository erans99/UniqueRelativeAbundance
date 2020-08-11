import numpy
import pandas
import os
import math
import glob
import sys
from scipy import stats
from multiprocessing import Pool
import matplotlib.pyplot as plt
import statsmodels.stats.multitest
import socket

ONLY_FULL = False #True

data_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/DataForRevision"
out_path = os.path.join(data_path, "Atlas")
if ONLY_FULL:
    out_path += "_FullLen"

if not os.path.exists(out_path):
    os.makedirs(out_path)

NUM_TH = 1
MIN_FILL = 0.0001


def get_top_occ(df_in, oc):
    df = df_in.copy()
    lst_species = []
    df[df == MIN_FILL] = numpy.nan
    for c in df_in.columns:
        if df[c].count() > oc * df.shape[0]:
            lst_species.append(c)

    return lst_species


def throw_outliers(df_meta, pheno, num_stds=4):
    if num_stds <= 0:
        print("WTF. Need to throw above positive number of STDs from mean")
    # calc median and std without outliers
    curr_005 = df_meta[pheno].quantile(0.05)
    curr_095 = df_meta[pheno].quantile(0.95)
    df_curr = df_meta[(df_meta[pheno] >= curr_005) & (df_meta[pheno] <= curr_095)]
    curr_median = df_curr[pheno].median()
    curr_stdev = df_curr[pheno].std()

    # take all data <= 4 stds from median
    return df_meta[(df_meta[pheno] >= (curr_median - num_stds * curr_stdev)) &
                   (df_meta[pheno] <= (curr_median + num_stds * curr_stdev))][pheno]

def get_spear_ttest(inp):
    df_ura, df_tar, pheno = inp
    df_res = pandas.DataFrame()
    dct = {}
    dct['coef'] = []
    dct['spear_pval'] = []
    dct['spear'] = []
    dct['pears_pval'] = []
    dct['pears'] = []
    dct['ttest_pval'] = []
    dct['ttest'] = []
    dct['species'] = []
    dct['occ'] = []
    dct['occ_zero_bac'] = []
    dct['occ_non_zero_bac'] = []
    dct['pheno'] = []

    for c in df_ura.columns:
        ura_zero = df_ura[df_ura[c] == MIN_FILL].index
        ura_non_zero = df_ura[df_ura[c] > MIN_FILL].index
        if len(ura_non_zero) > 0:
            spear = stats.spearmanr(list(df_ura[c]), list(df_tar))
            if (len(ura_non_zero) > 1) and (len(ura_zero) > 1):
                ttest = stats.ttest_ind(list(df_tar[ura_zero]), list(df_tar[ura_non_zero]), equal_var=False)
            else:
                ttest = [numpy.nan, numpy.nan]
            pears = stats.pearsonr(list(numpy.log(df_ura[c])), list(df_tar))
            x = numpy.log(df_ura[c])
            dct['coef'].append(((x - x.mean()) * (df_tar - df_tar.mean())).sum() /
                               ((x - x.mean()) ** 2).sum())
            if len(ura_zero) > 0:
                dct['ttest'].append(numpy.mean(list(df_tar[ura_zero])) - numpy.mean(list(df_tar[ura_non_zero])))
            else:
                dct['ttest'].append(numpy.nan)
        else:
            spear = [numpy.nan, numpy.nan]
            ttest = [numpy.nan, numpy.nan]
            pears = [numpy.nan, numpy.nan]
            dct['coef'].append(numpy.nan)
            dct['ttest'].append(numpy.nan)

        dct['occ_zero_bac'].append(len(ura_zero))
        dct['occ_non_zero_bac'].append(len(ura_non_zero))
        dct['spear_pval'].append(spear[1])
        dct['spear'].append(spear[0])

        dct['pears_pval'].append(pears[1])
        dct['pears'].append(pears[0])

        dct['ttest_pval'].append(ttest[1])
        dct['pheno'].append(pheno)
        dct['species'].append(c)
        dct['occ'].append(len(df_tar))

    for k in dct.keys():
        df_res[k] = dct[k]

    return df_res


def genes_spear(df_ura, df_pheno, pheno):
    lst = []
    step = math.ceil(len(df_ura.columns) / NUM_TH)
    if NUM_TH != 1:
        for i in range(0, len(df_ura.columns), step):
            curr_cols = df_ura.columns[i:min(i + step, len(df_ura.columns))]
            lst.append((df_ura[curr_cols], df_pheno, pheno))

        p = Pool(NUM_TH)

        res = p.map(get_spear_ttest, lst)
        p.close()

        spear = pandas.concat(res)
    else:
        spear = get_spear_ttest((df_ura, df_pheno, pheno))

    return spear


df = pandas.read_csv(os.path.join(data_path, 'ura_q_il.csv'), index_col=0)
df.index = ['D2_%d' % x for x in df.client_id]
for i, c in enumerate(df.columns):
    if not 'SGB' in c:
        break
df_ura_il = df[df.columns[:i]]
df_meta_il = df[df.columns[i:]]
df = pandas.read_csv(os.path.join(data_path, 'ura_q_us.csv'), index_col=0)
df.index = ['D2_%d' % x for x in df.client_id]
df_ura_us = df[df.columns[:i]]
df_meta_us = df[df.columns[i:]]
df = pandas.read_csv(os.path.join(data_path, 'ura_q_il_validation.csv'), index_col=0)
df.index = ['D2_%d' % x for x in df.client_id]
df_ura_il_val = df[df.columns[:i]]
df_meta_il_val = df[df.columns[i:]]

if ONLY_FULL:
    link = pandas.read_csv(os.path.join(data_path, 'link_clientID_to_LabData.csv'), index_col=1)

    good_inds = df_ura_il.index.intersection(link[((link.RawReadLength == 75) & (link.URAMapUsed == 8*10**6)) |
                                                ((link.RawReadLength == 100) & (link.URAMapUsed == 5*10**6))].index)
    print("For IL left with %d of %d" % (len(good_inds), len(df_ura_il)))
    df_ura_il = df_ura_il.loc[good_inds]
    df_meta_il = df_meta_il.loc[good_inds]
    good_inds = df_ura_us.index.intersection(link[((link.RawReadLength == 75) & (link.URAMapUsed == 8*10**6)) |
                                            ((link.RawReadLength == 100) & (link.URAMapUsed == 5*10**6))].index)
    print("For US left with %d of %d" % (len(good_inds), len(df_ura_us)))
    df_ura_us = df_ura_us.loc[good_inds]
    df_meta_us = df_meta_us.loc[good_inds]
    good_inds = df_ura_il_val.index.intersection(link[((link.RawReadLength == 75) & (link.URAMapUsed == 8*10**6)) |
                                                ((link.RawReadLength == 100) & (link.URAMapUsed == 5*10**6))].index)
    print("For ILval left with %d of %d" % (len(good_inds), len(df_ura_il_val)))
    df_ura_il_val = df_ura_il_val.loc[good_inds]
    df_meta_il_val = df_meta_il_val.loc[good_inds]

top_occ_spec = get_top_occ(df_ura_il, 0.005)
print(len(top_occ_spec))
if not os.path.exists(out_path):
    os.makedirs(out_path)

if not os.path.exists(os.path.join(out_path, 'IL_atlas.csv')):
    curr_df_il = pandas.DataFrame()
    cols = ['age', 'hba1c', 'bmi', 'bt__fasting_glucose', 'bt__fasting_triglycerides', 'bt__hdl_cholesterol']
    for pheno in cols:
        df_curr = throw_outliers(df_meta_il, pheno)
        print("Working on IL %s (%d values after throwing outliers and nans)" % (pheno, len(df_curr)))
        df_il = genes_spear(df_ura_il.loc[df_curr.index][top_occ_spec], df_curr, pheno)
        curr_df_il = curr_df_il.append(df_il[['pheno', 'species', 'spear_pval', 'spear', 'pears_pval', 'pears']],
                                       ignore_index=True)
    curr_df_il['spear_pass_fdr'] = None
    curr_df_il['pears_pass_fdr'] = None
    for pheno in cols:
        inds = curr_df_il.loc[curr_df_il.pheno == pheno].index
        curr_df_il.loc[inds, 'spear_pass_fdr'] = statsmodels.stats.multitest.multipletests(
            curr_df_il.loc[inds, 'spear_pval'].values, method='fdr_bh')[0]
        curr_df_il.loc[inds, 'pears_pass_fdr'] = statsmodels.stats.multitest.multipletests(
            curr_df_il.loc[inds, 'pears_pval'].values, method='fdr_bh')[0]
        print("For %s %d of %d pass FDR" % (pheno, sum(curr_df_il.loc[inds, 'spear_pass_fdr']), len(inds)))

    curr_df_il.to_csv(os.path.join(out_path, 'IL_atlas.csv'))
else:
    curr_df_il = pandas.read_csv(os.path.join(out_path, 'IL_atlas.csv'), index_col=0)

if not os.path.exists(os.path.join(out_path, 'US_atlas.csv')):
    curr_df_us = pandas.DataFrame()
    cols = ['age', 'hba1c', 'bmi']

    for pheno in cols:
        df_curr = throw_outliers(df_meta_us, pheno)
        print("Working on US %s (%d values after throwing outliers and nans)" % (pheno, len(df_curr)))
        df_us = genes_spear(df_ura_us.loc[df_curr.index][top_occ_spec], df_curr, pheno)
        curr_df_us = curr_df_us.append(df_us[['pheno', 'species', 'spear_pval', 'spear', 'pears_pval', 'pears']],
                                       ignore_index=True)

    curr_df_us.columns = ['pheno', 'species', 'spear_pval_us', 'spear_us', 'pears_pval_us', 'pears_us']
    curr_df_us.to_csv(os.path.join(out_path, 'US_atlas.csv'))
else:
    curr_df_us = pandas.read_csv(os.path.join(out_path, 'US_atlas.csv'), index_col=0)

import scipy

print("Figure 5a-c")
curr_df = pandas.merge(curr_df_us, curr_df_il, on=['pheno', 'species'])
cols = ['age', 'hba1c', 'bmi']

for pheno in cols:
    if os.path.exists(os.path.join(out_path, "s_corr_%s_IL_US.png" % pheno)):
        continue
    plt.figure()
    c = numpy.maximum(numpy.log10(curr_df[curr_df.pheno == pheno].spear_pval),
                      numpy.array([-20]*len(curr_df[curr_df.pheno == pheno])))
    plt.scatter(curr_df[curr_df.pheno == pheno].spear, curr_df[curr_df.pheno == pheno].spear_us, c=c)
    cbar = plt.colorbar()
    cbar.set_ticks([cbar.vmax, -20])
    cbar.set_ticklabels(['P=1', 'P<%g' % (10**-20)])
    R = scipy.stats.spearmanr(curr_df[curr_df.pheno == pheno].spear, curr_df[curr_df.pheno == pheno].spear_us,
                              nan_policy='omit')
    plt.title("%s for R=%g p_val=%g" % (pheno, R[0], R[1]))
    plt.xlabel("Train-IL\nSpearmann correlation")
    plt.ylabel("Test-US\nSpearmann correlation")
    plt.savefig(os.path.join(out_path, "s_corr_%s_IL_US.png" % pheno))
plt.close("all")

if not os.path.exists(os.path.join(out_path, 'il_us_spear_saturation.csv')):
    final_df1 = pandas.DataFrame()
    final_df1_0 = pandas.DataFrame()
    final_df1_fdr = pandas.DataFrame()
    cols = ['age', 'hba1c', 'bmi']

    for pheno in cols:
        df_final_il_pheno = throw_outliers(df_meta_il, pheno)
        df_final_us_pheno = throw_outliers(df_meta_us, pheno)
        print("Working on saturation %s (%d, %d values after throwing outliers and nans)" %
              (pheno, len(df_final_il_pheno), len(df_final_us_pheno)))
        max_us = df_final_us_pheno.shape[0]
        for size in [100, 200, 400, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, max_us]:
            if size <= max_us:
                print("Compare %s size %d" % (pheno, size))
                lst = []
                lst0 = []
                lstfdr = []
                for j in range(10):
                    samp_il = df_final_il_pheno.sample(n=size)
                    samp_us = df_final_us_pheno.sample(n=size)
                    df_il = genes_spear(df_ura_il.loc[samp_il.index][top_occ_spec], samp_il, pheno)
                    df_us = genes_spear(df_ura_us.loc[samp_us.index][top_occ_spec], samp_us, pheno)
                    df_il = df_il[['species', 'spear', 'spear_pval']]
                    df_il.columns = ['species', 'spear_il', 'spear_pval_il']
                    df_us = df_us[['species', 'spear', 'spear_pval']]
                    df_us.columns = ['species', 'spear_us', 'spear_pval_us']
                    df_all = pandas.merge(df_us, df_il, on='species')
                    fillna_0 = df_all.fillna(0)
                    drop_na = df_all.dropna()
                    lst.append(stats.spearmanr(list(drop_na.spear_il), list(drop_na.spear_us))[0])
                    lst0.append(stats.spearmanr(list(fillna_0.spear_il), list(fillna_0.spear_us))[0])
                    df_il.dropna(inplace=True)
                    df_il['fdr_il'] = statsmodels.stats.multitest.multipletests(df_il.spear_pval_il.values,
                                                                                method='fdr_bh')[0]
                    df_il = df_il[df_il['fdr_il']]
                    df_us.dropna(inplace=True)
                    df_us['fdr_us'] = statsmodels.stats.multitest.multipletests(df_us.spear_pval_us.values,
                                                                                method='fdr_bh')[0]
                    df_us.dropna(inplace=True)
                    df_us = df_us[df_us['fdr_us']]
                    if (len(df_il) > 0) & (len(df_us) > 0):
                        print("For %s found FDR at size %d (%d, %d in common %d)" %
                              (pheno, size, len(df_il), len(df_us), len(df_il.index.intersection(df_us.index))))
                        df_all = pandas.merge(df_us, df_il, how='outer', on='species')
                        if len(df_all) > 1:
                            df_all.fillna(0, inplace=True)
                            lstfdr.append(stats.spearmanr(list(df_all.spear_il), list(df_all.spear_us))[0])

                final_df1 = final_df1.append({'pheno': pheno, 'size': size, 'r2': numpy.mean(lst),
                                              'std': numpy.std(lst)}, ignore_index=True)
                final_df1_0 = final_df1_0.append({'pheno': pheno, 'size': size, 'r2': numpy.mean(lst0),
                                                'std': numpy.std(lst0)}, ignore_index=True)
                if len(lstfdr) > 5:
                    final_df1_fdr = final_df1_fdr.append({'pheno': pheno, 'size': size, 'r2': numpy.mean(lstfdr),
                                                          'std': numpy.std(lstfdr), 'num_spears': len(lstfdr)},
                                                         ignore_index=True)
    final_df1.to_csv(os.path.join(out_path, 'il_us_spear_saturation.csv'))
    final_df1_0.to_csv(os.path.join(out_path, 'il_us_spear_saturation_0.csv'))
    final_df1_fdr.to_csv(os.path.join(out_path, 'il_us_spear_saturation_FDR.csv'))
else:
    final_df1 = pandas.read_csv(os.path.join(out_path, 'il_us_spear_saturation.csv'), index_col=0)

print("figure 5d-f")
cols = ['age', 'hba1c', 'bmi']

for pheno in cols:
    if os.path.exists(os.path.join(out_path, "cohort_size_%s.png" % pheno)):
        continue
    plt.figure()
    plt.plot(final_df1[final_df1.pheno == pheno]['size'], final_df1[final_df1.pheno == pheno].r2, color='b')
    for s in final_df1[final_df1.pheno == pheno]['size']:
        r = final_df1[(final_df1.pheno == pheno) & (final_df1['size'] == s)].iloc[0].r2
        std = final_df1[(final_df1.pheno == pheno) & (final_df1['size'] == s)].iloc[0]['std']
        plt.vlines(s, r-std, r+std, color='b')
    plt.title("%s" % pheno)
    plt.xlabel("Cohort sizes (IL & US)")
    plt.ylabel("Spearmann correlation\n between cohorts atlas")
    plt.savefig(os.path.join(out_path, "cohort_size_%s.png" % pheno))
plt.close("all")

print("figure 5d-f fillna")
cols = ['age', 'hba1c', 'bmi']

for pheno in cols:
    if os.path.exists(os.path.join(out_path, "cohort_size_%s_fill0.png" % pheno)):
        continue
    plt.figure()
    plt.plot(final_df1_0[final_df1_0.pheno == pheno]['size'], final_df1_0[final_df1_0.pheno == pheno].r2, color='b')
    for s in final_df1_0[final_df1_0.pheno == pheno]['size']:
        r = final_df1_0[(final_df1_0.pheno == pheno) & (final_df1_0['size'] == s)].iloc[0].r2
        std = final_df1_0[(final_df1_0.pheno == pheno) & (final_df1_0['size'] == s)].iloc[0]['std']
        plt.vlines(s, r-std, r+std, color='b')
    plt.title("%s" % pheno)
    plt.xlabel("Cohort sizes (IL & US)")
    plt.ylabel("Spearmann correlation\n between cohorts atlas")
    plt.savefig(os.path.join(out_path, "cohort_size_%s_fill0.png" % pheno))
plt.close("all")

print("figure 5d-f fdr")
cols = ['age', 'hba1c', 'bmi']

for pheno in cols:
    if os.path.exists(os.path.join(out_path, "cohort_size_%s_fdr.png" % pheno)):
        continue
    plt.figure()
    plt.plot(final_df1_fdr[final_df1_fdr.pheno == pheno]['size'], final_df1_fdr[final_df1_fdr.pheno == pheno].r2, color='b')
    for s in final_df1_fdr[final_df1_fdr.pheno == pheno]['size']:
        r = final_df1_fdr[(final_df1_fdr.pheno == pheno) & (final_df1_fdr['size'] == s)].iloc[0].r2
        std = final_df1_fdr[(final_df1_fdr.pheno == pheno) & (final_df1_fdr['size'] == s)].iloc[0]['std']
        plt.vlines(s, r-std, r+std, color='b')
    plt.title("%s" % pheno)
    plt.xlabel("Cohort sizes (IL & US)")
    plt.ylabel("Spearmann correlation\n between cohorts atlas")
    plt.savefig(os.path.join(out_path, "cohort_size_%s_fdr.png" % pheno))
plt.close("all")

if not os.path.exists(os.path.join(out_path, 'FDR')):
    os.makedirs(os.path.join(out_path, 'FDR'))
cols = ['age', 'hba1c', 'bmi']

for pheno in cols:
    fs = glob.glob(os.path.join(out_path, 'FDR', 'single_bac_*_on_%s_variation.csv' % pheno))
    if len(fs) == 0:
        df_final_il_pheno = throw_outliers(df_meta_il, pheno)
        df_final_us_pheno = throw_outliers(df_meta_us, pheno)
        print("Working on FDR variation %s (%d, %d values after throwing outliers and nans)" %
              (pheno, len(df_final_il_pheno), len(df_final_us_pheno)))
        max_us = df_final_us_pheno.shape[0]
        all_bacs = []
        final_df1_fdr = {}
        for size in [100, 200, 400, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200]:
            if size <= max_us:
                print("Compare %s size %d" % (pheno, size))
                lstfdr = {'IL': {}, 'US': {}}
                for j in range(10):
                    samp_il = df_final_il_pheno.sample(n=size)
                    samp_us = df_final_us_pheno.sample(n=size)
                    df_il = genes_spear(df_ura_il.loc[samp_il.index][top_occ_spec], samp_il, pheno)
                    df_us = genes_spear(df_ura_us.loc[samp_us.index][top_occ_spec], samp_us, pheno)
                    df_il = df_il[['species', 'spear', 'spear_pval']]
                    df_il.columns = ['species', 'spear_il', 'spear_pval_il']
                    df_us = df_us[['species', 'spear', 'spear_pval']]
                    df_us.columns = ['species', 'spear_us', 'spear_pval_us']
                    df_il.dropna(inplace=True)
                    df_il['fdr_il'] = statsmodels.stats.multitest.multipletests(df_il.spear_pval_il.values,
                                                                                method='fdr_bh')[0]
                    df_us.dropna(inplace=True)
                    df_us['fdr_us'] = statsmodels.stats.multitest.multipletests(df_us.spear_pval_us.values,
                                                                                method='fdr_bh')[0]
                    df_us.dropna(inplace=True)
                    df_il.set_index('species', inplace=True)
                    df_us.set_index('species', inplace=True)
                    for b in df_il[df_il['fdr_il']].index:
                        all_bacs.append(b)
                        if not b in lstfdr['IL'].keys():
                            lstfdr['IL'][b] = []
                        lstfdr['IL'][b].append(df_il.loc[b, 'spear_il'])
                    for b in df_us[df_us['fdr_us']].index:
                        all_bacs.append(b)
                        if not b in lstfdr['US'].keys():
                            lstfdr['US'][b] = []
                        lstfdr['US'][b].append(df_us.loc[b, 'spear_us'])

                all_bacs = list(set(all_bacs))
                for b in all_bacs:
                    if not b in final_df1_fdr.keys():
                        final_df1_fdr[b] = pandas.DataFrame()
                    d = {'pheno': pheno, 'size': size}
                    for coh in ['IL', 'US']:
                        if b in lstfdr[coh].keys():
                            d['num_' + coh] = len(lstfdr[coh][b])
                            d['min_' + coh] = min(lstfdr[coh][b])
                            d['max_' + coh] = max(lstfdr[coh][b])
                        else:
                            d['num_' + coh] = numpy.nan
                            d['min_' + coh] = numpy.nan
                            d['max_' + coh] = numpy.nan
                    final_df1_fdr[b] = final_df1_fdr[b].append(d, ignore_index=True)
        for b in all_bacs:
            final_df1_fdr[b].to_csv(os.path.join(out_path, 'FDR', 'single_bac_%s_on_%s_variation.csv' %
                                                 (b.split("|")[-1], pheno)))

            print(pheno, b)
            print(final_df1_fdr[b])


for pheno in cols:
    fs = glob.glob(os.path.join(out_path, 'FDR', 'single_bac_*_on_%s_saturation.csv' % pheno))
    if len(fs) == 0:
        final_df1_fdr = None
        df_final_il_pheno = throw_outliers(df_meta_il, pheno)
        df_final_us_pheno = throw_outliers(df_meta_us, pheno)
        print("Working on FDR saturation %s (%d, %d values after throwing outliers and nans)" %
              (pheno, len(df_final_il_pheno), len(df_final_us_pheno)))
        max_us = df_final_us_pheno.shape[0]
        not_chosen = True
        list_bac = None
        for size in [1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, max_us]:
            if size <= max_us:
                print("Compare %s size %d" % (pheno, size))
                lstfdr = {'IL': {}, 'US': {}}
                if list_bac is not None:
                    for b in list_bac:
                        lstfdr['IL'][b] = []
                        lstfdr['US'][b] = []
                    if final_df1_fdr is None:
                        final_df1_fdr = {}
                        for b in list_bac:
                            final_df1_fdr[b] = pandas.DataFrame()
                for j in range(10):
                    samp_il = df_final_il_pheno.sample(n=size)
                    samp_us = df_final_us_pheno.sample(n=size)
                    df_il = genes_spear(df_ura_il.loc[samp_il.index][top_occ_spec], samp_il, pheno)
                    df_us = genes_spear(df_ura_us.loc[samp_us.index][top_occ_spec], samp_us, pheno)
                    df_il = df_il[['species', 'spear', 'spear_pval']]
                    df_il.columns = ['species', 'spear_il', 'spear_pval_il']
                    df_us = df_us[['species', 'spear', 'spear_pval']]
                    df_us.columns = ['species', 'spear_us', 'spear_pval_us']
                    df_il.dropna(inplace=True)
                    df_il['fdr_il'] = statsmodels.stats.multitest.multipletests(df_il.spear_pval_il.values,
                                                                                method='fdr_bh')[0]
                    df_us.dropna(inplace=True)
                    df_us['fdr_us'] = statsmodels.stats.multitest.multipletests(df_us.spear_pval_us.values,
                                                                                method='fdr_bh')[0]
                    df_us.dropna(inplace=True)
                    if not_chosen:
                        if list_bac is None:
                            list_bac = set(df_il[df_il['fdr_il']].species.values)
                        else:
                            list_bac = list_bac.intersection(df_il[df_il['fdr_il']].species.values)
                        list_bac = list_bac.intersection(df_us[df_us['fdr_us']].species.values)
                    else:
                        df_il.set_index('species', inplace=True)
                        df_us.set_index('species', inplace=True)
                        for b in list_bac:
                            if (b in df_il.index) and (not numpy.isnan(df_il.loc[b].spear_il)):
                                lstfdr['IL'][b].append(df_il.loc[b].spear_il)
                            else:
                                lstfdr['IL'][b].append(0)
                            if (b in df_us.index) and (not numpy.isnan(df_us.loc[b].spear_us)):
                                lstfdr['US'][b].append(df_us.loc[b].spear_us)
                            else:
                                lstfdr['US'][b].append(0)

                if not not_chosen:
                    for b in list_bac:
                        final_df1_fdr[b] = final_df1_fdr[b].append({'pheno': pheno, 'size': size,
                                                                    'mean_il': numpy.mean(lstfdr['IL'][b]),
                                                                    'std_il': numpy.std(lstfdr['IL'][b]),
                                                                    'mean_us': numpy.mean(lstfdr['US'][b]),
                                                                    'std_us': numpy.std(lstfdr['US'][b])},
                                                                   ignore_index=True)
            if not_chosen:
                if len(list_bac) > 0:
                    list_bac = list(list_bac)
                    print("Using size %d chose %d bacteria" % (size, len(list_bac)))
                    print(list_bac)
                    not_chosen = False
                else:
                    print("Couldn't choose with size %d. Trying next size" % size)
                    list_bac = None
        for b in list_bac:
            final_df1_fdr[b].to_csv(os.path.join(out_path, 'FDR', 'single_bac_%s_on_%s_saturation.csv' %
                                                 (b.split("|")[-1], pheno)))

            print(pheno, b)
            print(final_df1_fdr[b])
