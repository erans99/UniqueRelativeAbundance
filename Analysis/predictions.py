import pandas
import numpy
import os
import sys
import matplotlib.pyplot as plt
from scipy import stats
import socket


ONLY_FULL = False #True
NUM_TH = 1

data_path = "/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/DataForRevision"
out_path = os.path.join(data_path, "Predictions")
atlas_path = os.path.join(data_path, "Atlas")
if ONLY_FULL:
    out_path += "_FullLen"
    atlas_path += "_FullLen"

if not os.path.exists(out_path):
    os.makedirs(out_path)

if not os.path.exists(os.path.join(out_path, "all_preds")):
    os.makedirs(os.path.join(out_path, "all_preds"))


MIN_FILL = 0.0001
params = {'objective': 'reg:squarederror', 'colsample_bylevel': 0.075, 'max_depth': 6, 'learning_rate': 0.0025,
          'n_estimators': 4000, 'subsample': 0.6, 'min_child_weight': 20}
sub_params = {'objective': 'reg:squarederror', 'colsample_bylevel': 0.075, 'max_depth': 6, 'learning_rate': 0.01,
              'n_estimators': 1700, 'subsample': 0.6, 'min_child_weight': 80}


def throw_outliers(df_meta, pheno, return_df=False, num_stds=4):
    if num_stds <= 0:
        print("WTF. Need to throw above positive number of STDs from mean")
    curr_005 = df_meta[pheno].quantile(0.05)
    curr_095 = df_meta[pheno].quantile(0.95)
    df_curr = df_meta[(df_meta[pheno] >= curr_005) & (df_meta[pheno] <= curr_095)]
    curr_median = df_curr[pheno].median()
    curr_stdev = df_curr[pheno].std()
    print("For %s: med %g std %g. Thus throwing outsude %g %g." % (pheno, curr_median, curr_stdev,
                                                                   (curr_median - num_stds * curr_stdev),
                                                                   (curr_median + num_stds * curr_stdev)))
    # take all data <= num_stds stds from median
    if return_df:
        return df_meta[(df_meta[pheno] >= (curr_median - num_stds * curr_stdev)) &
                       (df_meta[pheno] <= (curr_median + num_stds * curr_stdev))]
    else:
        return df_meta[(df_meta[pheno] >= (curr_median - num_stds * curr_stdev)) &
                       (df_meta[pheno] <= (curr_median + num_stds * curr_stdev))][pheno]


def get_top_occ(df_in, oc):
    df = df_in.copy()
    lst_species = []
    df[df == MIN_FILL] = numpy.nan
    for c in df_in.columns:
        if df[c].count() > oc * df.shape[0]:
            lst_species.append(c)

    return lst_species

from sklearn.metrics import r2_score, roc_curve, auc
import xgboost
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import SGDClassifier

predictor = xgboost.XGBRegressor(nthread=NUM_TH, **params)
predictor_class = xgboost.XGBClassifier(nthread=NUM_TH, **params)
predictor_sub = xgboost.XGBRegressor(nthread=NUM_TH, **sub_params)

pred_ridge = RidgeCV(alphas=[0.1, 1, 10, 100, 1000], normalize=True)
pred_sgdc = SGDClassifier(loss='log')


def fit_and_validate(predictor, in_data, validations):
    df_in, df_out = in_data
    predictor.fit(df_in, list(df_out))

    res = []
    for validation in validations:
        df_vin, df_vout = validation
        pred = predictor.predict(df_vin)
        res.append(r2_score(list(df_vout), pred))
    return res


def pipline_cross_val_x(predictor, in_data, validation_list, cv, validate, is_classifier, is_log, save_pred=None,
                        save_name=None):
    df_in, df_out = in_data

    saved = [[], []]
    lst_clients = list(df_in.index)
    part = int(len(lst_clients) / cv)
    print("Part len %d" % part)
    lst_res = []
    for i in range(cv):
        start = i * part
        end = start + part
        test_clients = lst_clients[start:end]
        train_clients = df_in.index.difference(test_clients)
        predictor.fit(df_in.loc[train_clients], list(df_out.loc[train_clients]))
        if is_classifier:
            prob = predictor.predict_proba(df_in.loc[test_clients])
            saved[0] += list(df_out.loc[test_clients])
            saved[1] += list(prob[:, 0])
            fpr, tpr, threshold = roc_curve(list(df_out.loc[test_clients]), 1 - numpy.array(list(prob[:, 0])))
            lst_res.append(auc(fpr, tpr))
        else:
            pred = predictor.predict(df_in.loc[test_clients])
            saved[0] += list(df_out.loc[test_clients])
            saved[1] += list(pred)
            lst_res.append(r2_score(list(df_out.loc[test_clients]), pred))

    if save_pred is not None:
        if is_classifier:
            pandas.DataFrame(saved, index=["real_val", 'pred_prob']).T.to_csv(save_pred)
        else:
            pandas.DataFrame(saved, index=["real_val", 'pred_val']).T.to_csv(save_pred)
    lst_means = []
    lst_stds = []
    if validate:
        predictor.fit(df_in, list(df_out))
        if save_name is not None:
            print("Saving %s" % os.path.basename(save_name))
            predictor.save_model(save_name)
        for validation in validation_list:
            df_vin, df_vout = validation
            lst_val = []
            lst_clients = list(df_vin.index)
            part = int(len(lst_clients) / cv)
            print("Validation part len %d" % part)

            for i in range(cv):
                start = i * part
                end = start + part
                val_clients = lst_clients[start:end]
                pred = predictor.predict(df_vin.loc[val_clients])
                if len(pred) == 0:
                    print("WTF")
                    continue
                lst_val.append(r2_score(list(df_vout.loc[val_clients]), pred))
            lst_means.append(numpy.mean(lst_val))
            lst_stds.append(numpy.std(lst_val))
    return (numpy.mean(lst_res), numpy.std(lst_res)), (lst_means, lst_stds)


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
    for name, df, df_meta in [["IL", df_ura_il, df_meta_il], ["US", df_ura_us, df_meta_us],
                              ["IL_val", df_ura_il_val, df_meta_il_val]]:
        good_inds = df.index.intersection(link[((link.RawReadLength == 75) & (link.URAMapUsed == 8*10**6)) |
                                               ((link.RawReadLength == 100) & (link.URAMapUsed == 5*10**6))].index)
        print("For %s left with %d of %d" % (name, len(good_inds), len(df)))
        df = df.loc[good_inds]
        df_meta = df_meta.loc[good_inds]

top_occ_spec = get_top_occ(df_ura_il, 0.005)
print(len(top_occ_spec))
if not os.path.exists(out_path):
    os.makedirs(out_path)

pandas.Series(top_occ_spec).to_csv(os.path.join(out_path, "mdl_SGBs_list.csv"))

df_meta_il = df_meta_il[df_meta_il.gender.notnull()]
df_meta_il = df_meta_il[df_meta_il.age.notnull()]

df_meta_us = df_meta_us[df_meta_us.gender.notnull()]
df_meta_us = df_meta_us[df_meta_us.age.notnull()]

df_meta_il_val = df_meta_il_val[df_meta_il_val.gender.notnull()]
df_meta_il_val = df_meta_il_val[df_meta_il_val.age.notnull()]


# figure 4a
df = pandas.DataFrame()
if not os.path.exists(os.path.join(out_path, 'prediction_xgboostMB_gender.csv')):
    cols = ['age', 'hba1c', 'bt__fasting_glucose', 'bmi', 'bt__fasting_triglycerides', 'height',
            'bowel_movement_frequency', 'bt__hdl_cholesterol', 'bt__sgpt', 'bt__cholesterol',
            'bt__protein', 'bt__sgot', 'bt__albumin', 'bt__alkaline_phosphatase', 'bt__tsh', 'bt__inr']

    for pheno in cols:
        df_curr_il = throw_outliers(df_meta_il, pheno, True)
        print("Figure 4a %s (%d values after throwing outliers and nans)" % (pheno, len(df_curr_il)))
        df_curr_male = df_curr_il[df_curr_il.gender == 1]
        df_curr_female = df_curr_il[df_curr_il.gender == 0]
        pheno_segata, x = pipline_cross_val_x(predictor,
                                              [df_ura_il.loc[df_curr_il.index][top_occ_spec], df_curr_il[pheno]], [],
                                              10, False, False, False,
                                              save_pred=os.path.join(out_path, "all_preds",
                                                                     "xgb_MB_pred_%s.csv" % pheno))
        pheno_segata_female, x = pipline_cross_val_x(predictor, [df_ura_il.loc[df_curr_female.index][top_occ_spec],
                                                            df_curr_female[pheno]], [], 10, False, False, False)
        pheno_segata_male, x = pipline_cross_val_x(predictor, [df_ura_il.loc[df_curr_male.index][top_occ_spec],
                                                          df_curr_male[pheno]], [], 10, False, False, False)

        df = df.append({'gender': 'all', 'target': pheno, 'mean_pearson': pheno_segata[0],
                        'std_pearson': pheno_segata[1], 'cohort_size': df_curr_il.shape[0]}, ignore_index=True)
        df = df.append({'gender': 'female', 'target': pheno, 'mean_pearson': pheno_segata_female[0],
                        'std_pearson': pheno_segata_female[1], 'cohort_size': df_curr_female.shape[0]},
                       ignore_index=True)
        df = df.append({'gender': 'male', 'target': pheno, 'mean_pearson': pheno_segata_male[0],
                        'std_pearson': pheno_segata_male[1], 'cohort_size': df_curr_male.shape[0]},
                       ignore_index=True)
        df.to_csv(os.path.join(out_path, 'prediction_xgboostMB_gender.csv'))
    print("Done 4a")

# figure 4b
if not os.path.exists(os.path.join(out_path, 'classification_xgboostMB_phenotypes.csv')):
    df = pandas.DataFrame()
    cols = ['type_2_diabetes', 'gender', 'ever_smoked', 'currently_smokes']
    for pheno in cols:
        print("Classification %s" % pheno)
        df_pheno = df_meta_il[df_meta_il[pheno].notnull()]
        num_removed_nulls = df_meta_il.shape[0] - df_pheno.shape[0]
        df_ura = numpy.log(df_ura_il.loc[df_pheno.index][top_occ_spec])

        df_curr_male = df_pheno[df_pheno.gender == 1]
        df_curr_female = df_pheno[df_pheno.gender == 0]
        pheno_segata, x = pipline_cross_val_x(predictor_class, [df_ura.loc[df_pheno.index][top_occ_spec], df_pheno[pheno]],
                                              [], 10, False, True, False,
                                              save_pred=os.path.join(out_path, "all_preds",
                                                                     "xgb_MB_pred_%s.csv" % pheno))
        if pheno != 'gender':
            pheno_segata_female, x = pipline_cross_val_x(predictor_class, [df_ura.loc[df_curr_female.index][top_occ_spec],
                                                                df_curr_female[pheno]], [], 10, False, True, False)
            pheno_segata_male, x = pipline_cross_val_x(predictor_class, [df_ura.loc[df_curr_male.index][top_occ_spec],
                                                              df_curr_male[pheno]], [], 10, False, True, False)
        df = df.append({'gender': 'all', 'target': pheno, 'auc': pheno_segata[0], 'stdev': pheno_segata[1],
                        'n': df_pheno.shape[0]}, ignore_index=True)
        if pheno != 'gender':
            df = df.append({'gender': 'female', 'target': pheno, 'auc': pheno_segata_female[0],
                            'stdev': pheno_segata_female[1],
                            'n': df_curr_female.shape[0]}, ignore_index=True)
            df = df.append({'gender': 'male', 'target': pheno, 'auc': pheno_segata_male[0],
                            'stdev': pheno_segata_male[1],
                            'n': df_curr_male.shape[0]}, ignore_index=True)
    df.to_csv(os.path.join(out_path, 'classification_xgboostMB_phenotypes.csv'))
    print("Done 4b")

# figure 4c-f
cols = ['hba1c', 'bmi', 'age']
for pheno in cols:
    if not os.path.exists(os.path.join(out_path, 'prediction_xgboostMB_%s_ext_and_validate.csv' % pheno)):
        extra_cols = ['age', 'gender']
        extracols_name = 'age+sex'
        if pheno == 'age':
            extra_cols = ['gender']
            extracols_name = 'sex'
        df_curr_il = throw_outliers(df_meta_il, pheno, True)
        df_curr_us = throw_outliers(df_meta_us, pheno, True)
        df_curr_il_val = throw_outliers(df_meta_il_val, pheno, True)
        print("Figure 4c-f %s (%d %d %d) values after throwing outliers and nans)" %
              (pheno, len(df_curr_il), len(df_curr_us), len(df_curr_il_val)))

        pheno_segata, x = pipline_cross_val_x(predictor,
                                              [df_ura_il.loc[df_curr_il.index][top_occ_spec], df_curr_il[pheno]],
                                              [[df_ura_us.loc[df_curr_us.index][top_occ_spec], df_curr_us[pheno]],
                                               [df_ura_il_val.loc[df_curr_il_val.index][top_occ_spec],
                                                df_curr_il_val[pheno]]],
                                              10, True, False, False,
                                              save_pred=os.path.join(out_path, "xgb_MB_pred_%s.csv" % pheno),
                                              save_name=os.path.join(out_path, "xgb_MB_%s.mdl" % pheno))
        pheno_age_sex, x_only_age_sex = pipline_cross_val_x(predictor, [df_curr_il[extra_cols], df_curr_il[pheno]],
                                                      [[df_curr_us[extra_cols], df_curr_us[pheno]],
                                                       [df_curr_il_val[extra_cols], df_curr_il_val[pheno]]],
                                                      10, True, False, False)
        tmp = [df_ura_il.loc[df_curr_il.index][top_occ_spec].merge(df_curr_il[extra_cols],
                                                                   left_index=True, right_index=True),
               df_ura_us.loc[df_curr_us.index][top_occ_spec].merge(df_curr_us[extra_cols],
                                                                   left_index=True, right_index=True),
               df_ura_il_val.loc[df_curr_il_val.index][top_occ_spec].merge(df_curr_il_val[extra_cols],
                                                                   left_index=True, right_index=True),
               ]
        pheno_segata_age_sex, x_age_sex = pipline_cross_val_x(predictor, [tmp[0], df_curr_il[pheno]],
                                              [[tmp[1], df_curr_us[pheno]], [tmp[2], df_curr_il_val[pheno]]],
                                                              10, True, False, False,
                                                        save_name=os.path.join(out_path, "xgb_MBext_%s.mdl" % pheno))

        df_xgb_res = pandas.DataFrame()

        dct = {}
        dct['cohort'] = 'il'
        dct['features'] = 'sgbs+' + extracols_name
        dct['is_validation'] = 0
        dct['mean_pearson'] = pheno_segata_age_sex[0]
        dct['std_pearson'] = pheno_segata_age_sex[1]
        dct['cohort_size'] = df_curr_il.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_il'
        dct['features'] = 'sgbs+' + extracols_name
        dct['is_validation'] = 1
        dct['mean_pearson'] = x_age_sex[0][1]
        dct['std_pearson'] = x_age_sex[1][1]
        dct['cohort_size'] = df_curr_il_val.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_us'
        dct['features'] = 'sgbs+' + extracols_name
        dct['is_validation'] = 1
        dct['mean_pearson'] = x_age_sex[0][0]
        dct['std_pearson'] = x_age_sex[1][0]
        dct['cohort_size'] = df_curr_us.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'il'
        dct['features'] = 'sgbs'
        dct['is_validation'] = 0
        dct['mean_pearson'] = pheno_segata[0]
        dct['std_pearson'] = pheno_segata[1]
        dct['cohort_size'] = df_curr_il.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_il'
        dct['features'] = 'sgbs'
        dct['is_validation'] = 1
        dct['mean_pearson'] = x[0][1]
        dct['std_pearson'] = x[1][1]
        dct['cohort_size'] = df_curr_il_val.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_us'
        dct['features'] = 'sgbs'
        dct['is_validation'] = 1
        dct['mean_pearson'] = x[0][0]
        dct['std_pearson'] = x[1][0]
        dct['cohort_size'] = df_curr_us.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'il'
        dct['features'] = extracols_name
        dct['is_validation'] = 0
        dct['mean_pearson'] = pheno_age_sex[0]
        dct['std_pearson'] = pheno_age_sex[1]
        dct['cohort_size'] = df_curr_il.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_il'
        dct['features'] = extracols_name
        dct['is_validation'] = 1
        dct['mean_pearson'] = x_only_age_sex[0][1]
        dct['std_pearson'] = x_only_age_sex[1][1]
        dct['cohort_size'] = df_curr_il_val.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        dct = {}
        dct['cohort'] = 'validation_us'
        dct['features'] = extracols_name
        dct['is_validation'] = 1
        dct['mean_pearson'] = x_only_age_sex[0][0]
        dct['std_pearson'] = x_only_age_sex[1][0]
        dct['cohort_size'] = df_curr_us.shape[0]
        df_xgb_res = df_xgb_res.append(dct, ignore_index=True)

        df_xgb_res.to_csv(os.path.join(out_path, 'prediction_xgboostMB_%s_ext_and_validate.csv' % pheno))
print("Done 4c-f")

# figure 4g-i - old CV
cols = ['hba1c', 'bmi', 'age']
for pheno in cols:
    print("Saturation XGB %s" % pheno)
    if not os.path.exists(os.path.join(out_path, 'saturation_xgboostMB_%s.csv' % pheno)):
        df = pandas.DataFrame()
        df_pheno = throw_outliers(df_meta_il, pheno)

        for size in [200, 500, 1000, 2000, 3000, 4000, 6000, 8000, 12000, 16000, 20000, 24000]:
            if size > (0.9*len(df_pheno)):
                print("for pheno %s, %d is too high (of %d)" % (pheno, size, len(df_pheno)))
                continue
            print(pheno, size)
            lst = []
            lst_std = []
            for i in range(10):
                df_sub = df_pheno.sample(n=size)
                pheno_segata, x = pipline_cross_val_x(predictor_sub, [df_ura_il.loc[df_sub.index][top_occ_spec], df_sub],
                                                      [], 10, False, False, False)

                lst.append(pheno_segata[0])
                lst_std.append(pheno_segata[1])
            df = df.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst), 'std_pearson': numpy.std(lst),
                            'mean_std': numpy.mean(lst_std)}, ignore_index=True)
        df.to_csv(os.path.join(out_path, 'saturation_xgboostMB_%s.csv' % pheno))

    if not os.path.exists(os.path.join(out_path, 'saturation_RidgeMB_%s.csv' % pheno)):
        df = pandas.DataFrame()
        df_pheno = throw_outliers(df_meta_il, pheno)
        df_ura = numpy.log(df_ura_il.loc[df_pheno.index][top_occ_spec])

        for size in [200, 500, 1000, 2000, 3000, 4000, 6000, 8000, 12000, 16000, 20000, 24000]:
            if size > (0.9*len(df_pheno)):
                print("for pheno %s, %d is too high (of %d)" % (pheno, size, len(df_pheno)))
                continue
            print(pheno, size, "linear")
            lst = []
            lst_std = []
            for i in range(10):
                df_sub = df_pheno.sample(n=size)
                pheno_segata, x = pipline_cross_val_x(pred_ridge, [df_ura.loc[df_sub.index][top_occ_spec], df_sub],
                                                      [], 10, False, False, False)
                lst.append(pheno_segata[0])
                lst_std.append(pheno_segata[1])
            df = df.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst), 'std_pearson': numpy.std(lst),
                            'mean_std': numpy.mean(lst_std)}, ignore_index=True)
        df.to_csv(os.path.join(out_path, 'saturation_RidgeMB_%s.csv' % pheno))
print("Done 4g-i")

# figure 4g-i & S4a-c - IL_val, US
print("4g-i & S4a-c - IL_val/US")
cols = ['age', 'hba1c', 'bmi']
for pheno in cols:
    if not os.path.exists(os.path.join(out_path, 'saturation_xgboostMB_%s_vs_ILval.csv' % pheno)):
        dfXGB_US = pandas.DataFrame()
        dfXGB_IL = pandas.DataFrame()
        df_pheno = throw_outliers(df_meta_il, pheno)
        us_pheno = throw_outliers(df_meta_us, pheno)
        il_pheno = throw_outliers(df_meta_il_val, pheno)
        print("Saturation XGB %s (%d %d %d) values after throwing outliers and nans)" %
              (pheno, len(df_pheno), len(us_pheno), len(il_pheno)))
        sizes = [200, 500, 1000, 2000, 3000, 4000, 6000, 8000, 12000, 16000, 20000, 24000]

        for size in sizes:
            if size > (0.9*len(df_pheno)):
                print("for pheno %s, %d is too high (of %d)" % (pheno, size, len(df_pheno)))
                continue
            print(pheno, size)
            lst = [[], []]
            for i in range(10):
                df_sub = df_pheno.sample(n=size)
                R2s = fit_and_validate(predictor_sub, [df_ura_il.loc[df_sub.index][top_occ_spec], df_sub],
                                       [[df_ura_us.loc[us_pheno.index][top_occ_spec], us_pheno],
                                        [df_ura_il_val.loc[il_pheno.index][top_occ_spec], il_pheno]])

                lst[0].append(R2s[0])
                lst[1].append(R2s[1])
            dfXGB_US = dfXGB_US.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst[0]),
                                        'std_pearson': numpy.std(lst[0])}, ignore_index=True)
            dfXGB_IL = dfXGB_IL.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst[1]),
                                        'std_pearson': numpy.std(lst[1])}, ignore_index=True)
        dfXGB_US.to_csv(os.path.join(out_path, 'saturation_xgboostMB_%s_vs_US.csv' % pheno))
        dfXGB_IL.to_csv(os.path.join(out_path, 'saturation_xgboostMB_%s_vs_ILval.csv' % pheno))
    else:
        dfXGB_US = pandas.read_csv(os.path.join(out_path, 'saturation_xgboostMB_%s_vs_US.csv' % pheno), index_col=0)
        dfXGB_IL = pandas.read_csv(os.path.join(out_path, 'saturation_xgboostMB_%s_vs_ILval.csv' % pheno), index_col=0)

    if not os.path.exists(os.path.join(out_path, 'saturation_RidgeMB_%s_vs_ILval.csv' % pheno)):
        dfRig_US = pandas.DataFrame()
        dfRig_IL = pandas.DataFrame()
        df_pheno = throw_outliers(df_meta_il, pheno)
        us_pheno = throw_outliers(df_meta_us, pheno)
        il_pheno = throw_outliers(df_meta_il_val, pheno)
        print("Saturation Ridge %s (%d %d %d) values after throwing outliers and nans)" %
              (pheno, len(df_pheno), len(us_pheno), len(il_pheno)))
        df_ura = numpy.log(df_ura_il.loc[df_pheno.index][top_occ_spec])
        us_ura = numpy.log(df_ura_us.loc[us_pheno.index][top_occ_spec])
        il_ura = numpy.log(df_ura_il_val.loc[il_pheno.index][top_occ_spec])

        for size in sizes:
            if size > (0.9*len(df_pheno)):
                print("for pheno %s, %d is too high (of %d)" % (pheno, size, len(df_pheno)))
                continue
            print(pheno, size, "linear")
            lst = [[], []]
            for i in range(10):
                df_sub = df_pheno.sample(n=size)
                R2s = fit_and_validate(pred_ridge, [df_ura.loc[df_sub.index], df_sub],
                                       [[us_ura, us_pheno], [il_ura, il_pheno]])

                lst[0].append(R2s[0])
                lst[1].append(R2s[1])
            dfRig_US = dfRig_US.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst[0]),
                                        'std_pearson': numpy.std(lst[0])}, ignore_index=True)
            dfRig_IL = dfRig_IL.append({'cohort_size': size, 'mean_pearson': numpy.mean(lst[1]),
                                        'std_pearson': numpy.std(lst[1])}, ignore_index=True)
        dfRig_US.to_csv(os.path.join(out_path, 'saturation_RidgeMB_%s_vs_US.csv' % pheno))
        dfRig_IL.to_csv(os.path.join(out_path, 'saturation_RidgeMB_%s_vs_ILval.csv' % pheno))
    else:
        dfRig_US = pandas.read_csv(os.path.join(out_path, 'saturation_RidgeMB_%s_vs_US.csv' % pheno), index_col=0)
        dfRig_IL = pandas.read_csv(os.path.join(out_path, 'saturation_RidgeMB_%s_vs_ILval.csv' % pheno), index_col=0)
    plt.figure()
    plt.plot(dfXGB_US.cohort_size.values, dfXGB_US.mean_pearson.values, color="b", label="GBDT US")
    plt.plot(dfXGB_IL.cohort_size.values, dfXGB_IL.mean_pearson.values, color="k", label="GBDT IL")
    plt.xlabel("Sample Size")
    plt.ylabel("%s R2 (%%)" % pheno)
    ylim = plt.ylim()
    plt.title("Predict %s from small cohort, on validation data" % pheno)
    # plt.legend()
    # plt.savefig(os.path.join(out_path, 'saturation_XGB_%s_vs_US.png' % pheno))
    plt.plot(dfRig_US.cohort_size.values, dfRig_US.mean_pearson.values, color="g", label="Ridge US")
    plt.plot(dfRig_IL.cohort_size.values, dfRig_IL.mean_pearson.values, color="y", label="Ridge IL")
    plt.ylim(ylim)
    plt.legend()
    plt.savefig(os.path.join(out_path, 'saturation_%s_vs_US.png' % pheno))
    plt.close("all")

print("Done new 4g-i")

# supplamentary figure S1a
df = pandas.DataFrame()
if not os.path.exists(os.path.join(out_path, 'prediction_Ridge_linear.csv')):
    cols = ['age', 'hba1c', 'bt__fasting_glucose', 'bmi', 'bt__fasting_triglycerides', 'height',
            'bowel_movement_frequency', 'bt__hdl_cholesterol', 'bt__sgpt', 'bt__cholesterol',
            'bt__protein', 'bt__sgot', 'bt__albumin', 'bt__alkaline_phosphatase', 'bt__tsh', 'bt__inr']

    for pheno in cols:
        df_curr_il = throw_outliers(df_meta_il, pheno, True)
        df_curr_us = throw_outliers(df_meta_us, pheno, True)
        df_curr_il_val = throw_outliers(df_meta_il_val, pheno, True)
        print("Figure S1a %s (%d %d %d) values after throwing outliers and nans)" %
              (pheno, len(df_curr_il), len(df_curr_us), len(df_curr_il_val)))

        pheno_segata, x = pipline_cross_val_x(pred_ridge,
                                              [df_ura_il.loc[df_curr_il.index][top_occ_spec], df_curr_il[pheno]],
                                              [], 10, False, False, False,
                                              save_pred=os.path.join(out_path, "all_preds",
                                                                     "Ridge_MB_pred_%s.csv" % pheno))
        df = df.append({'target': pheno, 'mean_pearson': pheno_segata[0],
                        'std_pearson': pheno_segata[1], 'cohort_size': df_curr_il.shape[0]}, ignore_index=True)
        df.to_csv(os.path.join(out_path, 'prediction_Ridge_linear.csv'))

# supplamentary figure S1b
if not os.path.exists(os.path.join(out_path, 'classification_sgdMB_phenotypes.csv')):
    df = pandas.DataFrame()
    cols = ['type_2_diabetes', 'gender', 'ever_smoked', 'currently_smokes']
    for pheno in cols:
        print("Figure S1b %s" % pheno)
        df_pheno = df_meta_il[df_meta_il[pheno].notnull()]
        num_removed_nulls = df_meta_il.shape[0] - df_pheno.shape[0]
        df_ura = numpy.log(df_ura_il.loc[df_pheno.index][top_occ_spec])
        # print(pheno, num_removed_nulls)

        df_curr_male = df_pheno[df_pheno.gender == 1]
        df_curr_female = df_pheno[df_pheno.gender == 0]
        pheno_segata, x = pipline_cross_val_x(pred_sgdc, [df_ura.loc[df_pheno.index][top_occ_spec], df_pheno[pheno]],
                                              [], 10, False, True, False,
                                              save_pred=os.path.join(out_path, "all_preds",
                                                                     "Ridge_MB_pred_%s.csv" % pheno))
        df = df.append({'target': pheno, 'auc': pheno_segata[0], 'stdev': pheno_segata[1],
                        'n': df_pheno.shape[0]}, ignore_index=True)
    df.to_csv(os.path.join(out_path, 'classification_sgdMB_phenotypes.csv'))
    print("Done Sup1b")

print("For Atlas")
if not os.path.exists(os.path.join(out_path, 'linear_coeffs_new.csv')):
    linear = pandas.DataFrame()
    df_ura_il_curr = numpy.log(df_ura_il[top_occ_spec])
    cols = ['age', 'hba1c', 'bmi', 'bt__fasting_glucose', 'bt__fasting_triglycerides', 'bt__hdl_cholesterol']
    for pheno in cols:
        df_curr = throw_outliers(df_meta_il, pheno)

        pred_ridge.fit(df_ura_il_curr.loc[df_curr.index], df_curr)
        df_temp = pandas.DataFrame()
        df_temp['species'] = top_occ_spec
        df_temp['lr_coeff'] = pred_ridge.coef_
        df_temp['intercept'] = pred_ridge.intercept_
        df_temp['pheno'] = pheno
        linear = linear.append(df_temp, ignore_index=True)
    linear.to_csv(os.path.join(out_path, 'linear_coeffs_new.csv'))
else:
    linear = pandas.read_csv(os.path.join(out_path, 'linear_coeffs_new.csv'), index_col=0)

base_atlas = pandas.read_csv(os.path.join(atlas_path, 'IL_atlas.csv'), index_col=0)
atlas = {}
intercept = {}
cols = ['age', 'hba1c', 'bmi', 'bt__fasting_glucose', 'bt__fasting_triglycerides', 'bt__hdl_cholesterol']
for pheno in cols:
    atlas[pheno] = pandas.merge(base_atlas[base_atlas.pheno == pheno], linear[linear.pheno == pheno], on='species')
    intercept[pheno] = atlas[pheno].iloc[0]["intercept"]
    atlas[pheno].drop(["pheno_x", "pheno_y", "intercept"], axis=1, inplace=True)


import shap

df_all_shaps = pandas.DataFrame()
cols = ['age', 'hba1c', 'bmi', 'bt__fasting_glucose', 'bt__fasting_triglycerides', 'bt__hdl_cholesterol']
for pheno in cols:
    if not os.path.exists(os.path.join(out_path, 'top_phenotypes_xgb_shaps_%s.pkl' % pheno)):
        print("Shap for %s" % pheno)
        df_curr = throw_outliers(df_meta_il, pheno)

        shap.initjs()
        predictor.fit(df_ura_il.loc[df_curr.index][top_occ_spec], list(df_curr))
        predictor.save_model(os.path.join(out_path, "pred_%s.mdl" % pheno))
        print("fit")
        shap_sub = df_ura_il.loc[df_curr.index][top_occ_spec].sample(n=1000)
        shap_sub.columns = [x.split("|")[-1] for x in top_occ_spec]
        explainer = shap.TreeExplainer(predictor)
        print("explain")
        shap_values = explainer.shap_values(shap_sub)
        print("get")
        shap.summary_plot(shap_values, shap_sub, max_display=20, show=False)
        plt.title("Top shap for %s" % pheno)
        plt.tight_layout()
        plt.savefig(os.path.join(out_path, 'shap_%s.png' % pheno))
        plt.close()

        lst_features = []
        lst_shap_val = []
        lst_spear = []
        for i in range(len(top_occ_spec)):
            lst_features.append(top_occ_spec[i])
            lst_spear.append(stats.spearmanr(shap_values[:, i], shap_sub[top_occ_spec[i].split("|")[-1]])[0])
            lst_shap_val.append((numpy.sum(numpy.abs(shap_values[:, i]), axis=0)) / float(len(shap_values[:, i])))
        df_shap = pandas.DataFrame()
        df_shap['species'] = lst_features
        df_shap['shap'] = lst_shap_val
        df_shap['shap_spearman'] = lst_spear
        df_shap['pheno'] = pheno
        df_shap.to_pickle(os.path.join(out_path, 'top_phenotypes_xgb_shaps_%s.pkl' % pheno))
    else:
        df_shap = pandas.read_pickle(os.path.join(out_path, 'top_phenotypes_xgb_shaps_%s.pkl' % pheno))
    df_all_shaps = df_all_shaps.append(df_shap, ignore_index=True)

new_cols = ['Spearman correlation P value', 'Spearman correlation coefficient',
            'Pearson correlation P value', 'Pearson correlation coefficient',
            'Mean absolute shap value', 'Ridge regression coefficient', 'Species']

cols = ['age', 'hba1c', 'bmi', 'bt__fasting_glucose', 'bt__fasting_triglycerides', 'bt__hdl_cholesterol']
for pheno in cols:
    if not os.path.exists(os.path.join(atlas_path, 'final_atlas_%s.csv' % pheno)):
        atlas[pheno] = pandas.merge(atlas[pheno], df_all_shaps[df_all_shaps.pheno == pheno], on='species')
        atlas[pheno].drop(["pheno"], axis=1, inplace=True)
        atlas[pheno].sort_values("spear_pval", inplace=True)
        atlas[pheno] = atlas[pheno][['spear_pval', 'spear', 'pears_pval', 'pears', 'shap', 'lr_coeff', 'species']]
        atlas[pheno].sort_values("spear_pval", inplace=True)
        atlas[pheno].columns = new_cols
        atlas[pheno].to_csv(os.path.join(atlas_path, 'final_atlas_%s.csv' % pheno))
        open(os.path.join(atlas_path, 'intercept_%s.txt' % pheno), "w").write(str(intercept[pheno]))
