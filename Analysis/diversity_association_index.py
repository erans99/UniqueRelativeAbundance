from numpy.random import choice
import os
import pandas as pd
from statsmodels.formula.api import ols,logit
import numpy as np
from lib.queue.qp import qp, fakeqp
from lib.addloglevels import sethandlers
import glob
from math import pi
basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/Alpha-Diversity'
jobs_output = os.path.join(basepath,'jobs')

keephenotypes = {'age':'Age',
          'bmi':'BMI',
          'hba1c':'HbA1C%',
          'bt__fasting_glucose':'Glucose',
          'bt__fasting_triglycerides':'Triglycerides',
          'bt__hdl_cholesterol': 'HDL cholesterol',
          'height':'Height',
          'bt__cholesterol': 'Total cholesterol',
          'bt__albumin':"Albumin",
          'bt__sgot': "SGOT",
          'bt__alkaline_phosphatase':"Alkaline phosphatase",
          'bt__tsh':'TSH',
          'bt__inr':'INR',
          'bt__gtt':'GTT',
          'type_2_diabetes':'Type 2 diabetes',
          'currently_smokes':'Currently smokes',
          'ever_smoked':'Ever smoked',
          'stress': 'Stress',
          'bt__bilirubin':'Bilirubin'
          }.keys()

def explained_variance_binary(data,y_name,x_name):
    data['alpha'] = (data['alpha'] - data['alpha'].mean()) / data['alpha'].std()
    res = logit(formula=y_name + ' ~ '+x_name,data=data).fit()
    beta = res.params[x_name]
    alpha_explained_variance = beta**2 / (beta**2 + pi**(2./3))
    return alpha_explained_variance

def explained_variance(data,y_name,x_name,binary=False):
    if binary:
        return explained_variance_binary(data, y_name, x_name)

    data = (data - data.mean()) / data.std()
    try:
        res = ols(formula=y_name + ' ~ '+x_name, data=data).fit()
    except Exception:
        print data.head()
        raise Exception('%s - Bad phenotype for some reason '%y_name)
    if res.params['Intercept'] > 0.01:
        raise Exception('%s - Interception no zero ' % y_name)
    residual = data[y_name] - res.params[x_name] * data[x_name]
    alpha_explained_variance = (data.shape[0] - np.sum(residual * residual)) / data.shape[0]
    return alpha_explained_variance

def bootstrap_explained_variance(data,phenotype,clean_age_sex=True,times = 1000,binary=False):
    if binary or clean_age_sex:
        res = ols(formula='alpha ~ age + gender ', data=data).fit()
        alpha_clean = data['alpha'] - res.params['age'] * data['age'] - \
                         res.params['gender'] * data['gender']- res.params['Intercept']
        clean_data = pd.DataFrame(index=data.index, data={phenotype: data[phenotype], 'alpha': alpha_clean})
    else:
        res = ols(formula='alpha ~ gender ',data=data).fit()
        alpha_clean = data['alpha'] - res.params['gender']*data['gender'] - res.params['Intercept']
        clean_data = pd.DataFrame(index=data.index, data={phenotype: data[phenotype], 'alpha': alpha_clean})

    estimates = []
    for i in range(times):
        selected_ids = choice(clean_data.index, clean_data.shape[0])
        selected_data = clean_data.loc[selected_ids].copy()
        estimates.append(explained_variance(selected_data, phenotype,'alpha',binary))
    alpha_explained_variance = explained_variance(clean_data, phenotype, 'alpha',binary)
    estimates.sort()
    CI = [estimates[int(times*0.025)],estimates[int(times*0.975)]]
    return alpha_explained_variance,CI

def calc_alpha_explained_variance(phenotype,phenofile,alphafile,outputfile):
    phenotypes_df = pd.read_csv(phenofile).set_index('client_id')
    alpha = pd.Series.from_csv(alphafile)
    assert (alpha.index == phenotypes_df.index).all()
    res = pd.DataFrame(index=[phenotype], columns=['Explained variance', 'Low', 'High'])
    res.index.name = 'Phenotype'
    phenotype_data = phenotypes_df[phenotype].dropna()
    if len(set(phenotype_data))>2:
        is_outlier = np.abs(phenotype_data - phenotype_data.mean()) > 3 * phenotype_data.std()
        if (np.sum(is_outlier) > 0): phenotype_data.loc[is_outlier] = np.nan
        phenotype_data = phenotype_data.dropna()
        binary = False
    else:
        binary = True
    if phenotype=='age':
        clean_age_gender=False
        data = pd.DataFrame(index=phenotype_data.index,
               data={phenotype: phenotype_data, 'alpha': alpha[phenotype_data.index],
                     'gender': phenotypes_df.loc[phenotype_data.index, 'gender']})
    elif phenotype == 'gender':
        return
    else:
        clean_age_gender=True
        data = pd.DataFrame(index = phenotype_data.index,
                data = {phenotype: phenotype_data, 'age': phenotypes_df.loc[phenotype_data.index, 'age'],
                'gender': phenotypes_df.loc[phenotype_data.index, 'gender'],
                'alpha': alpha[phenotype_data.index]})
    alpha_explained_variance,CI=bootstrap_explained_variance(data,phenotype,clean_age_gender,binary=binary)
    res.loc[phenotype, 'H2'] = alpha_explained_variance
    res.loc[phenotype, 'CI_low'] = CI[0]
    res.loc[phenotype, 'CI_high'] = CI[1]
    res.to_csv(outputfile)

def main(cohort):
    phenofile = os.path.join(basepath, 'all_phenotypes_with_age_gender__%s.csv'%cohort)
    alphafile = os.path.join(basepath, 'Alpha-Shannon__%s.csv'%cohort)
    outputfile = os.path.join(basepath, 'Alpha-explained-variance-%s'%cohort + '-%s.csv')
    phenotypes_df = pd.read_csv(phenofile).set_index('client_id')
    alpha = pd.Series.from_csv(alphafile)
    assert (alpha.index == phenotypes_df.index).all()

    sethandlers()
    os.chdir(jobs_output)
    print ("Starting")
    with fakeqp(jobname='alpha2_', q=['himem7.q'], tryrerun=False, mem_def='1G', trds_def=2, \
            qworker='~/DevelopePycharm/lib/SegalQueue/qworker.py') as q:
        q.startpermanentrun()
        waiton = []
        for phenotype in keephenotypes:
            waiton.append(q.method(calc_alpha_explained_variance,
                            (phenotype,phenofile,alphafile,outputfile%phenotype)))
        q.wait(waiton)
    print ("Finished!")

if __name__=='__main__':
    for cohort in ['us','il__il_validation']:
        main(cohort)
        results=[]
        for phenotype in glob.glob(os.path.join(basepath, 'Alpha-explained-variance-%s*.csv'%cohort)):
            results.append(pd.read_csv(phenotype))
        pd.concat(results).set_index('Phenotype').sort_values('Explained variance',ascending=False).to_csv(os.path.join(basepath, 'Alpha-explained-variance-%s.csv'%cohort))
