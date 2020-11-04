import numpy as np
import pandas as pd
import os,subprocess,sys,glob
import scipy.linalg as la
from PNPChip.PNPChip.ForPaper import kernel_utils
from PNPChip.PNPChip.albi.albi_lib import fiesta_lib
import datetime
from mne.stats import fdr_correction
from lib.SegalQueue.qp import qp, fakeqp
from lib.addloglevels import sethandlers
presence_absence_th = 0.0001
presence_absence_mode = False

cohort = 'all-il' ##all-il/us
input = 'us' if cohort == 'us' else 'il__il_validation'
output_dir='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/Biome-Association-Index-%s'%cohort
name='SGB_association_index'
grm_file = os.path.join(output_dir, name + '_01%s'%presence_absence_mode)
GCTA_EXE = '/net/mraid08/export/genie/Bin/gcta/gcta64'
GRM_DIR = output_dir
abundances_f =os.path.join(output_dir,'sgb__%s.csv'%input)
covariates_f = os.path.join(output_dir,'covariates__%s.csv'%input)
pheno_f = os.path.join(output_dir,'all_phenotypes_plus_high_hba1c__%s.csv'%input) #il__il_validation
gcta_output = os.path.join(output_dir,'gcta_%s.txt')
jobs_output = os.path.join(output_dir,'jobs')
phenotypes_thresholds=os.path.join(output_dir,'Figures - phenotype_ranges.csv')

def correct_binary(cohort):
    from scipy.stats import norm
    if cohort == 'IL':
        currently_smokes=0.09
        ever_smoked=0.4
        type_2_diabetes=0.12
        hba1c_high = 0.12
    if cohort == 'US':
        currently_smokes=0.02
        ever_smoked = 0.26
        type_2_diabetes = 0.05
        hba1c_high = 0.05
    Ks=[currently_smokes,ever_smoked,type_2_diabetes,hba1c_high]
    estimates=pd.read_csv(os.path.join(output_dir,'LMM_results.csv'),index_col=0)
    for i,estimate in enumerate(['currently_smokes','ever_smoked','type_2_diabetes','hba1c_high']):
        t = norm(0, 1).isf(Ks[i])
        phi_t = norm(0, 1).pdf(t)
        corrected_h2 = Ks[i] * (1 - Ks[i]) / phi_t ** 2 * (estimates.loc[estimate,'H2'])
        corrected_h2_ci_low= Ks[i] * (1 - Ks[i]) / phi_t ** 2 * (estimates.loc[estimate,'CI_low'])
        corrected_h2_ci_high= Ks[i] * (1 - Ks[i]) / phi_t ** 2 * (estimates.loc[estimate,'CI_high'])
        estimates.loc[estimate,'H2']=corrected_h2
        estimates.loc[estimate, 'CI_low'] = corrected_h2_ci_low
        estimates.loc[estimate, 'CI_high'] = corrected_h2_ci_high
    estimates.to_csv(os.path.join(output_dir,'LMM_results_corrected_binary_%s.csv'%cohort))
def _parse_GCTA_results(stdout, phenotype_name):
    intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std = None, None, None, None, None, None, None
    fixed_effects = []
    for line in stdout.split('\n'):
        line_split = line.split()
        if (len(line_split) == 0): continue
        if (line_split[0] == 'mean'): intercept = float(line_split[1])
        if (line_split[0][:2] == 'X_'): fixed_effects.append(float(line_split[1]))
        if (line_split[0] == 'V(G)'): sig2g = float(line_split[1])
        if (line_split[0] == 'V(e)'): sig2e = float(line_split[1])
        if ('error' in line.lower()):
            print ('GCTA', phenotype_name, line)
            return intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects
        if line.startswith('V(G)/Vp') and 'V(G)/Vp_L' not in line: h2_reml, h2_reml_std = float(line_split[1]), float(line_split[2])
        if line.startswith('V(G)/Vp_L'): h2_reml_l, h2_reml_l_std = float(line_split[1]), float(line_split[2])
    return intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects

def norm_func(x):
    return (x - x.mean()) / x.std()

def run_LMM():
    df_pheno = pd.read_csv(pheno_f,index_col=0)
    df_covariates = pd.read_csv(covariates_f,index_col=0)
    # remove individules with missing values
    assert (df_covariates.index==df_pheno.index).all()
    if not os.path.exists(grm_file+'.grm.gz') or \
       not os.path.exists(grm_file+'.grm.id') :
        print ("Calculating grm matrix %s"%datetime.datetime.now())

        df = pd.read_csv(abundances_f,index_col=0)
        # remove SGB in case 0 std
        df = df.apply(np.log10)
        keep = df.columns[df.std() != 0]
        df = df[keep]
        keep = df.index[df.std(axis=1) != 0]
        df = df.loc[keep]
        df_norm = df.apply(norm_func)
        K_mb = df_norm.dot(df_norm.T) / df_norm.shape[1]
        df_K = pd.DataFrame(K_mb, index=df.index, columns=df.index)
        df_K.to_csv(grm_file + '.csv')
        kernel_utils.write_grm(df_K.values, grm_file + '.grm.gz')
        df_K['IID'] = df_K.index
        df_K['IID'].to_csv(grm_file+ '.grm.id', sep='\t', index_label='IID')
        print ("Done Calculating grm matrix %s" %datetime.datetime.now())

    sethandlers()
    os.chdir(jobs_output)
    print ("Starting")
    with qp(jobname='b2_', q=['himem7.q'], tryrerun=False, mem_def='10G', trds_def=2, \
            qworker='~/DevelopePycharm/lib/SegalQueue/qworker.py') as q:
        q.startpermanentrun()
        waiton = []
        for phenotype_name in df_pheno.columns:
            waiton.append(q.method(run_biome_association_index,(phenotype_name,)))
        q.wait(waiton)
    print ("Finished!")

def get_outliers_by_phenotype(phenotype_name,df_pheno):
    phenotype_threshold = pd.read_csv(phenotypes_thresholds,index_col=0)
    if phenotype_name in phenotype_threshold.index:
        if phenotype_threshold.loc[phenotype_name,'Selection']=='Log':
            log=True
            lower=10**phenotype_threshold.loc[phenotype_name,'Log lower']
            upper = 10**phenotype_threshold.loc[phenotype_name, 'Log upper']
        else:
            log=False
            lower=phenotype_threshold.loc[phenotype_name,'Linear lower']
            upper = phenotype_threshold.loc[phenotype_name, 'Linear upper']
        is_outlier = (df_pheno[phenotype_name] < lower).values | (df_pheno[phenotype_name] > upper).values
    else:
        log = False
        is_outlier = [False]*len(df_pheno[phenotype_name])
    return is_outlier,log

def run_biome_association_index(phenotype_name):
    print ("Staring working on %s - time %s"%(phenotype_name,datetime.datetime.now()))
    df_pheno = pd.read_csv(pheno_f, index_col=0)
    df_covariates = pd.read_csv(covariates_f, index_col=0)
    if phenotype_name=='age':
        df_covariates=df_covariates.loc[:,['gender']]
    df_K = pd.read_csv(grm_file + '.csv',index_col=0)
    df_K.columns = df_K.columns.astype(int)
    df_K.index = df_K.index.astype(int)
    df_K['IID'] = df_K.index
    assert (df_K.index == df_pheno.index).all()

    print ("Finished loading GRM matrix")
    pheno_file = os.path.join(output_dir, '%s.phe'%phenotype_name)
    covariates_file = os.path.join(output_dir, '%s.cov'%phenotype_name)

    df_covariates['ID'] = df_covariates.index.values
    df_covariates = df_covariates[['ID'] + [c for c in df_covariates.columns if c != 'ID']]
    df_covariates.to_csv(covariates_file, sep='\t', header=False, na_rep='NA')
    del df_covariates['ID']

    with open(os.path.join(output_dir, '%s_LMM_results.txt'%phenotype_name), 'w') as handle:
        handle.write('\t'.join(['Phenotype', 'microbiome-association index', '95% CI', 'P value',
                                'Sample size', 'V(G)', 'V(e)', 'mean'] + list(df_covariates.columns)) + '\n')

    is_outlier,log=get_outliers_by_phenotype(phenotype_name,df_pheno)
    if (np.sum(is_outlier) > 0): df_pheno.loc[is_outlier, phenotype_name] = np.nan
    # check if there are enough samples
    if df_pheno[phenotype_name].notnull().sum() < 500:
        print (("Not enough samples", str(phenotype_name)))
        return

    if log:
        df_pheno_temp = np.log10(df_pheno[[phenotype_name]].copy())
    else:
        df_pheno_temp = df_pheno[[phenotype_name]].copy()
    df_pheno_temp['ID'] = df_pheno_temp.index.values
    df_pheno_temp[['ID', phenotype_name]].to_csv(pheno_file, sep='\t', header=False, na_rep='NA', index=True)
    del df_pheno_temp['ID']

    #invoke GCTA
    print ("Running GCTA for phenotype %s, time %s" % (phenotype_name,datetime.datetime.now()))
    cmdLine = [GCTA_EXE, '--grm-gz', grm_file, '--reml', '--pheno', pheno_file,
               '--qcovar', covariates_file, '--reml-est-fix', '--thread-num', '1']
    proc = subprocess.Popen(cmdLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    with open(gcta_output%phenotype_name,'w') as fname:
        fname.write(stdout)
    print ("Done Running GCTA for phenotype %s, time %s" % (phenotype_name, datetime.datetime.now()))

    if (stderr is not None):
        print ('GCTA reml error:')
        print (stderr)
        raise Exception()
    with open(gcta_output%phenotype_name, 'r') as fname:
        stdout = fname.read()

    intercept, sig2g, sig2e, h2_reml, h2_reml_std, h2_reml_l, h2_reml_l_std, fixed_effects = _parse_GCTA_results(stdout, phenotype_name)

    print ("Running RL-SKAT and fiesta for phenotype %s, time %s" % (phenotype_name, datetime.datetime.now()))
    if h2_reml is None:
        print (stdout)
        raise Exception('h2_reml is none')
    print ("Phenotype %s - b2 - %s"%(phenotype_name,h2_reml))
    #compute a p-value with RL-SKAT
    df_merged = df_K.merge(df_covariates.dropna(), left_index=True, right_index=True).copy()
    df_merged = df_merged.merge(df_pheno_temp[[phenotype_name]], left_index=True, right_index=True)
    df_merged = df_merged.loc[df_merged[phenotype_name].notnull()]
    covariates = df_merged[[c for c in df_covariates.columns if c != 'ID']].copy()
    covariates['intercept'] = 1.0
    covariates = covariates.values
    pheno_vec = df_merged[phenotype_name].values
    df_K_subset = df_K.loc[df_K.index.isin(df_merged.index), df_K.columns.isin(df_merged.index)]
    df_K_subset = df_K_subset.loc[df_merged.index, df_merged.index]
    assert (df_K_subset.index == df_K_subset.columns).all()
    assert (df_K_subset.shape[0] == df_merged.shape[0])
    assert (df_K_subset.index.isin(df_merged.index).all())
    assert (df_K_subset.index == df_merged.index).all()
    kernel = df_K_subset.values
    #pvalues = rl_skat.RL_SKAT_Full_Kernel(kernel, covariates, add_intercept=False).test(np.row_stack(pheno_vec))
    pvalues=[1]
    #compute CIs with FIESTA
    s,U = la.eigh(kernel)
    ind = np.argsort(s)[::-1]; s=s[ind]; U=U[:,ind]
    CI = fiesta_lib.calculate_cis_general(np.array([h2_reml]), s, U, covariates,
                                          iterations=100, alpha=0.05,
                                          tau=0.4, use_convergence_criterion=False,
                                          use_progress_bar=False) # iterations=1000
    line =  '\t'.join([str(c) for c in [
                     phenotype_name,
                     '%0.3g'%(h2_reml), '%0.3g - %0.3g'%(CI[0][0], CI[0][1]),
                     '%0.3g'%(pvalues[0]), pheno_vec.shape[0], '%0.3g'%(sig2g), '%0.3g'%(sig2e),
                     '%0.3f'%(intercept)] + ['%0.3g'%(f) for f in fixed_effects]]) + '\n'

    with open(os.path.join(output_dir, '%s_LMM_results.txt'%phenotype_name), 'a') as handle:
        handle.write(line)
    sys.stdout.flush()
    print ("Done Running RL-SKAT and fiesta for phenotype %s, time %s" % (phenotype_name, datetime.datetime.now()))

def prepare_data():
    objUS=pd.read_pickle(os.path.join(output_dir,'segata_q_us.pkl'))
    obj2=pd.read_pickle(os.path.join(output_dir,'segata_q_il.pkl'))
    obj3=pd.read_pickle(os.path.join(output_dir,'segata_q_il_validation.pkl'))
    objIL = pd.concat([ obj2, obj3])
    cohort=['IL','US']
    for idx_cohort,obj in enumerate([objIL,objUS]):
        obj.drop_duplicates(subset='client_id', keep='first', inplace=True)
        objbac = obj.loc[:, obj.columns.map(lambda x: ('|sSGB__' in x) | (x == 'client_id'))].set_index('client_id')
        objbac[objbac<presence_absence_th] = presence_absence_th
        objbac.to_csv(os.path.join(output_dir,'sgb__us.csv'))
        obj.loc[:,['age','gender','client_id']].set_index('client_id').to_csv(os.path.join(output_dir,
                                                                            'covariates__%s.csv'%cohort[idx_cohort]))
        obj.loc[:,['hba1c','bmi','client_id']].set_index('client_id').to_csv(os.path.join(output_dir,
                                                                            'hba1c_bmi__%s.csv'%cohort[idx_cohort]))
        obj.loc[:,
        [u'age',u'hba1c', u'bmi', u'bt__hdl_cholesterol', u'bt__cholesterol', u'bt__fasting_glucose', u'bt__fasting_triglycerides',
        u'bowel_movement_frequency', u'bt__inr', u'bt__protein', u'bt__sgot', u'bt__albumin',
        u'bt__alkaline_phosphatase', u'bt__bilirubin', u'bt__tsh', u'currently_smokes',
        u'height', u'weight', u'client_id']].set_index('client_id').to_csv(os.path.join(output_dir,
            'all_phenotypes__us.csv'))

def JoinAndParse():
    results = glob.glob(os.path.join(output_dir,'*_LMM_results.txt'))
    objs = [pd.read_csv(res, sep='\t') for res in results]
    obj = pd.concat(objs).set_index('Phenotype').sort_values('microbiome-association index', ascending=False)
    obj['H2'] = obj['microbiome-association index']
    obj['CI_low'] = obj['95% CI'].apply(lambda x: np.float(x.split(' - ')[0]))
    obj['CI_high'] = obj['95% CI'].apply(lambda x: np.float(x.split(' - ')[1]))
    obj['Q value'] = fdr_correction(obj['P value'])[1]
    obj[['H2', 'CI_low', 'CI_high', 'P value', 'Q value',
                       'Sample size', 'V(G)', 'V(e)', 'mean', 'age', 'gender']].to_csv(os.path.join(output_dir,
                                                                                'LMM_results.csv'))

if __name__=="__main__":
    run_LMM()
    JoinAndParse()
    correct_binary()