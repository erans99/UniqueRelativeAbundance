from skbio.diversity.alpha import shannon, simpson, berger_parker_d  
import pandas as pd
import os
import numpy as np
from scipy.stats import spearmanr,pearsonr,linregress
from Unicorn.Figures import nature_guidline_utils
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.stats import ranksums
from mne.stats import fdr_correction
sns.set_style("ticks", {'axes.edgecolor': 'black'})
# pd.set_option('display.width', 1000)
# np.set_printoptions(precision=4, linewidth=200)

params = {'axes.labelsize': 10,
          'axes.titlesize':8,
          'text.fontsize': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8}

basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Analyses/Alpha-Diversity'
phenotypes_thresholds=os.path.join(basepath,'Figures - phenotype_ranges.csv')
presence_absence_th=0.0001

matplotlib.rcParams.update(params)
clean_age=False
def alpha_div(x, method ):
    return method((x.values*10**7).astype(int)) #TODO: 5M?

def compute_alpha_div( M_df, method ):
    alpha = M_df.apply( alpha_div, 1, args = ( method, ))
    return alpha

def prepare_data():
    # obj1 = pd.read_pickle(os.path.join(basepath,'segata_q_us.pkl'))
    obj = pd.read_pickle(os.path.join(basepath,'segata_q_il_validation.pkl'))
    # obj3 = pd.read_pickle(os.path.join(basepath,'segata_q_il_validation.pkl'))  # Take from real data
    # obj = pd.concat([obj2, obj3])#([obj1, obj2, obj3])
    # obj=pd.read_pickle(os.path.join(basepath,'segata_q_il.pkl'))
    obj.drop_duplicates(subset='client_id', keep='first', inplace=True)
    objbac = obj.loc[:, obj.columns.map(lambda x: ('|sSGB__' in x) | (x == 'client_id'))].set_index('client_id')
    objbac[objbac<presence_absence_th] = 0
    objbac.to_csv(os.path.join(basepath,'sgb__il_validation.csv'))
    obj.loc[obj['quality'] == True, 'quality'] = 1
    obj.loc[obj['quality'] == False, 'quality'] = 0
    obj.loc[:,
    [u'age',u'gender',u'hba1c', u'bmi', u'bt__hdl_cholesterol', u'bt__cholesterol', u'bt__fasting_glucose', u'bt__fasting_triglycerides',
    u'bowel_movement_frequency', u'bt__gtt', u'morbid_obesity', u'high_triglycerides', u'type_1_diabetes',
    u'type_2_diabetes', u'pre_diabetes', u'impaired_glucose_tolerance_or_impaired_fasting_glucose', u'gestational_diabetes',
    u'inflammatory_bowel_disease', u'crohns_disease', u'ulcerative_colitis', u'undetermined_colitis',
    u'pancreatic_disease(s)', u'bt__inr', u'bt__protein', u'bt__sgot', u'bt__albumin',
    u'bt__alkaline_phosphatase', u'bt__bilirubin', u'bt__sgpt', u'bt__tsh', u'currently_smokes', u'evening_hunger',
    u'ever_smoked', u'general_hunger', u'height', u'midday_hunger', u'morning_hunger', u'physical_activity_freq',
    u'physical_activity_mins', u'quality', u'regular_defecation', u'sleep_quality', u'sleeping_time',
    u'stress', u'wake_up_time', u'weight', u'work_activity',u'pancreatic_disease(s)',
     u'client_id']].set_index('client_id').to_csv(os.path.join(basepath,'all_phenotypes_with_age_gender__il_validation.csv'))

def calcMeanWindow(vector,s=1000,jump=500):
    current=0
    res = []
    while current<len(vector):
        res.append(np.mean(vector[current:current+s]))
        current = current + jump
    return res

def plot_clouds(phenotype,pheno_df=None):
    if pheno_df is None:
        pheno_df= pd.read_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__il__il_validation.csv'),index_col=0)
    pheno_df=pheno_df[['alpha',phenotype]]
    pheno_df=pheno_df.dropna()
    r2 = spearmanr(pheno_df['alpha'],pheno_df[phenotype])#(mean_alpha[phenotype], mean_alpha['alpha'])
    plt.scatter(pheno_df['alpha'],pheno_df[phenotype],color='gray',s=1)#mean_alpha[phenotype], mean_alpha['alpha'])
    window_alpha=calcMeanWindow(pheno_df['alpha'])
    window_phenotype=calcMeanWindow(pheno_df[phenotype])
    plt.plot(window_alpha,window_phenotype,color='b')
    plt.ylabel(phenotype)
    plt.xlabel('Alpha diveristy')# % m)
    plt.title(' Spearmann r=%.3f%% p=%.4f'%(r2[0]*100,r2[1]))
    plt.show()
    plt.savefig(os.path.join(basepath,'%s_alpha_cloud.png'%phenotype))
    plt.close()


limits = {'age':[18,90],
          'bmi':[10,40],
          'hba1c':[3,10],
          'bowel_movement_frequency':[0,5],
          'bt__fasting_glucose':[60,150],
          'bt__fasting_triglycerides':[50,200],
          'bt__hdl_cholesterol':[30,400],
          'alpha':[3,7]}

def params_for_subplots(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.tick_params(top=False, right=False, pad=2)
    ymin, ymax = ax.get_ylim()
    x1, x2 = 0, 9
    y, h, col = ymin + (1.05 + 2 / 100.) * (ymax - ymin), 4.0 / 100 * (ymax - ymin), 'k'
    ll = plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c=col)
    plt.text((x1 + x2) * .5, y + h, "***", ha='center', va='bottom', color=col)
    ll[0].set_clip_on(False)

def plot_alpha_deciles_vs_pheontypes(pheno_df=None):
    if pheno_df is None:
        pheno_df= pd.read_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__il__il_validation.csv'),index_col=0)
    pheno_df['alpha-decile']=pd.qcut(pheno_df['alpha'],10,labels=[str(x) for x in range(1,11)])
    pheno_df=pheno_df[['alpha-decile', 'age', 'bmi', 'hba1c', 'bt__fasting_glucose',
              'bt__fasting_triglycerides', 'bt__hdl_cholesterol', 'alpha']]
    # plt.figure()
    plt.figure(figsize=(nature_guidline_utils.two_columns()/3,
                               nature_guidline_utils.full_page()), dpi=300)
    ax_alpha = plt.gca()
    ax_alpha = plt.subplot2grid((pheno_df.shape[1]-1,1), (pheno_df.shape[1]-2, 0), rowspan=1, colspan=1)
    # plt.ylabel('alpha', rotation = 0, labelpad = 35)
    phenotype='alpha'
    ax_alpha.set_ylabel('Alpha\ndiversity')

    # decile_df = pheno_df[['alpha-decile',phenotype]].pivot_table(values=phenotype,
    #    index=pheno_df[['alpha-decile', phenotype]].index,
    #    columns='alpha-decile', aggfunc='first')
    # decile_df.boxplot()
    ax_alpha = sns.boxplot(x=pheno_df['alpha-decile'], y=pheno_df[phenotype],
                           color='white',fliersize=0)
    ax_alpha.set_xlabel('Alpha diversity decile')
    # ax_alpha.set_ylabel(phenotype, rotation=0, labelpad=35)
    ax_alpha.set_ylabel('Shannon\nindex')
    # params_for_subplots(ax_alpha)
    ax_alpha.spines['right'].set_visible(False)
    ax_alpha.spines['top'].set_visible(False)
    ax_alpha.set_title('')
    ax_alpha.tick_params(top=False, right=False, pad=2)
    pvals=[]
    for i,phenotype in enumerate(pheno_df.columns):
        if phenotype == 'alpha-decile' or phenotype=='alpha':
            continue
        # plt.figure(figsize=(nature_guidline_utils.two_columns() / 3,
        #                         nature_guidline_utils.full_page()), dpi=300)
        # plt.figure()
        # ax_p = plt.gca()
        ax_p = plt.subplot2grid((pheno_df.shape[1]-1,1), (i-1, 0), rowspan=1, colspan=1)
        decile_df=pheno_df[['alpha-decile', phenotype]].pivot_table(values=phenotype,
               index=pheno_df[['alpha-decile', phenotype]].index,
               columns='alpha-decile', aggfunc='first')
        print decile_df[['1','10']].describe()
        # decile_df.boxplot()
        res = ranksums(decile_df['1'].dropna(),decile_df['10'].dropna())
        pvals.append(res[1])
        print phenotype, res[1]
        ax_p = sns.boxplot(x=pheno_df['alpha-decile'], y=pheno_df[phenotype],
                           color='white',fliersize=0)
        # ax_p.set_ylabel(phenotype, rotation = 0, labelpad = 35)
        ax_p.set_ylabel(phenotype.replace('bt__','').replace('_','\n'))
        params_for_subplots(ax_p)
        # ax_p.spines['right'].set_visible(False)
        # ax_p.spines['top'].set_visible(False)
        # ax_p.set_title('')
        # ax_p.set_xlabel('')
        # ax_p.set_xticklabels([])
        # ax_p.tick_params(top=False, right=False, pad=2)
        #
        # ymin, ymax =  ax_p.get_ylim()
        # x1, x2 = 0, 9
        # y, h, col = ymin + (1.05 + 2 / 100.) * (ymax - ymin), 4.0 / 100 * (ymax - ymin), 'k'
        # ll = plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c=col)
        # plt.text((x1 + x2) * .5, y + h, "***", ha='center', va='bottom', color=col)
        # ll[0].set_clip_on(False)

    # for phenotype in res.keys():
    #     res[phenotype].axes.autoscale(enable=False, axis='y')
    #     res[phenotype].axes.spines['right'].set_visible(False)
    #     res[phenotype].axes.spines['top'].set_visible(False)
    #     res[phenotype].axes.set_ylim(limits[phenotype])
    #     res[phenotype].axes.set_ylabel(phenotype,rotation=0,labelpad=35,fontsize=6)
    #     res[phenotype].axes.set_title('')
    #     res[phenotype].axes.tick_params(top=False, right=False, pad=2,labelsize=6)
        plt.subplots_adjust(left=0.3)
        plt.savefig(os.path.join(basepath,'phenotype_alpha_decile.png'))
    qvals=fdr_correction(pvals)
    print qvals
    #pheno_df.groupby('alpha-decile').unstack().plot(type='bar')
    # plt.show()

def get_outliers_by_phenotype(phenotype_name,phenotypes_df):
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
        is_outlier = (phenotypes_df[phenotype_name] < lower).values | (phenotypes_df[phenotype_name] > upper).values
    else:
        log = False
        is_outlier = [False]*len(phenotypes_df[phenotype_name])
    return is_outlier,log
    # is_outlier = np.abs(df_pheno[phenotype_name] - df_pheno[phenotype_name].mean()) > 3*df_pheno[phenotype_name].std()



def report_alpha_diverity():
    phenotypes_df = pd.read_csv(os.path.join(basepath,'all_phenotypes_with_age_gender__il.csv')).set_index('client_id')
    MB_df = pd.read_csv(os.path.join(basepath,'sgb__il.csv')).set_index('client_id')
    assert (MB_df.index==phenotypes_df.index).all()
    for m in [ shannon]:#, simpson, berger_parker_d ]:
        alpha = compute_alpha_div(MB_df, m)
        print ("Using %s" % m.func_name )
        alpha.to_csv(os.path.join(basepath,'Alpha-Shannon__il.csv'))
        pheno_df = pd.DataFrame(index=phenotypes_df.index,
                                columns=['alpha']+ phenotypes_df.columns.values)
        pheno_df['alpha'] = alpha[pheno_df.index]
        pheno_df = pheno_df.sort_values('alpha')
        for phenotype in phenotypes_df.columns:
            phenotype_data=phenotypes_df[phenotype].dropna()

            is_outlier, log = get_outliers_by_phenotype(phenotype, phenotypes_df)
            if (np.sum(is_outlier) > 0): phenotypes_df.loc[is_outlier, phenotype] = np.nan
            phenotype_data=phenotype_data.dropna()
            if log:
                phenotype_data = np.log10(phenotype_data)
            pheno_df.loc[phenotype_data.index,phenotype]=phenotype_data
        pheno_df = pheno_df.dropna(axis=1,how='all')
        pheno_df =  pheno_df.loc[:,pheno_df.std(axis=0)>0]
        pheno_df.to_csv(os.path.join(basepath, 'Phenotype-Alpha-Shannon__il.csv'))
    return pheno_df


            # bins_phenotype=np.linspace(phenotype_data.min(),phenotype_data.max(),100)
            # bin_number=np.digitize(phenotype_data,bins_phenotype)
            # binned_df=pd.DataFrame(index=phenotype_data.index,columns=[phenotype,'bin_number','alpha'])
            # binned_df[phenotype]=phenotype_data.values
            # binned_df['bin_number']=bin_number
            # binned_df['alpha']=alpha[phenotype_data.index].values
            # mean_alpha=binned_df.groupby('bin_number').mean()
            # r1 = pearsonr(mean_alpha[phenotype], mean_alpha['alpha'])

if __name__=="__main__":
    # prepare_data()
    report_alpha_diverity()
    #plot_alpha_deciles_vs_pheontypes()
