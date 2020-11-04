from scipy.spatial.distance import squareform,pdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu
import seaborn as sns
def BrayCurtis(X_orig, is_log_abundance=True, zero_min_value=True):
    if is_log_abundance:
        X = 10 ** X_orig
    else:
        X = X_orig.copy()
    if zero_min_value: X[X_orig == np.min(X_orig)] = 0

    D = squareform(pdist(X, metric='braycurtis'))
    return D


def plot_distributions(IL,US,COMB,save_fig_path,fontsize=12):
    illustraitor_il = (34. / 255, 181. / 255, 115. / 255)
    illustraitor_il_validation = (41. / 255, 171. / 255, 226. / 255)
    illustraitor_us = (102. / 255, 45. / 255, 145. / 255)
    colors_rgb = [illustraitor_il, illustraitor_il_validation, illustraitor_us]
    #data=pd.read_csv(save_fig_path + '_df.csv')
    # IL=data.loc[data['x']==0]
    # US=data.loc[data['x']==1]
    # COMB=data.loc[data['x']==2]
    labels=["IL-IL","IL-US","US-US"]
    sns.distplot(IL, hist=False, kde=True, label=labels[0], color=colors_rgb[0], kde_kws={'linewidth': 1})
    sns.distplot(COMB, hist=False, kde=True, label=labels[1], color=colors_rgb[1], kde_kws={'linewidth': 1})
    sns.distplot(US, hist=False, kde=True, label=labels[2], color=colors_rgb[2], kde_kws={'linewidth': 1})

    plt.xlabel('Bray-Curtis dissimilarity', fontsize=fontsize - 2)
    plt.ylabel('Frequency', fontsize=fontsize - 2)
    #plt.yticks([0, 5])
    #plt.xticks([2, 4, 6, 8, 10])
    ax= plt.gca()
    #ax.set_yticklabels([0, 0.9])
    #ax.yaxis.set_label_coords(-0.33, 0.45)
    plt.legend(labels, #loc = 'upper right', bbox_to_anchor = (1.57, 1.39),
    frameon = False, fontsize = fontsize - 2)#, labelspacing = 0.2, handletextpad = 0.2)
    plt.savefig(save_fig_path+'_distribution4.png',format='png')
    plt.savefig(save_fig_path+'_distribution4.pdf',format='pdf')


def compare_us_il(train_il_path,test_il_path,test_us_path,save_bc_path=None,save_fig_path=None,toplot=False,link_path=None):
    if True:
        ils=pd.read_pickle(save_fig_path+'_ils2.pkl')
        uss=pd.read_pickle(save_fig_path+'_uss2.pkl')
        comb=pd.read_pickle(save_fig_path+'_comb2.pkl')
        description=pd.concat([pd.Series(ils).describe(), pd.Series(uss).describe(), pd.Series(comb).describe()],1)
        description.columns=['IL-IL','US-US','IL-US']
        description.to_csv(save_fig_path + '_describe2.csv')
    else:
        keep_ids = pd.read_csv(link_path).set_index('client_id')
        train_il = pd.read_csv(train_il_path).set_index('client_id')
        print(train_il.shape)
        train_il = train_il.loc[train_il.index.isin(keep_ids.index)]
        print(train_il.shape)
        test_il = pd.read_csv(test_il_path).set_index('client_id')
        print(test_il.shape)
        test_il = test_il.loc[test_il.index.isin(keep_ids.index)]
        print(test_il.shape)
        test_us = pd.read_csv(test_us_path).set_index('client_id')
        print(test_us.shape)
        test_us = test_us.loc[test_us.index.isin(keep_ids.index)]
        print(test_us.shape)
        print('Loaded data')
        all_il = pd.concat([train_il, test_il])
        #matched = pd.read_csv(os.path.join(os.path.dirname(save_fig_path),'matched_isr_age_gender_bmi.csv'),index_col=1)
        #all_il = all_il.loc[matched.index]
        # all_il =  test_il
        il_ids = all_il.index
        us_ils = test_us.index
        all_data = pd.concat([all_il, test_us])
        print('Concated data')
        all_data_sgb = all_data.loc[:, all_data.columns.str.startswith('k__')]
        bc_all = BrayCurtis(all_data_sgb,False)
        bc_all = pd.DataFrame(data=bc_all,index=all_data.index,columns=all_data.index)
        print('Calculated BC')
        #bc_all.to_csv(save_bc_path)
        print('Saved BC')
        ils_square=bc_all.loc[il_ids,il_ids].values#.flatten()
        ils=ils_square[np.triu_indices(ils_square.shape[0], 1)]
        print(ils)
        print('Done flattening IL')
        uss_square=bc_all.loc[us_ils,us_ils].values#.flatten()
        uss=uss_square[np.triu_indices(uss_square.shape[0],1)]
        print('Done flattening US')
        comb=bc_all.loc[il_ids,us_ils].values.flatten()
        print('Done flattening COMB')

    # x=[0 for i in ils]+[1 for i in uss]+[2 for i in comb]
    # y=list(np.concatenate([ils,uss,comb]))
    # data = pd.DataFrame({'x':x,'y':y})
        pd.Series(ils).to_pickle(save_fig_path+'_ils2.pkl')
        pd.Series(uss).to_pickle(save_fig_path + '_uss2.pkl')
        pd.Series(comb).to_pickle(save_fig_path + '_comb2.pkl')
    #data.to_csv(save_fig_path+'_df.csv')
        print('Done saving DF')
        stats1 = mannwhitneyu(ils, uss)
        print("should not necessarily be significantly different ", stats1)
        stats2 = mannwhitneyu(ils, comb)
        print("should hopefully be significant ", stats2)
        stats3 = mannwhitneyu(uss, comb)
        print("should hopefully be significant ", stats3)
        stats = [[np.mean(ils), np.std(ils), [stats1, stats2]], [np.mean(uss), np.std(uss), [stats1, stats3]],
                 [np.mean(comb), np.std(comb), [stats2, stats3]]]
        pd.concat([pd.Series(ils).describe(),pd.Series(uss).describe(),pd.Series(comb).describe()]).to_csv(
            save_fig_path + '_describe2.csv')
        pd.DataFrame(data=stats, columns=['mean', 'std', 'stats1,2'], index=['ILs', 'USs', 'comb']).to_csv(
            save_fig_path + '_stats2.csv')

    # ax =sns.boxplot(data['x'],data['y'],color='white',fliersize=0,whis=[5, 95],width=0.5)
    # ax.set_xticklabels(['IL-IL','US-US','IL-US'])
    # plt.ylabel('Bray Curtis diversity')
    # plt.savefig(save_fig_path+'.png')
    # plt.savefig(save_fig_path+'.pdf',format='pdf')
    #stats1=mannwhitneyu(data.loc[data['x']==0,'y'],data.loc[data['x']==1,'y'])

    if toplot:
        plot_distributions(ils,uss,comb,save_fig_path)

if __name__=='__main__':
    data_path='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/DataForRevision'
    train_il_path = os.path.join(data_path, 'ura_q_il.csv')
    test_il_path = os.path.join(data_path, 'ura_q_il_validation.csv')
    test_us_path = os.path.join(data_path, 'ura_q_us.csv')
    output_basepath='/net/mraid08/export/jafar/Microbiome/Analyses/Unicorn/Cohort_Paper/revision_Analyses/beta_diversity/'
    output_basepath=os.path.join(output_basepath,'full_cohort')
    bc_output_path=os.path.join(output_basepath,'bc_all.csv')
    fig_output_path=os.path.join(output_basepath,'bc_fig')
    keep_ids_path=os.path.join(data_path,'link_clientID_to_LabData.csv')
    compare_us_il(train_il_path, test_il_path, test_us_path,bc_output_path,fig_output_path,link_path=keep_ids_path)
    #plot_distributions(fig_output_path)



