import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection


# rel_p - pathways/modules lists
# rel_kos_df  - ko's X pathways/modules matrix
#                        (1- ko is in the pathway/module; 0 - ko is not in pathway/module)
# temp_K0_res - a matrix that shows ko's rank of p-value associated with a phenotype
#               (i.e association between phenotype and a ko)
def test_enrichment(temp_K0_res, rel_kos_df, rel_p):
    ko_enrichment_res = pd.DataFrame(index=rel_p, columns=['pval', 'natural_pval', 'fdr', 'in_fraction_of_KOs'])

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

        ko_enrichment_res.loc[pathway, 'in_fraction_of_KOs'] = rel_kos_df[pathway].sum() / rel_kos_df.shape[0]
    ko_enrichment_res['fdr'] = fdrcorrection(ko_enrichment_res['natural_pval'])[1]
    ko_enrichment_res = ko_enrichment_res.sort_values('fdr')
    return ko_enrichment_res


if __name__ == "__main__":
    phenos = ['hba1c', 'bmi', 'age', 'bt__hdl_cholesterol', 'bt__fasting_triglycerides', 'bt__fasting_glucose']
    
    base_path = ""

    all_modules_list = pd.read_csv(os.path.join(base_path, 'kegg_db', 'list_of_modules.csv'))[
        ['module_no', 'name', 'super_module', 'sub_super_module']]
    all_modules_list.set_index('module_no', inplace=True)
    
    # Creating a binary df from kegg genes and pathways/modules
    gene_to_path = pd.read_csv(os.path.join(base_path, 'kegg_db', 'gene_to_path.csv'))
    pathways = pd.get_dummies(gene_to_path.set_index('kegg_gene')['path']).max(level=0).reset_index()

    gene_to_module = pd.read_csv(os.path.join(base_path, 'kegg_db', 'gene_to_module.csv'))
    modules = pd.get_dummies(gene_to_module.set_index('kegg_gene')['module']).max(level=0).reset_index()

    cols = {'KEGG_KO': 'kegg_gene', 'Spearman correlation P value': 'spear_pval',
            'Spearman correlation coefficient': 'spear'}
    # RankSum Mann-Whitney Test for pathway/module enrichment for different phenotypes
    all_path = {}
    all_module = {}
    for pheno in phenos:
        kegg = pd.read_csv(os.path.join(base_path, 'final_atlas_5_%s.csv' % pheno), skipinitialspace=True, index_col=0)
        kegg.rename(columns=cols, inplace=True)

        kegg['pval_log_sign'] = -(np.log(kegg.spear_pval) * np.sign(kegg.spear))
        kegg.sort_values('pval_log_sign', ascending=False, inplace=True)

        kegg['rank'] = range(kegg.shape[0])
        kegg = kegg[['kegg_gene', 'rank', 'pval_log_sign']].set_index('kegg_gene')

        rel_kos_pathways_df = pathways.copy()
        rel_kos_modules_df = modules.copy()

        # partial pathways that appear in converted data - note its the same set for all phenotypes
        rel_kos_pathways_df = rel_kos_pathways_df[rel_kos_pathways_df.kegg_gene.isin(kegg.index)]
        rel_kos_pathways_df.set_index('kegg_gene', inplace=True)

        pathway_ko_enrichment_res = test_enrichment(kegg, rel_kos_pathways_df, pathways.columns[1:])
        print("For %s got %d pathways" % (pheno, len(pathway_ko_enrichment_res[pathway_ko_enrichment_res.fdr <= 0.05])))
        pathway_ko_enrichment_res.to_csv(os.path.join(base_path, 'insights_max', '%s_paths.csv' % pheno))
        all_path[pheno] = pathway_ko_enrichment_res[['fdr', 'pval']].T

        rel_kos_modules_df = rel_kos_modules_df[rel_kos_modules_df.kegg_gene.isin(kegg.index)]
        rel_kos_modules_df.set_index('kegg_gene', inplace=True)

        module_ko_enrichment_res = test_enrichment(kegg, rel_kos_modules_df, modules.columns[1:])
        print("For %s got %d modules" % (pheno, len(module_ko_enrichment_res[module_ko_enrichment_res.fdr <= 0.05])))
        module_ko_enrichment_res = pd.merge(module_ko_enrichment_res, all_modules_list, left_index=True,
                                            right_index=True)
        module_ko_enrichment_res.to_csv(os.path.join(base_path, 'insights_max', '%s_module.csv' % pheno))
        all_module[pheno] = module_ko_enrichment_res[['fdr', 'pval']].T
all_path = pd.concat(all_path).T
for pheno in phenos:
    all_path[(pheno, 'd_pass')] = (all_path[(pheno, 'fdr')] < 0.05) * (-1)**(all_path[(pheno, 'pval')] < 0)
all_path = all_path.T.sort_index().T
all_path['num_phenos'] = (all_path[all_path.columns[all_path.columns.get_level_values(1) == 'fdr']] < 0.05).sum(1)
all_path.sort_values("num_phenos", ascending=False, inplace=True)
all_path.to_csv(os.path.join(base_path, 'insights_max', 'merged_pathways.csv'))
all_path[all_path.num_phenos >= 3][all_path.columns[all_path.columns.get_level_values(1) == 'd_pass']].to_csv(
    os.path.join(base_path, 'insights_max', 'pathways_for_figure_6.csv'))

all_module = pd.concat(all_module).T
for pheno in phenos:
    all_module[(pheno, 'd_pass')] = (all_module[(pheno, 'fdr')] < 0.05) * (-1)**(all_module[(pheno, 'pval')] < 0)
all_module = all_module.T.sort_index().T
all_module['num_phenos'] = (all_module[all_module.columns[all_module.columns.get_level_values(1) == 'fdr']] < 0.05).\
    sum(1)
all_module.sort_values("num_phenos", ascending=False, inplace=True)
all_module[all_module.num_phenos >= 3][all_module.columns[all_module.columns.get_level_values(1) == 'd_pass']].to_csv(
    os.path.join(base_path, 'insights_max', 'modules_for_figure_6.csv'))

all_module['function'] = all_modules_list.loc[all_module.index, 'name'].values
all_module['super_module'] = all_modules_list.loc[all_module.index, 'super_module'].values
all_module['sub_super_module'] = all_modules_list.loc[all_module.index, 'sub_super_module'].values
all_module.to_csv(os.path.join(base_path, 'insights_max', 'merged_modules.csv'))

print("Done")
