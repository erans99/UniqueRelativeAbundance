import configparser
import os
import pandas as pd

def eatOrKeepDistantGenus(GenusInfo,mash_path,rep_to_sgb_path,SGBdescription_path):
    allGenomes=pd.read_csv(SGBdescription_path,index_col=0)
    GenusInfo = pd.read_csv(mash_path,index_col=0)
    mash_path = os.path.join(GenusInfo)
    rep_to_sgb = pd.read_csv(rep_to_sgb_path, delim_whitespace=True, index_col=0, header=None)[1]
    mash_df = pd.read_csv(mash_path, index_col=0)
    new_names = [rep.replace('.fa', '') for rep in mash_df.index]
    mash_df.columns = new_names
    mash_df.index = new_names
    mash_df=mash_df.loc[rep_to_sgb.index,rep_to_sgb.index]
    new_names = [rep_to_sgb.loc[rep] for rep in mash_df.index]
    mash_df.columns = new_names
    mash_df.index = new_names
    sgbCounts=allGenomes.loc[new_names,'# Reconstructed genomes']
    largeSGBs=sgbCounts[sgbCounts>4].index.values
    print ("We have %s large SGBs"%len(largeSGBs))
    genusOfSGBs=GenusInfo.loc[rep_to_sgb.index]
    genusOfSGBs=genusOfSGBs.set_index('SGB')
    generaOfLargeSGBs=genusOfSGBs.loc[largeSGBs,'genus']
    smallSGBcounts=sgbCounts.drop(largeSGBs)
    tooClose=[]
    count=0
    for sgb in smallSGBcounts.index:
        if genusOfSGBs.loc[sgb,'genus'] in generaOfLargeSGBs.values:
            tooClose.append(sgb)
        else:
            count+=1
    print ("Added %s SGBs as they come from new genera"%count)
    smallSGBcounts = smallSGBcounts.drop(tooClose)
    annotation={k:[allGenomes.loc[k,'Estimated taxonomy'].split('|')[-1]] for k in smallSGBcounts.index}
    res = pd.DataFrame(index=smallSGBcounts.index,columns=['Assigned genomes','genus','annotation','minMashFromLarge','groupedSGBs','Grouped annotation'])
    res['Assigned genomes']=smallSGBcounts.values
    res['genus']=genusOfSGBs.loc[smallSGBcounts.index,'genus']
    res['annotation'] = [annotation[sgb] for sgb in res.index]
    res['minMashFromLarge']=mash_df.loc[smallSGBcounts.index, largeSGBs].min(axis=1)
    keepSGBs=[]
    keepSGBsAnnotations=[]
    groupedSGBs=[]
    for genus,genus_df in res.groupby('genus'):
        maxAssignedGenomes=genus_df['Assigned genomes'].max()
        maxOptions=genus_df[genus_df['Assigned genomes']==maxAssignedGenomes].index.values
        location = res.loc[maxOptions, 'minMashFromLarge'].argmax()
        keepID = res.loc[maxOptions,'minMashFromLarge'].index[location]
        keepSGBsAnnotations.append(genus_df['annotation'].values)
        keepSGBs.append(keepID)
        groupedSGBs.append(genus_df.index.values)
    res = res.loc[keepSGBs]
    res.loc[keepSGBs,'Grouped annotation'] = keepSGBsAnnotations
    res.loc[keepSGBs,'groupedSGBs'] = groupedSGBs
    return res,list(res.index)+list(largeSGBs)

def run(configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    LargeOrNew = config['LargeOrNewGenus']
    SelectedSGBs,allSGBs = eatOrKeepDistantGenus(LargeOrNew['representatives'],LargeOrNew['SGBdescription'],
                        LargeOrNew['GenusInfo'],LargeOrNew['mash_path'])
    pd.Series(allSGBs).to_csv(LargeOrNew['all_large_or_new_sgbs'],header=False,index=False)

if __name__=='__main__':
    run()