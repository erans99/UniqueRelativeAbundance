import glob
import os
import pandas
import pickle
import configparser

def createBinsAndScores(st,unite_fict_maps):
    try:
        if not os.path.exists(unite_fict_maps['scores_path']):
            os.makedirs(unite_fict_maps['scores_path'])
    except FileExistsError:
        print('WTF')
    dictOutput=os.path.join(unite_fict_maps['bins_dict_path']%st)
    scoresOuput=os.path.join(unite_fict_maps['scores_dict_path']%st)
    redo = False
    if os.path.exists(dictOutput) and os.path.exists(scoresOuput):
        try:
            pickle.load(open(dictOutput))
            pickle.load(open(scoresOuput))
        except Exception:
            print ("Redoing %s" % st )
            redo = True
    if redo or (not os.path.exists(dictOutput) or \
       not os.path.exists(scoresOuput)):
        bins_dict = {}
        scores_dict = {}
        score_st = pandas.read_csv(os.path.join(unite_fict_maps['strains_info_path'], st + '.csv'), index_col=0)
        c = -1
        score_c = []
        try:
            for c in range(score_st.contig.max() + 1):
                score_c = score_st[score_st.contig == c]
                bins_dict[c] = [-1] + list(score_c.end_pos.values)
                scores_dict[c] = list(score_c.num_in_strain.values)
            pickle.dump(bins_dict, open(dictOutput, "wb"))
            pickle.dump(scores_dict, open(scoresOuput, "wb"))
        except TypeError:
            print ("Failed reading info file %s %s %s %s %s %s" % ( st, \
                        len(score_st), score_st.contig.max(),type(score_st.contig.max()), c, len(score_c)))
            raise Exception("Failed reading info file %s %s %s %s %s %s" % ( st, \
                        len(score_st), score_st.contig.max(),type(score_st.contig.max()), c, len(score_c)))

def run_sample(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    unite_fict_maps = config['unite_fict_maps']
    info=glob.glob(os.path.join(unite_fict_maps['strains_info_path'], 'SGB_%s_*.csv'%sgb))
    if len(info)!=1:
        raise Exception('WTF. Info file not found for %s'%sgb)
    info = info[0]
    name = os.path.basename(info)[:-4] #remove .csv
    createBinsAndScores(name,unite_fict_maps)
    return
