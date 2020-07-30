import pysam
import glob
import os
import time
import pandas
import configparser

def read_map_files( fs ):
    for f in fs:
        read_map_file( f[0], f[1] )
    return

def add_group( ID, reads ):
    index = ["None", "No perfect", "Strain perfect", "Species perfect", "Many perfect" ]
    old_score = 0
    contigs = []
    ids = []
    poses = []
    grps = []
    found_real = False
    for read in reads:
        aligned_id = read.reference_name
        if aligned_id is None:
            if len(reads) != 1:
                print ("WTF. no name in group of size %d for %s" % ( len(reads), ID ))
            return False, [index[0], [], [], [], [] ]
#        aligned_group = df_sp[df_sp.GCA == aligned_id[:aligned_id.find('_c_')]].iloc[0].SGB
        aligned_group = aligned_id[4:aligned_id.find('_c_')]
        score = read.get_tag("AS")
        if score <= old_score:
            old_score = score
        else:
            print ("WTF score went up")
        if score == 0:
            pos1 = aligned_id.find('_c_')
            contigs.append ( int(aligned_id[pos1 + 3:].split('_')[0] ))
            ids.append( aligned_id[:pos1])
            poses.append( read.reference_start )
            grps.append( aligned_group )
            if ( not found_real ) and ( ID == ( aligned_id + "|pos%d" % poses[-1] )):
                found_real = True
    if not found_real:
        if not ( 'N' in read.seq ):
            print ("WTF. Real not found %s (not because of unknowns)" % ID )
        return False, [ index[0], [], [], [], [] ]
    if len( grps ) == 0:
        return False, [ index[1], [], [], [], [] ]
    if len( grps ) == 1:
        return True, [index[2], grps[0], ids[0], contigs[0], poses[0]]
    if len(set(grps)) == 1:
        return False, [index[3], grps[0], ids, contigs, poses ]
    return False, [index[4], grps, ids, contigs, poses ]

def read_map_file( mappingFile, outputf,analyse_fict_maps ):
    print ("Starting work on %s at" % mappingFile, time.ctime())
    mapiter = pysam.AlignmentFile( mappingFile, 'r')
    cnt = [ 0, 0, 0, 0 ]
    previous_id = None
    reads = []
    stats = []
    res = []
    start_pos = 0
    end_pos = 0
    try:
        for read in mapiter.fetch():
            cnt[0] += 1
            current_id = read.query_name
            if ( current_id != previous_id ) and ( previous_id is not None ):
                cnt[1] += 1
                tmp = add_group( previous_id, reads )
                cnt[2] += tmp[0]
                cnt[3] += tmp[0]
                stats.append( tmp[1] )
                reads = []
                if ( cnt[3] == eval(analyse_fict_maps['base_part_len'])) or \
                        ( current_id[-5:] == '|pos0' ):
                    res.append( add_stats_entry( stats, start_pos, previous_id ) )
                    start_pos = int(current_id[current_id.find("|pos") + 4:])
                    stats = []
                    cnt[3] = 0
            previous_id = current_id
            reads.append( read )
    except IOError:
        raise Exception('Mapping file %s is corrupted'%mappingFile)

    cnt[1] += 1
    if len(stats) != 0:
        res.append( add_stats_entry( stats, start_pos, previous_id ))
    res = pandas.DataFrame(res, columns= [ 'strain', 'contig', 'start_pos', 'end_pos', 'last_read', \
                                           'tot_len', 'num_in_strain'] )
    res.to_csv( outputf[0] )
    res.to_pickle(outputf[1])
    return res

def getOutput(sample,out_path):
    outputfile = [os.path.join(out_path, sample + '.csv'),
                  os.path.join(out_path, sample + '.pkl')]
    return outputfile

def add_stats_entry( stats, start_pos, previous_id ):
    stats = pandas.DataFrame(stats, columns=['type', 'grps', 'IDs', 'contigs', 'poses'])
    pos1 = previous_id.find('_c_')
    pos2 = previous_id.find('|pos')
    contig = int(previous_id[pos1 + 3:pos2].split('_')[0])
    id = previous_id[:pos1]
    end_pos = int(previous_id[previous_id.find("|pos") + 4:])
    return [ id, contig, start_pos, end_pos, previous_id, \
                len(stats), len(stats[stats['type'] == 'Strain perfect'])]

def run_sample(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    analyse_fict_maps=config['analyse_fict_maps']
    try:
        if not os.path.isdir(analyse_fict_maps['out_path']):
            os.makedirs(analyse_fict_maps['out_path'])
    except:
        pass
    inputfile = glob.glob(os.path.join(analyse_fict_maps['map_path'],
                        'SGB_%s_'%sgb+analyse_fict_maps['bowtie_endings']))
    if len(inputfile)!=1:
        raise Exception('Did not find bowtie any bowtie file %s'%os.path.join(analyse_fict_maps['map_path'],
                        'SGB_%s_'%sgb+analyse_fict_maps['bowtie_endings']))
    inputfile = inputfile[0]

    outputfiles=getOutput(os.path.basename(inputfile)[:-4], analyse_fict_maps['out_path'])
    if not os.path.exists(outputfiles[1]):
        read_map_file(inputfile, outputfiles,analyse_fict_maps)
if __name__=="__main__":
    pass
