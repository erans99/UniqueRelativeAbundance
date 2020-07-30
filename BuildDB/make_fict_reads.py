import os
import time
import glob
import configparser

def run_create_read( f, read_len ):
    fin = open( f[0], "r")
    fout = open( f[1], "w")
    seq = ""
    save_name = ""
    cnt = [0,0]
    while True:
        l = fin.readline()
        if (len(l) == 0) or (l[0] == '>'):
            if seq != "":
                out_all_read_to_file( seq, fout, save_name, read_len )
                cnt[1] += 1
            if len(l) == 0:
                break
            save_name = l[1:-1]
            seq = ""
        else:
            cnt[0] += 1
            seq += l[:-1]
    fin.close()
    fout.close()
    print ("Finished genome %s" % ( f[0] ), time.ctime())
    return

def out_all_read_to_file( seq, f_out, save_name, len_read ):
    for i in range(len(seq)-len_read+1):
        f_out.write(( ">%s|pos%d\n" % ( save_name.split(' ')[0],i )))
        f_out.write( seq[i:i+len_read] +'\n')

def clean(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    make_fict_reads=config['make_fict_reads']
    fictfile=glob.glob(os.path.join(make_fict_reads['fict_path'],'SGB_%s_*.fa'%(sgb)))
    if len(fictfile)==1:
        os.system('rm -rf %s'%fictfile[0])
    else:
        print ('Did not find ficticitios file for sgb %s in %s'%(sgb,make_fict_reads['fict_path']))
    return

def run_sample(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    make_fict_reads=config['make_fict_reads']
    try:
        if not os.path.exists(make_fict_reads['fict_path']):
            os.makedirs(make_fict_reads['fict_path'])
    except:
        pass
    inputfasta=glob.glob(os.path.join(make_fict_reads['singles_path'],
                        'SGB_%s_*.fa'%sgb))
    if len(inputfasta)!=1:
        raise Exception('Did not find core fasta for %s in %s'%(sgb,make_fict_reads['singles_path']))
    inputfasta=inputfasta[0]
    outputfasta = os.path.join(make_fict_reads['fict_path'],os.path.basename(inputfasta))
    run_create_read([inputfasta,outputfasta], eval(make_fict_reads['read_len']))

if __name__=='__main__':
    pass

