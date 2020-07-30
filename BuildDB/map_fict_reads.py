import os
import glob
import time
import configparser

def run_bowtie(sample, inf,outf, bowtie_exe,ind_path,ind_name,
               flags = "-a -f" ):
    print ("starting %s  at %s" % ( sample, time.ctime()))
    ferr = outf.replace( '.sam', '.stats' )
    com = "%s %s -x %s/%s -U %s -S %s >& %s " % ( bowtie_exe, flags, ind_path, \
                                                  ind_name, inf, outf, ferr )
    print("\nRunning:")
    print (com)
    print ("at", time.ctime())
    os.system(com)
    print ('Finish to call Bowtie!!!')
    print ("ended %s at %s" % (sample, time.ctime()))

def getOutput(f,map_res_path):
    return os.path.join(map_res_path, f.split('/')[-1].replace(".fastq","")\
                        .replace(".fa","").replace(".gz", "") + ".sam")

def ToRun(f,map_res_path):
    outputInitial = f.split('/')[-1].replace(".fa", "")
    fout = os.path.join(map_res_path, outputInitial + ".sam")
    fstats = os.path.join(map_res_path, outputInitial + ".stats")
    if os.path.isfile(fout) and os.path.isfile(fstats):
        stats = open(fstats).read()
        if "overall alignment rate" in stats:
            return False
    return True

def clean(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    map_fict_reads = config['map_fict_reads']
    samFile=os.path.join(map_fict_reads['map_res_path'],'SGB_%s_*.sam'%sgb)
    mapfile=glob.glob(samFile)
    if len(mapfile)==1:
        os.system('rm -rf %s'%mapfile[0])
    else:
        print('Did not find sam file for sgb %s in %s'%(sgb,map_fict_reads['map_res_path']))
    return

def run_sample(sgb,configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    map_fict_reads = config['map_fict_reads']
    run_on_sgb(sgb,**map_fict_reads)

def run_on_sgb(sgb,ind_path,ind_name,map_res_path,
               fict_path,bowtie_exe,**kwargs):
    if os.path.isdir(map_res_path):
        print ("Map results directory exists. Writing to it")
    else:
        os.makedirs(map_res_path)
    f = glob.glob(os.path.join(os.path.join(fict_path, 'SGB_%s_'%sgb+'*')))[0]#".fastq.gz" ))
    if not ToRun(f,map_res_path):
        print ("Mapping file exists - not running on %s"%sgb)
        return False

    run_bowtie(sgb, f, getOutput(f,map_res_path), bowtie_exe, ind_path, ind_name)
    return True

if __name__=='__main__':
    pass