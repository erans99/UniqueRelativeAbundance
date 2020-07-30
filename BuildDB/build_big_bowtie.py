import os
import subprocess
import time
import configparser
def build_index( ind_path,input_index,bowtie_exe,offrate,threads,ind_name):

    buildCommand='%s --offrate %d --threads %d %s %s/%s' % ( bowtie_exe, offrate, \
                                                             threads, input_index, \
                                                             ind_path, ind_name)
    print ("building")
    print ("%s" % buildCommand )
    status, output = subprocess.getstatusoutput(buildCommand)
    if status != 0:
        print(output)
        print("status:", status)
        raise Exception("failed")
    return

def run(configFile):
    print ("Initializing")
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    build_representatives = config['build_bowtie']
    print ("Starting")
    print (time.ctime())
    
    ind_path = build_representatives['ind_path']
    input_index = build_representatives['input_index']
    bowtie_exe = build_representatives['bowtie_exe']
    offrate = eval(build_representatives['offrate'])
    threads = eval(build_representatives['threads'])
    ind_name = build_representatives['ind_name']
    if os.path.isdir(ind_path):
        print ("Indexing directory exists. Writing to it")
    else:
        os.mkdir(ind_path)
    build_index( ind_path,input_index,bowtie_exe,offrate,threads,ind_name)
    print(time.ctime())

if __name__=="__main__":
    run()