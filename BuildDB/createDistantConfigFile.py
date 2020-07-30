import configparser
import os
def createDistantConfigFile(forceRewrite=False):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(os.path.join(os.getcwd(), 'config.txt'))
    if not os.path.exists(config['run_pipeline']['qp_base_dir']):
        os.makedirs(config['run_pipeline']['qp_base_dir'])
    outputConfigFile=os.path.join(config['run_pipeline']['qp_base_dir'], 'config.txt')
    if not os.path.exists(outputConfigFile) or forceRewrite:
        with open(outputConfigFile, 'w') as configfile:
            config.write(configfile)
        ###TWF! DO NOT ERASE THE FOLLOWING
        # config.txt FILE IS OVERWRITTEN WITH QP if the following line is in __main__
        # del configfile
    else:
        print ("Using existing config file %s"%outputConfigFile)
    return outputConfigFile
