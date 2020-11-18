import configparser
import os
def createDistantConfigFile(forceRewrite=False,expandVars=True,configname='config.txt'):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config_file=os.path.join(os.getcwd(), configname)
    if expandVars:
        with open(config_file, 'r') as cfg_file:
            cfg_txt = os.path.expandvars(cfg_file.read())
        config.read_string(cfg_txt)
    else:
        config.read(config_file)
    if not os.path.exists(config['run_pipeline']['path']):
        os.makedirs(config['run_pipeline']['path'])
    outputConfigFile=os.path.join(config['run_pipeline']['path'], 'config.txt')
    if not os.path.exists(outputConfigFile) or forceRewrite:
        with open(outputConfigFile, 'w') as configfile:
            config.write(configfile)
    else:
        print("Using existing config file %s"%outputConfigFile)
    return outputConfigFile

