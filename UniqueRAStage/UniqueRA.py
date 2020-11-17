import configparser
import logging
import os
import sys
import glob
import time

from BuildDB.createDistantConfigFile import createDistantConfigFile
from UniqueRAStage import map_reads, analyse_maps, score_maps

DB_PATH = '/net/mraid08/export/genie/Microbiome/Data/Databases/URA_IndicesAndScores'

_log = logging.getLogger('URA')
CONCAT_LIM = 500


def _tryrm(f):
    try:
        os.remove(f)
    except:
        pass


def main(finput, out_dir, read_len, num_mapped_to_subsample=None, min_mapped_to_retain=1000000,
         databases_path=DB_PATH, run_type='LargeOrNewGenusSGBs', only_perfect=False, exp_th=20,
         min_partial=0.5, min_sc_best=-40, keep_intermediate_files=False, min_abund=10**-4,
         bowtie_mapper='bowtie2',samtools_exe='samtools',unique_n=100):
    """
    Computes relative abundance by unique reads
    :param finput: SE fastq file
    :param out_dir: output directory, for temporary files and final output
    :param read_len: length of sequencing reads, int
    :param num_mapped_to_subsample: If specified, we subsample this number of reads out of the total reads mapped, int
    :param min_mapped_to_retain: Minimum number of mapped reads in order to retain the sample, int
    :param databases_path: Location of the database of files for unique relative abundance, str
    :param run_type: Name of parent database folder, str
    :param only_perfect: Use only perfect maps of reads to bacteria (no errs), bool
    :param exp_th: Minimal value for expected bin fill (stack bins together to get to this value), float
    :param min_partial: Minimal proportion of bin to use, at ends of contigs, float
    :param min_sc_best: Minimum score of map top use (if not only perfect maps), float
    :param keep_intermediate_files: Whether to keep, or remove, intermidiate files used for estimation process, bool
    :param min_abund: threshold to cut minimum abundance estimation, float
    :param bowtie_mapper bowtie2 executable path
    :param samtools_exe samtools executable path
    :return: map_perc - the percenta of reads of the samples mapped to index
             num_reads - if negative, minus the number of reads that mapped, which did not pass the minimal needed
                         if positive, the number of reads used for computing relative abundances
    """
    if (num_mapped_to_subsample is not None) and (min_mapped_to_retain > num_mapped_to_subsample):
        min_mapped_to_retain = num_mapped_to_subsample
    index_path = os.path.join(databases_path, run_type)
    index_name = os.path.join(index_path, 'bowtie_index', 'index')
    index_names2length = os.path.join(index_path, 'all_contigs.txt')

    sample = os.path.splitext(os.path.basename(finput))[0]

    _log.info('Mapping {} reads to DB'.format(sample))
    if only_perfect:
        min_sc_best = 0
    out_sam_file, map_perc, num_mapped = \
        map_reads.run_sample(sample, finput, out_dir, index_name, min_sc_best,bowtie_mapper=bowtie_mapper,samtools_exe=samtools_exe)
    _log.info('Aligned {}% for {}'.format(map_perc, sample))

    if num_mapped < min_mapped_to_retain:
        num_used = 0
    elif num_mapped_to_subsample is None:
        num_used = num_mapped
    else:
        num_used = min(num_mapped_to_subsample, num_mapped)

    if num_used == 0:
        _log.info('File %s had too few mapped reads, removed from RA analysis.' % sample)
        _tryrm(out_sam_file)
        _tryrm(out_sam_file.replace('.sam', '.stats'))
        return map_perc, -num_mapped

    _log.info('Analyzing maps')
    analyse_maps.run_sample(sample, out_sam_file, index_names2length, num_mapped_to_subsample, num_mapped, out_dir,
                            keep_intermediate_files)

    _log.info('Scoring maps')
    out_score_file = os.path.join(out_dir, 'abundance_{}.csv.gz'.format(sample))
    if len(glob.glob(os.path.join(index_path, 'scores_%d' % read_len))) == 0:
        raise Exception("read_len %d is not supported" % read_len)
    scores_dict_files = os.path.join(index_path, 'scores_%d' % read_len, 'dict_scores_%s')
    bins_dict_files = os.path.join(index_path, 'scores_%d' % read_len, 'dict_bins_%s')
    score_maps.run_sample(sample, out_score_file, out_dir, scores_dict_files, exp_th, min_partial, bins_dict_files,
                          min_abund,unique_n)

    if not keep_intermediate_files:
        _log.info('Removing intermediate analysis file')
        _tryrm(out_sam_file)
        _tryrm(out_sam_file.replace('.sam', '.stats'))
        _tryrm(os.path.join(out_dir, sample + '.pkl'))

    _log.info('Done with sample {}'.format(sample))
    return map_perc, num_used

def get_parameters_for_main(configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    base_path = config['user_parameters']['base_path']
    bowtie_mapper = config['user_parameters']['bowtie_mapper']
    samtools_exe = config['user_parameters']['samtools_exe']
    exp_th=int(config['run_pipeline']['exp_th'])
    run_type = config['run_pipeline']['run_type']
    keep_intermediate_files=eval(config['run_pipeline']['keep_intermediate_files'])
    unique_n=int(config['run_pipeline']['unique_n'])
    read_len = int(config['user_parameters']['read_len'])
    in_path=config['user_parameters']['in_path']
    in_file_ext=config['user_parameters']['in_file_ext']
    samples = glob.glob(os.path.join(in_path, "*.%s" % in_file_ext))
    out_path = config['run_pipeline']['path']
    num_mapped_to_subsample = eval(config['run_pipeline']['num_mapped_to_subsample'])
    min_mapped_to_retain=int(config['run_pipeline']['min_mapped_to_retain'])

    for samp in samples:
        _log.info("Started %s at" % samp, time.ctime())
        map_perc, num_used = main(samp, out_path, read_len, num_mapped_to_subsample, min_mapped_to_retain,databases_path=base_path, run_type=run_type,
                                  bowtie_mapper=bowtie_mapper,samtools_exe=samtools_exe,exp_th=exp_th,keep_intermediate_files=keep_intermediate_files,
                                  unique_n=unique_n)
        if num_used > 0:
            _log.info("%s had %g percent mapping. %d reads where used to asses URA" % (samp, map_perc, num_used))
        else:
            _log.info("%s had %g percent mapping. %d reads where too few asses URA" % (samp, map_perc, -num_used))

if __name__=="__main__":
    configFile = createDistantConfigFile(forceRewrite=False)
    get_parameters_for_main(configFile)