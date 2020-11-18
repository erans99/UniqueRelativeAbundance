import logging
import os
import sys
import glob
import time

import map_reads, analyse_maps, score_maps
import config

_log = logging.getLogger('URA')


def _tryrm(f):
    try:
        os.remove(f)
    except:
        pass


def main(finput, out_dir, read_len, num_mapped_to_subsample=None, min_mapped_to_retain=1000000,
         databases_path=config.DB_PATH, run_type=config.DB, num_uniq=config.NUM_UNIQ_IN_PART, only_perfect=False,
         exp_th=20, min_partial=0.5, min_sc_best=-40, keep_intermediate_files=False, min_abund=10**-4):
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
    :param num_uniq: number of uniq reads in base part of dictionaries, int
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
        map_reads.run_sample(sample, finput, out_dir, index_name, min_sc_best)
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
                          min_abund, num_uniq)

    if not keep_intermediate_files:
        _log.info('Removing intermediate analysis file')
        _tryrm(out_sam_file)
        _tryrm(out_sam_file.replace('.sam', '.stats'))
        _tryrm(os.path.join(out_dir, sample + '.pkl'))

    _log.info('Done with sample {}'.format(sample))
    return map_perc, num_used


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("python3 UniqueRAStage.py <in_dir> <in_file_ext> <read_length> <out_dir> <max_reads> <min_reads>")
        print("<in_dir> input directory, for fastq files")
        print("<in_file_ext> fastq file extention (usually 'fq' or 'fastq')")
        print("<read_length> length of reads in fastq (works with 75, 100 or 150)")
        print("<out_dir> output directory, for temporary files and final output")
        print("<max_reads> maximal number of mapped reads to use when computing abundances, if 0 use all, default 5e6")
        print("<min_reads> minimal number of mapped reads to use when computing abundances, ")
        print("            will not compute if there are not enough reads, default 1e6")
        sys.exit(0)
    else:
        in_path = sys.argv[1]
        samples = glob.glob(os.path.join(in_path, "*.%s" % sys.argv[2]))
        read_len = int(sys.argv[3])
        out_path = sys.argv[4]
        if len(sys.argv) > 5:
            num_mapped_to_subsample = int(eval(sys.argv[5]))
            if num_mapped_to_subsample == 0:
                num_mapped_to_subsample = None
                print("Will not subsample")
        else:
            num_mapped_to_subsample = 5 * 10**6
        if len(sys.argv) > 6:
            min_mapped_to_retain = int(eval(sys.argv[6]))
        else:
            min_mapped_to_retain = 10**6

    for samp in samples:
        _log.info("Started %s at" % samp, time.ctime())
        map_perc, num_used = main(samp, out_path, read_len, num_mapped_to_subsample, min_mapped_to_retain)
        if num_used > 0:
            _log.info("%s had %g percent mapping. %d reads where used to asses URA" % (samp, map_perc, num_used))
        else:
            _log.info("%s had %g percent mapping. %d reads where too few asses URA" % (samp, map_perc, -num_used))
