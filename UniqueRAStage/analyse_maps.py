import os
import time

import pandas
import pysam
import logging

import map_reads

_log = logging.getLogger('URA')


def _tryrm(f):
    try:
        os.remove(f)
    except:
        pass


def add_group(map_of_read):
    diff_sc = -1000
    aligned_id = map_of_read.reference_name
    if aligned_id is None:
        return None
    aligned_group = aligned_id[4:aligned_id.find('_c_')]
    score = map_of_read.get_tag("AS")
    pos1 = aligned_id.find('_c_')
    pos2 = aligned_id.find('_', pos1 + 3)
    try:
        contig = int(aligned_id[pos1 + 3:pos2])
        if '|' in aligned_id:
            st = aligned_id[:pos1] + aligned_id[pos2:aligned_id.find('|')]
        else:
            st = aligned_id[:pos1] + aligned_id[pos2:]
    except ValueError:
        contig = int(aligned_id[pos1 + 3:])
        st = aligned_id[:pos1]
    pos = map_of_read.reference_start
    grp = aligned_group
    if map_of_read.has_tag('XS'):
        diff_sc = map_of_read.get_tag('XS') - score
        if diff_sc == 0:
            return None
        if diff_sc > 0:
            _log.warning("WTF. of %s %s XS %d > AS %d" % (map_of_read.reference_name, map_of_read.query_name,
                                                    map_of_read.get_tag('XS'), score))
    strain_dsc = 'Strain perfect' if score==0 else 'Strain best'
    return [strain_dsc, grp, st, contig, pos, score, diff_sc]


def read_map_file(mapping_file, outputf):
    _log.info(("Starting work on %s at" % mapping_file) + time.ctime())
    try:
        mapiter = pysam.AlignmentFile(mapping_file, 'rb', check_header=False, check_sq=False)
    except ValueError:
        raise ValueError("mapping file is invalid %s" % mapping_file)
    res = []
    for map_of_read in mapiter.fetch(until_eof=True):
        r = add_group(map_of_read)
        if r is not None:
            res.append(r)
    res = pandas.DataFrame(res, columns=['type', 'grps', 'IDs', 'contigs', 'poses', 'best_score', 'diff_score'])
    _log.info('Got %d uniquely best mapped reads' % len(res))
    res.to_pickle(outputf)
    return res


def get_output(sample, out_path):
    return os.path.join(out_path, sample + '.pkl')


def run_sample(sample, insam_file, index_names2length, num_mapped_to_subsample, num_mapped, out_dir,
               keep_intermediate_files):
    out_bam_f = os.path.join(out_dir, "part_" + sample + ".bam")
    # -F 256 means that only a single copy (the best, or one of the best by the score field 'AS' will be outputted.
    #    From bowtie manual: "Each reported read or pair alignment beyond the first has the SAM secondary bit
    #                         (which equals 256) set in its FLAGS field"
    if num_mapped_to_subsample is not None and num_mapped >= num_mapped_to_subsample:
        prop = num_mapped_to_subsample / num_mapped
        map_reads.sam2bam(insam_file, bam_outfile=out_bam_f, sort_bam=False, remove_sam_infile=False,
                          samtools_view_args='-F 256 -s {} -t {}'.format(prop, index_names2length))
    else:
        map_reads.sam2bam(insam_file, bam_outfile=out_bam_f, sort_bam=False, remove_sam_infile=False,
                          samtools_view_args='-F 256 -t {}'.format(index_names2length))

    assert os.path.exists(out_bam_f)

    output_file = get_output(sample, out_dir)
    if not (os.path.exists(output_file) and os.stat(output_file).st_size > 0):
        read_map_file(out_bam_f, output_file)
    _log.info('Done analysing sample %s' % sample)

    if not keep_intermediate_files:
        _log.info('Removing intermediate map file')
        _tryrm(out_bam_f)
