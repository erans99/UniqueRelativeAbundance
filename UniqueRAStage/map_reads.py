import glob
import os
import logging
import config

_log = logging.getLogger('URA')


def _tryrm(f):
    try:
        os.remove(f)
    except:
        pass


def sam2bam(sam_infile, bam_outfile=None, sort_bam=False, index_bam=False, remove_sam_infile=False,
            samtools_view_args=''):
    """
    Uses samtools to create a BAM file from a SAM file
    :param sam_infile: str, input SAM file
    :param bam_outfile: str, optional. If not provided, replaces .sam with .bam in sam_infile
    :param sort_bam: bool. Whether to sort the output BAM file
    :param index_bam: bool. Whether to create an index for the bam file
    :param remove_sam_infile: bool. Whether to remove the input sam file
    :param samtools_view_args: str. Additional arguments to samtools view
    :return: Name of bam output file
    """
    samtools_exe = config.SAM_TOOLS
    if bam_outfile is None:
        bam_outfile = sam_infile.replace('.sam', '') + '.bam'
    sort_bam = '| {} sort - '.format(samtools_exe) if sort_bam else ''
    cmd = '{} view -bS {} {} {} > {}'.format(samtools_exe, samtools_view_args, sam_infile, sort_bam, bam_outfile)
    os.system(cmd)
    if index_bam:
        cmd = '{} index {}'.format(samtools_exe, bam_outfile)
        os.system(cmd)
    if remove_sam_infile:
        _tryrm(sam_infile)
    return bam_outfile


def sam_fname_to_bam_fname(sam_infile):
    return (sam_infile[:-4] if sam_infile.endswith('.sam') else sam_infile) + '.bam'


def map(in_index_fname, in_fastq, out_sam_fname=None, out_std_err_fname=None, flags = "",
        output_bam=False, sort_bam=False, sam2bam_args='', remove_sam=True, inp_is_fasta=False):
    """
    Maps an input fastq file to a given index
    :param in_index_fname: pre-built bowtie index to map to
    :param in_fastq: str. Input fastq file
    :param out_sam_fname: str, optional. Output alignment SAM file. Default: in_fastq_1.sam (replace .fastq)
    :param out_std_err_fname: str, optional. Outout file for std output and error
    :param flags: str, optional. flags to add to bowtie mapping.
    :param output_bam: bool. If true output BAM instead of SAM
    :param sam2bam_args: str. Additional arguments in case of SAM to BAM conversion
    :param remove_sam: bool. Whether to remove sam file in case of SAM to BAM conversion
    :param inp_is_fasta: bool. Input file is fasta, not fastq
    :return: SAM / BAM alignment file name
    """
    bowtie_mapper = config.BOWTIE_MAPPPER
    if out_sam_fname is None:
        out_sam_fname = os.path.join(os.path.splitext(in_fastq)[0]) + '.sam'

    output_file = out_sam_fname
    if out_sam_fname.endswith('.bam') and output_bam:
        out_sam_fname = out_sam_fname[:-4] + '.sam'

    if out_std_err_fname is None:
        out_std_err_fname = os.path.join(os.path.splitext(in_fastq)[0]) + '.stats'
        _tryrm(out_std_err_fname)

    in_fastq = '-U {}'.format(in_fastq)
    out = '' if out_std_err_fname is None else '>& {}'.format(out_std_err_fname)
    if inp_is_fasta:
        cmd = '{} -f --threads {:d} {} -x {} {} -S {} {}'.format(bowtie_mapper, 2, flags,
                                                                 in_index_fname, in_fastq, out_sam_fname, out)
    else:
        cmd = '{} -q --threads {:d} {} -x {} {} -S {} {}'.format(bowtie_mapper, 2, flags,
                                                                 in_index_fname, in_fastq, out_sam_fname, out)
    os.system(cmd)

    if output_bam:
        output_file = sam2bam(out_sam_fname, bam_outfile=sam_fname_to_bam_fname(out_sam_fname), sort_bam=sort_bam,
                              remove_sam_infile=remove_sam, samtools_view_args=sam2bam_args)

    return output_file


def getOutput(sample, out_dir):
    return os.path.join(out_dir, sample + ".bam")


def getSamOutput(sample, out_dir):
    return os.path.join(out_dir, sample + ".sam")


def getStats(sample, out_dir):
    return os.path.join(out_dir, sample + ".stats")


def getStatsRate(sample, out_dir):
    try:
        txt = open(getStats(sample, out_dir)).read()
        pos2 = txt.find("% overall alignment rate")
        pos1 = txt[:pos2].rfind('\n')
        return float(txt[pos1 + 1:pos2])
    except:
        return 0.


def getNumMapped(sample, out_dir):
    try:
        txt = open(getStats(sample, out_dir)).read()
        pos2 = txt.find("aligned exactly 1 time")
        pos1 = txt[:pos2].rfind('\n')
        num_mapped = int(txt[pos1 + 1:pos2].strip().split()[0])
        pos2 = txt.find("aligned >1 times")
        pos1 = txt[:pos2].rfind('\n')
        num_mapped += int(txt[pos1 + 1:pos2].strip().split()[0])
        return num_mapped
    except:
        return 0.


def run_sample(sample, finput, out_dir, map_index_fname, min_sc_best):
    if os.path.isdir(out_dir):
        _log.info("Output results directory exists. Writing to it")
    else:
        try:
            os.makedirs(out_dir)
        except OSError:
            pass

    try:
        open(finput) # in order to raise exception if the file does not exist, of is not readable
        out_sam_file = map(map_index_fname, finput, out_sam_fname=getSamOutput(sample, out_dir),
                           out_std_err_fname=getStats(sample, out_dir),
                           flags="-a --no-unal --no-sq --no-hd --score-min L,%d,0" % min_sc_best)
        map_perc = getStatsRate(sample, out_dir)
        num_mapped = getNumMapped(sample, out_dir)
        return out_sam_file, map_perc, num_mapped
    except IOError:
        _log.warning("can't open/map input file for %s" % sample)
        return False, 0., 0
