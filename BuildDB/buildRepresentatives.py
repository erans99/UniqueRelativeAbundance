import os
import pandas as pd
from Bio import SeqIO
import time
import glob
import configparser
from multiprocessing import Pool


def renameInputFastaFiles(inp):
    SelectedSGBs, output_cores_dir, genomes_dir,base_path = inp
    print(SelectedSGBs)
    all_contig_lens = []
    all_contig_ids = []
    for idx in SelectedSGBs.index:
        sgb = SelectedSGBs.loc[idx, 'SGB']
        genome = SelectedSGBs.loc[idx, 'nameOfFile']
        print(genome)
        print(genomes_dir)
        matchingFasta = glob.glob(os.path.join(genomes_dir, '%s.fa' % genome))
        coreContigs = []
        if len(matchingFasta) != 1:
            print("no core found for %s under %s" % (sgb, genomes_dir))
            continue
        matchingFasta = matchingFasta[0]
        coreFasta = os.path.join(output_cores_dir, 'SGB_%s_%s.fa' % (sgb, genome))
        if not os.path.exists(coreFasta):
            count = 0
            for record in SeqIO.parse(matchingFasta, "fasta"):
                record.id = 'SGB_%s_c_%s_%s' % (sgb, count, genome)
                all_contig_ids.append(record.id)
                count += 1
                coreContigs.append(record)
                all_contig_lens.append(len(record))
            SeqIO.write(coreContigs, coreFasta, 'fasta')
    pd.Series(index=all_contig_ids, data=all_contig_lens).to_csv( \
        os.path.join(base_path, 'all_contigs.txt'), sep='\t', header=False)
    print("Done writing core chunk fasta")


def buildByCore(SelectedSGBs, output_fasta, genomes_dir):
    if os.path.exists(output_fasta):
        os.system("rm -f %s" % output_fasta)
    first = True
    assert len(SelectedSGBs) == len(set(SelectedSGBs))
    for sgb in SelectedSGBs:
        matchingFasta = glob.glob(os.path.join(genomes_dir, 'SGB_%s_*.fa' % sgb))
        if len(matchingFasta) != 1:
            print("no core found for %s under %s" % (sgb, genomes_dir))
            continue
        matchingFasta = matchingFasta[0]
        if first:
            os.system('cat %s > %s' % (matchingFasta, output_fasta))
            first = False
        else:
            os.system('cat %s >> %s' % (matchingFasta, output_fasta))
    print("Done writing one big fasta for bowtie index")


def run(SelectedSGBs, configFile):
    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(configFile)
    build_representatives = config['build_representatives']
    print("Initializing")

    if not os.path.exists(build_representatives['output_cores_dir']):
        os.makedirs(build_representatives['output_cores_dir'])
    num_threads = eval(build_representatives['num_threads'])
    basedir = build_representatives['base_path']
    os.chdir(basedir)
    print("Starting")
    print(time.ctime())
    p = Pool(num_threads)
    waiton = []
    chunk_size = eval(build_representatives['chunksize'])
    for chunk in range(0, len(SelectedSGBs), chunk_size):
        waiton.append((SelectedSGBs[chunk:chunk + chunk_size],
                       build_representatives['output_cores_dir'],
                       build_representatives['genomes_dir'],
                       build_representatives['base_path']))
    p.map(renameInputFastaFiles, waiton)

    buildByCore(SelectedSGBs['SGB'], build_representatives['output_fasta'],
                build_representatives['output_cores_dir'])

    print(time.ctime())


if __name__ == '__main__':
    run()