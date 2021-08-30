# UniqueRelativeAbundance
standalone version of URA 
This Git includes 4 directories:
* Analysis - code used for predictions and atlas building, given URA output and phenotype data.
* BuildDB - code to build your own URA indices and scores for running the UniqueRAStage, of your own set of species.
* TestingsDBs - the DBs for running and debugging UniqueRAStage, for a subset of 6 species.
* UniqueRAStage - the actual code of computing URA for a single end fastq file.


# Installation requirements
In this repository, the following are needed:
* Tools: python 3, bowtie2, samtools
* python packages for URA: logging, os, sys, glob, time, pandas, pysam, functools, pickle, math, configparser, multiprocessing, Bio, subprocess
* further packages for scripts under Analysis folder: skbio, numpy, scipy, fiesta, mne, datetime, statsmodels, socket, matplotlib

# URA example usage
Running URA requires two steps:
1) Building a reference DB - to be done only once
2) Running URA - for each sample

# Building a reference DB
You can build your own DB given a set of genomes.
Edit UniqueRelativeAbundance/BuildDB/config.txt
The config.txt file contains documentation for the different parameters there.
* base_path = *A DIRECTORY WHERE URA BUILD FILES WILL BE WRITTEN*
* bowtie_path = *BOWTIE EXECUTABLE PATH*
* genomes_dir= *PATH FOR FOLDER CONTAINING REPRESENTATIVE GENOMES FASTA FILES*

Under *genomes_dir* folder, you need to create a file named *"Representatives_to_SGBs.txt"*
This file links the names of the genomes to number indetifiers that you choose.
It needs to be tab seperated without a header

For example:
```cat Representatives_to_SGBs.txt```

```
GENOME_ID1    1
GENOME_ID1    2
GENOME_ID1    3
```

Once the *config.txt* file is ready, run from BuildDB folder:
```python pipeline_run.py```

# Running URA 
Edit *UniqueRelativeAbundance/UniqueRAStage/config.py* - this requires a URA reference DB
```
DB_PATH = './'  # put path to reference DBs here - <base_path> from UniqueRelativeAbundance/BuildDB/config.txt
DB = 'URA_DB'  # URA_DB is the default name if you build a new DB 
```    

Run *UniqueRA.py* on your sample. Here is a command example:

```python UniqueRAStage.py input_fastq_dir fastq_file read_length output_dir```
* <input_fastq_dir> input directory, where fastq file is located
* <fastq_file> can be a specific fastq file or an extention (usually 'fq' or 'fastq')
* <read_length> fastq read length used for building the reference DB - default DB build read length is 75
* <output_dir> output directory, for temporary files and final output
      
