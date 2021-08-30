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
