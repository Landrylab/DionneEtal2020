# load required libraries
# (see ./requirements.txt for the required python packages.)
from rohan.dandage.io_dict import read_dict
from rohan.dandage.align.demultiplex import run_demupliplex
from rohan.dandage.align.deepseq_dms import get_mutmat
# load configuration file containing parameters for the analysis 
cfg=read_dict('configuration.yml')
# demultiplex samples
run_demupliplex(cfg)
# obtain sample-wise mutation matrices
get_mutmat(cfg)
