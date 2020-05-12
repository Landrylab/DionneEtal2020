from rohan.dandage.io_dict import read_dict
from rohan.dandage.align.demultiplex import run_demupliplex
from rohan.dandage.align.deepseq_dms import get_mutmat
cfg=read_dict('all_pooled_global_input_cfg.yml')
run_demupliplex(cfg)
get_mutmat(cfg)
