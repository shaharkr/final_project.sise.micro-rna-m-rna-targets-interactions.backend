import os
import sys
os.environ['PATH'] = '/home/efrco/.conda/envs/my_env/bin:/storage/modules/packages/anaconda3/bin:/storage/modules/bin:/storage/modules/packages/anaconda3/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/storage/modules/packages/matlab/R2019B/bin:/home/efrco/.local/bin:/home/efrco/bin'

from logbook import Logger, RotatingFileHandler, StreamHandler

from consts.global_consts import log_file


class PipelineLogger(Logger):
    def __init__(self, log_file):
        super().__init__()
        self.handlers.append(RotatingFileHandler(log_file, bubble=True))
        self.handlers.append(StreamHandler(sys.stdout))

logger = PipelineLogger(log_file)
