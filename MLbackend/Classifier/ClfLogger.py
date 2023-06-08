import sys
from pathlib import Path
from logbook import Logger, RotatingFileHandler, StreamHandler
from consts.global_consts import ROOT_PATH

log_file = ROOT_PATH / "Results/clf_log.log"

class ClfLogger(Logger):
    def __init__(self, log_file):
        super().__init__()
        self.handlers.append(RotatingFileHandler(log_file, bubble=True))
        # self.handlers.append(StreamHandler(sys.stdout))

logger = ClfLogger(log_file)
