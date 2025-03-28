#!/usr/bin/python3.13

import os
import logging
from datetime import datetime

class Logger():
    """ Basic logging class """

    def __init__(self):
        pass

    def get_logger(log_level: str='INFO', logfile: bool=False, current_time: str=datetime.now().isoformat()):
        """ Returns a logger at a given level, option to write to a file instead of to console """
        LOG_LEVEL = logging.NOTSET
        match log_level.upper():
            case 'DEBUG':
                LOG_LEVEL = logging.DEBUG
            case 'INFO':
                LOG_LEVEL = logging.INFO
            case 'WARNING':
                LOG_LEVEL = logging.WARNING
            case 'ERROR':
                LOG_LEVEL = logging.ERROR
            case 'CRITICAL':
                LOG_LEVEL = logging.CRITICAL
            case _:
                pass

        logger = logging.getLogger(__name__)
        if logfile:
            filename = f'URAProteomeMatcher_LOG_{current_time}.log'
            logging.basicConfig(filename=filename, encoding='utf-8', level=LOG_LEVEL)
        else:
            logging.basicConfig(format='%(levelname)s:%(message)s', level=LOG_LEVEL)

        return logger
