"""
Function of this scripts
"""

import logging
import coloredlogs

from iga.apps.base import emain

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

if __name__ == "__main__":
    emain()
