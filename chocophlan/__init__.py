from functools import partial
from lxml import etree 
from urllib.parse import urlparse
import datetime
import ftplib
import hashlib
import math
import multiprocessing as mp
import multiprocessing.dummy as mpdummy
import os
import pickle
import re
import requests
import sys
import tarfile
import tempfile
import time
import traceback
import glob
import itertools
import gzip
import bz2
import traceback
import shutil
import random
import shutil
import logging

def init_parse(terminating_):
    global terminating
    terminating = terminating_

def grouper(iterable, n, fillvalue=None):
    #"Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

