#!/usr/bin/env python
"""
Gene Cluster Analysis Tools
===========================

    Gene Cluster Analysis Tools are a set of methods to robustly perform 
    general queries of gene cluster libraries.
"""
#    Gene Cluster Tools for Python
#    Copyright (C) 2014 by
#    Thomas E. Gorochowski <tom@chofski.co.uk>
#    All rights reserved.
#    OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

import sys
if sys.version_info[:2] < (2, 6):
    m = "Python version 2.6 or later is required for Gene Cluster Tools (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

import csv
import gene_cluster_library as gcl

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'OSI Non-Profit OSL 3.0'
__version__ = '1.0'



