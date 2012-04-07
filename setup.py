# -*- coding: utf-8 -*-
"""
setup.py

Distutils setup script for geofit
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'


import os
import sys
import glob

from distutils.core import setup

data_files =[ (os.path.join(sys.prefix, 'share', 'geofit') ,
               glob.glob('data/*'))
              ]
setup(
    name='geofit',
    version = str(1.0),
    description = 'Package to fit ellipses',
    author = 'Rudolf Hoefler',
    author_email = 'rudolf.hoefler@gmail.com',
    maintainer = 'Rudolf Hoefler',
    maintainer_email = 'rudolf.hoefler@gmail.com',
    package_dir = {'geofit': 'src'},
    data_files = data_files,
    packages = ['geofit'],
    )
