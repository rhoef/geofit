# -*- coding: utf-8 -*-
"""
setup.py

Distutils setup script for ellipse_fit
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'
__svn_id__ = '$Id$'

import os
import sys
import glob

from distutils.core import setup, Extension

data_files =[ (os.path.join(sys.prefix, 'share', 'ellipse_fit') ,
               glob.glob('data/*'))
              ]
setup(
    name='ellipse_fit',
    version = str(1.0),
    description = 'Package to fit ellipses',
    author = 'Rudolf Hoefler',
    author_email = 'rudolf.hoefler@gmail.com',
    maintainer = 'Rudolf Hoefler',
    maintainer_email = 'rudolf.hoefler@gmail.com',
    package_dir = {'ellipse_fit': 'src'},
    data_files = data_files,
    packages = ['ellipse_fit'],
    )
