#!/usr/bin/env python3
from setuptools import setup, Extension
from distutils import sysconfig

setup(name = 'bcl2fastq_pipeline',
       version = '0.3.1',
       description = 'bcl2fastq_pipeline',
       author = "Devon P. Ryan",
       author_email = "ryan@ie-freiburg.mpg.de",
       scripts = ['bin/bfq.py'],
       packages = ['bcl2fastq_pipeline','flowcell_manager'],
       include_package_data = False,
       install_requires = ['configparser',
                           'numpy',
                           'matplotlib',
                           'bioblend',
                           'gcf-tools'],
       dependency_links = ['git+https://github.com/gcfntnu/gcf-tools#egg=gcf-tools'],
    )
