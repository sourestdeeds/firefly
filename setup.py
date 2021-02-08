#!/usr/bin/env python3

import os
import codecs
import pathlib
import setuptools

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()
        
def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            # __version__ = "0.9"
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")
    
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup( name='firefly',
                  version=get_version('firefly/__init__.py'),
                  description='For use with TransitFit',
                  author='Steven Charles-Mindoza',
                  license='MIT',
                  long_description=long_description,
                  long_description_content_type='text/markdown',
                  packages=setuptools.find_packages(),
                  classifiers=[
                    "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License"
                  ],
                  install_requires=['numpy',
                                    'pandas', 
                                    'tabulate',
                                    'fuzzywuzzy',
                                    'python-Levenshtein',
                                    'natsort',
                                    'seaborn',
                                    'astroquery',
                                    'transitfit'],
                  python_requires='>=3.6',
                  include_package_data=True )
