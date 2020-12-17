#!/usr/bin/env python3

import os
import pathlib
from setuptools import setup, find_packages

here = pathlib.Path(__file__).parent.resolve()

directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
  long_description = f.read()

setup(name='target',
      version='0.5.9',
      description='For use with TransitFit',
      author='Steven Charles-Mindoza',
      license='MIT',
      long_description=long_description,
      long_description_content_type='text/markdown',
      packages=find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
      ],
      install_requires=['numpy',
                        'pandas', 
                        'lightkurve',
                        'datetime',
                        'shutil',
                        'multiprocessing',
                        'functools',
                        'transitfit'],
      python_requires='>=3.6',
      include_package_data=True)
