#!/usr/bin/env python3

import os
import pathlib
from setuptools import setup

here = pathlib.Path(__file__).parent.resolve()

directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
  long_description = f.read()

setup(name='retrieval',
      version='1.0.0',
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
      install_requires=['batman-package', 'dynesty', 'numpy', 'matplotlib',
                      'pandas', 'ldtk', 'lightkurve', 'transitfit'],
      python_requires='>=3.6',
      include_package_data=True)
