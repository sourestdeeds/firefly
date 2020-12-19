#!/usr/bin/env python3

import os
import pathlib
import setuptools

here = pathlib.Path(__file__).parent.resolve()

directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup( name='firefly',
                  version=get_version('firefly/__init__.py'),
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
                                    'transitfit'],
                  python_requires='>=3.6',
                  include_package_data=True )
