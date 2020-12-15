from setuptools import setup
setup(name='retriever',
version='1.0.0',
description='For use with TransitFit',
url='TBA',
author='Steven Charles-Mindoza',
author_email='stevemindoza@gmail.com',
license='Open GL 3.0',
packages=['retriever'],
python_requires='>=3.6',
install_requires=['batman-package', 'dynesty', 'numpy', 'matplotlib',
                      'pandas', 'ldtk', 'lightkurve', 'transitfit'],
zip_safe=False)
