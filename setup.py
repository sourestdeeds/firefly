from setuptools import setup
setup(name='retriever',
version='1.0.0',
description='For use with TransitFit',
url='https://github.com/sourestdeeds/retrieval',
author='Steven Charles-Mindoza',
author_email='stevemindoza@gmail.com',
license='GPL-3.0 Licence',
packages=['retriever'],
python_requires='>=3.6',
install_requires=['batman-package', 'dynesty', 'numpy', 'matplotlib',
                      'pandas', 'ldtk', 'lightkurve', 'transitfit'],
zip_safe=False)
