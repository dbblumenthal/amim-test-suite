"""
Setup module adapted from setuptools code. See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
	name='NetCore',
	version='0.1',
	description='NetCore: A Network Propagation approach using node coreness',
	url='https://github.molgen.mpg.de/barel/NetCore',
	author='Gal Barel',
	author_email='barel@molgen.mpg.de',
	license='MIT',
	classifiers=[
	  'Development Status :: 2 - Pre-Alpha',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Programming Language :: Python :: 3.6',
      'Programming Language :: Python :: 3.7'
	],
	packages=find_packages(exclude=['copy', 'itertools', 'os', 'sys', 'math', 'random','pickle']),
	install_requires=[
        'argparse>=1.1',
        'networkx>=2.3',
        'numpy>=1.16.3',
        'matplotlib>=3.1.1',
        'pandas>=0.24.2',
        'scipy>=1.2.1',
        'seaborn>=0.9.0']
)