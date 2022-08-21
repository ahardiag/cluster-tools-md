#!/usr/bin/env python

from distutils.core import setup

requirements = ['matplotlib',
                'MDAnalysis>2',
                'scikit-learn',
                'scipy',
                'pandas',
                'numba'
                ]

setup(name='cluster_tools_md',
      version='1.0.0',
      description='A command-line interface for analysis routines of Molecular Dynamics data with aggregating molecules.',
      author='Arthur Hardiagon',
      author_email='ahardiag@gmail.com',
      url='https://github.com/ahardiag/cluster_tools_md',
      packages=['src','src.multirdf'],
      install_requires = requirements,
      entry_points = {
              'console_scripts': [
                  'clustsize    = src.cli_clustsize:main',  
                  'eccentricity = src.cli_eccentricity:main',
                  'multirdf     = src.multirdf.cli_multirdf:main',
                  'radialdens   = src.cli_radialdens:main',
                  'plotclustersize =  src.cli_plotclustersize:main'
              ], 

          },
     )