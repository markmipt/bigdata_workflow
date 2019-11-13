#!/usr/bin/env python

'''
setup.py file for Bigdata_workflow
'''
from setuptools import setup
version = open('VERSION').readline().strip()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name             = 'Bigdata_workflow',
    version          = version,
    description      = '''Proteomics search wrapper''',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    author           = 'Mark Ivanov & Lev Levitsky',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'https://github.com/markmipt/bigdata_workflow',
    packages         = ['bigdata_workflow', ],
    install_requires = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 3',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    license          = 'License :: OSI Approved :: Apache Software License',
    entry_points     = {'console_scripts': ['bigdata_workflow_convert = bigdata_workflow.convert:run',
                                            # 'scav2diffacto = scavager.scav2diffacto:run',
                                            # 'scav2nsaf = scavager.scav2nsaf:run'
                                            ]}
    )