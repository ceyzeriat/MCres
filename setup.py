#!/usr/bin/env python
from sys import argv, exit

if "upl" in argv[1:]:
    import os
    os.system("python setup.py register -r pypi")
    os.system("python setup.py sdist upload -r pypi")
    exit()


from distutils.core import setup
setup(
    name = 'MCres',
    packages = ['MCres'],
    version = '0.4beta',
    description = 'Provides simple filtering and visualization tools for Monte-Carlo simulation data',
    author = 'Guillaume Schworer',
    author_email = 'guillaume.schworer@obspm.fr',
    url = 'https://github.com/ceyzeriat/MCres/',
    download_url = 'https://github.com/ceyzeriat/MCres/tree/master/dist',
    keywords = ['emcee', 'processing','data','big data',''], # arbitrary keywords
    classifiers = ['Development Status :: 4 - Beta','Environment :: Console','Intended Audience :: Education','Intended Audience :: Science/Research','License :: OSI Approved :: MIT License','Operating System :: Unix','Programming Language :: Python :: 2.7','Topic :: Documentation :: Sphinx', 'Topic :: Scientific/Engineering :: Physics'],
)

# http://peterdowns.com/posts/first-time-with-pypi.html
