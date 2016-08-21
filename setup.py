#! /usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv, exit
import os, re

if "upl" in argv[1:]:
    os.system("python setup.py register -r pypi")
    os.system("python setup.py sdist upload -r pypi")
    exit()

m = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "MCres", "_version.py")).read()
version = re.findall(r"__version__ *= *\"(.*?)\"", m)[0]

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name = 'MCres',
    packages = ['MCres'],
    version = version,
    description = 'Provides simple filtering and visualization tools for Monte-Carlo simulation data',
    long_description = open("README.rst").read() + "\n\n"
                    + "Changelog\n"
                    + "---------\n\n"
                    + open("HISTORY.rst").read(),
    license = "GNU General Public License v3 or later (GPLv3+)",
    author = 'Guillaume Schworer',
    author_email = 'guillaume.schworer@obspm.fr',
    url = 'https://github.com/ceyzeriat/MCres/',
    download_url = 'https://github.com/ceyzeriat/MCres/tree/master/dist',
    keywords = ['emcee','processing','data','big data','corner','cornerplot','triangle','plot'],
    package_data={"": ["README.rst", "LICENSE", "HISTORY.rst"]},
    include_package_data=True,
    install_requires=["corner"],
    classifiers = [
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Documentation :: Sphinx',
            'Topic :: Scientific/Engineering :: Physics']
)
