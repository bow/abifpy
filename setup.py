#!/usr/bin/env python

from setuptools import setup, find_packages

from abifpy import __version__


version = __version__
long_description = open("README").read()

setup(
    name = "abifpy",
    version = version,
    description = "abifpy is a module for reading ABI Sanger sequencing trace files.",
    long_description = long_description,
    author = "Wibowo Arindrarto",
    author_email = "bow@bow.web.id",
    py_modules = ['abifpy'],
    url = "http://github.com/bow/abifpy/",
    license = "MIT",
    zip_safe = False,
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
