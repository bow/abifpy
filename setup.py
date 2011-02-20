#!/usr/bin/env python

import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "abifpy",
    version = "0.1",
    author = "Wibowo Arindrarto",
    author_email = "w.arindrarto@gmail.com"
    description = ("Python module to read .ab1 trace files."),
    license = "MIT",
    keywords = "ab1 trace sequencing"
    url = "http://github.com/warindrarto/abifpy",
    packages = ['abifpy'],
    long_description = read('README.rst'),
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
