#!/usr/bin/env python

from setuptools import setup
import glob

scripts = glob.glob("*.p*")

version = "0.0.0"

setup(
    name='PatchPolish',
    version=version,
    description='Polish assembly patches from raw ONT reads.',
    author='Michael Alonge',
    author_email='malonge11@gmail.com',
    packages=['patchpolish_utilities'],
    package_dir={'patchpolish_utilities': 'patchpolish_utilities/'},
    install_requires=[
              'pysam'
          ],
    scripts=scripts,
    zip_safe=True
)