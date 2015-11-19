#!/usr/bin/env python

from distutils.core import setup

setup(name='sdmreader',
      version='1.2',
      description='Read/parse NRAO Science Data Model files',
      author='Casey Law',
      author_email='caseyjlaw@gmail.com',
      url='http://github.com/caseyjlaw/sdmreader',
      packages=['sdmreader'],
      requires=['pwkit'],
     )
