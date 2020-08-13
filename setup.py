#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:28:40 2019

@author: nico
"""
from distutils.core import setup
setup(
    name='CCDatabase',
    version='0.0.1a',
    description='Neat tools to manage a computational chemistry database',
    author='Niccolo Ricardi, Cristina Gonzalez',
    author_email='Niccolo.Ricardi@unige.ch, Cristina.Gonzalez@unige.ch',
    package_dir={'CCDatabase': 'CCDatabase'},
    packages=['CCDatabase'],
)
