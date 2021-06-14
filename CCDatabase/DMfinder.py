#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:58:40 2021

@author: nico
"""
import CCDatabase.utils as ut
import os
import glob as gl


def find_matrix(fname, wildcardlist=[], jsonfile="", prop_key="ext_DM"):
    path = os.path.split(fname)[0]
    if not wildcardlist:
        wildcardlist = ["",""]  # TODO put real options
    elif type(wildcardlist) not in [tuple, list]:
        wildcardlist = [wildcardlist]
    if not jsonfile:
        suspects = ["CCParser.json", "qcep.json", "parser.json", "parsed.json"]
        present = [ i for i in suspects if os.path.exists(os.path.join(path, i))]
        if len(present) == 0:
            raise FileNotFoundError("could not find any json parser file")
        elif len(present) != 1:
            raise FileNotFoundError("Too many json parser file, specify one")
        else:
            jsonfile = present[0]
    json_filepath = jsonfile if os.path.split(jsonfile)[0] else os.path.join(path,jsonfile)
    old = ut.load_js(json_filepath) if os.path.exists(json_filepath) else {}
    matching = [os.path.split(fp)[1] for wc in wildcardlist for fp in gl.glob(os.path.join(path,wc))]
    old[prop_key] = matching  # always overwrites values because if file is no longer there, pointer is useless
    ut.dump_js(old, json_filepath)
    return matching
    
    
