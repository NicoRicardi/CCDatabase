#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:15:05 2020

@author: nico
"""
import os
import CCDatabase.utils as ut
import numpy as np
import copy as cp

def raw_to_complex(path=None, rawfile="CCParser.json", raw_key="", n=0, linenumbers=True, first_only=False, group_values=False):
    """
    Note
    ----
    Simply copies a raw quantity as a complex quantity. Handles change of key, picks required "state" (n-th value of raw_quantity)
    
    Parameters
    ----------
    path: str
        path in which to look for raw quantity
    rawfile: str
        the json file with raw quantities. default is ccp
    raw_key: str
        the key of the raw quantity in the json dictionary
    n: int
        which value to take. Uses fortran/human counting, default is 1
    linenumbers: bool
        Whethere the quantity in the json file will have linenumbers. Default is True.
        e.g. json_dict[raw_key]=[[val1,line1],[val2,line2], ...] => use True
             json_dict[raw_key]=[[val1,val2,val3, ...] => use False
     Returns
     -------
     dictionary item
         generally a float, potentially str.
    """
    if raw_key == "":
        raise ValueError("you didn't specify the raw quantity's key")
    if n < 0:
        print("CRITICAL - raw_to_complex: negative index. Most likely you didn't specify that a quantity is in ex_qs")
        raise ValueError  #message in print to appear also when in "try"  statements
    if first_only and n > 0:
        raise ValueError("You want a value only, but are asking for n > 0")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile)
    vals = group_values(raws[raw_key]) if group_values else raws[raw_key]
    to_return = vals[n][0] if linenumbers else vals[n]
    return to_return

def group_values(vals):
    """
    Note
    ----
    Under construction
    """
    decimals = lambda x: len(str(x).split(".")[-1])
    options = {n: [] for n in range(len(vals))}
    decs = {n: decimals(vals[n]) for n in range(len(vals))}
    max_d = max(decs.values())
    if min(decs.values()) == max_d:
        return vals
    for n1,v1 in enumerate(vals):
        for n,v2  in enumerate(vals[n1+1:]):  # avoid double iteration (e.g. (0,1),(1,0))
            n2 = n + n1 + 1 # so that v2 = vals[n2]
            if decs[n1] > decs[n2]:
                if round(v1, decs[n2]) == v2:
                    options[n2].append(v1)
            elif decs[n1] < decs[n2]:
                if round(v2, decs[n1]) == v1:
                    options[n1].append(v2)
        if len(set(options[n1])) == 1:
            vals[n1] = options[n1][0]
        elif len(options[n1]) > 1:  # if len(options) == 0 continue
            dlist = [decimals(i) for i in options[n1]]
            dd = {d: dlist.count(d) for d in range(1,max(dlist)+1)}
            to_use = False
            for d in dd.keys():
                rounded = [round(o,d) for o in options[n1] if decimals(o) >= d]
                if len(set(rounded)) != 1:  # all equal if rounded with d decimals
                    break
                to_use = cp.copy(rounded[0])
            vals[n1] = to_use
    return vals

"""
for quantity functions it is often useful to use some general function with several parameters,
define all of them but "path" and "n" into a lambda function, which is then saved as quantity function.
e.g. 
"ex_en": lambda path, n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="exc_energies_rel")
"""
ccp_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="exc_energy_rel"),
        "osc": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="osc_str"),
        "SCF": lambda path,n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="scf_energy", first_only=True)}

qcep_ccp_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="exc_energy"),
        "osc": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="osc_strength"),
        "SCF": lambda path,n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="scf_energy", first_only=True)}

qcep_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="qcep.json", raw_key="exc_energy"),
        "osc": lambda path,n: raw_to_complex(path=path, n=n-1, rawfile="qcep.json", raw_key="osc_strength"),
        "SCF": lambda path,n: raw_to_complex(path=path, n=n, rawfile="qcep.json", raw_key="scf_energy", first_only=True)}