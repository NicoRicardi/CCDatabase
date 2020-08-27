#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:15:05 2020

@author: nico
"""
import os
import CCDatabase.utils as ut

def raw_to_complex(path=None, rawfile="CCParser.json", raw_key="", n=0, linenumbers=True):
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
#    n -= 1  # parameter uses human counting, code uses python counting
    path = ut.deal_with_type(path,condition=None,to=os.getcwd)
    if raw_key == "":
        raise ValueError("you didn't specify the raw quantity's key")
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile)
    to_return = raws[raw_key][n][0] if linenumbers else raws[raw_key][n]
    return to_return

"""
for quantity functions it is often useful to use some general function with several parameters,
define all of them but "path" and "n" into a lambda function, which is then saved as quantity function.
e.g. 
"ex_en": lambda path,n: raw_to_complex(path=path,n=n,rawfile="CCParser.json",raw_key="exc_energies_rel")
"""
ccp_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path,n=n,rawfile="CCParser.json",raw_key="exc_energy_rel")}

qcep_ccp_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path,n=n,rawfile="CCParser.json",raw_key="exc_energy")}

qcep_funcs = {
        "ex_en": lambda path,n: raw_to_complex(path=path,n=n,rawfile="qcep.json",raw_key="exc_energy")}