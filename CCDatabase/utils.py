#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:21:01 2020

@author: nico
"""
import json as js
import numpy as np
import os
import logging

def mkdif(a):
    """
    Note
    ----
    makes directory if it does not exist yet

    Parameters
    ----------
    a: str
        directory name
    """
    if not os.path.exists(a):
        os.makedirs(a)

def split_path(path):
    """
    Note
    ----
    same as os.path.split(path) but avoids issues with trailing slash.
    
    Parameters
    ----------
    path: str
        the path to split
    
    Returns
    -------
    list
        the split list
    """
    splt = os.path.split(path)
    if splt[-1] == "":
        splt = os.path.split(splt[0])
    return splt

def load_js(fname):
    """
    Note
    ----
    Some datatype may change when load=>dump=>load (e.g. tuple=> list, {1:"one"}=>{'1':'one'})
    
    Parameters
    ----------
    fname: str
        filepath to load
    
    Returns
    -------
    dict
        the dictionary in the json file
    """
    with open(fname,"r") as f:
        return js.load(f)

def dump_js(obj,fname):
    """
    Note
    ----
    Some datatype may change when load=>dump=>load (e.g. tuple=> list, {1:"one"}=>{'1':'one'})
    
    Parameters
    ----------
    obj: dict
        the dictionary to dump
    fname: str
        filepath to dump in
    
    Does
    ----
    dumps obj in fname
    """
    with open(fname,"w") as f:
        js.dump(obj,f)
        
def deal_with_type(obj,condition=False,to=None):
    """
    Note
    ----
    Changes type to handle different inputs
    
    Parameters
    ----------
    obj: any
        the object whose type you want to change
    condition: type
        only change if of "type". Give None for NoneType
    to: type/func
        type to convert to, function to return. e.g. list, os.getcwd
    
    Returns
    -------
    obj
        the object with the desired type (or processed as desired)
    """
    if condition != False:
        if type(obj) == condition:
            try:
                return to(obj)
            except:
                print("utils.deal_with_type: failed to return to(obj), returning obj")
                return obj
        elif condition == None and obj == None:
            try:
                return to()
            except:
                print("utils.deal_with_type: failed to return to(), returning obj")
                return obj
        else:
            return obj
    else:
        return to(obj)
    
def rq_in_keys(d,q, nvals=1):
    """
    Parameters
    ----------
    d : dict
        the dictionary to check
    q : str
        the quantity to check
    nvals : int, optional
        How many items should be available for d[q]. The default is 1.

    Returns
    -------
    bool
        Whethere the raw quantity is in d the desired number of times

    """
    if q not in d.keys():
        return False
    else:
        if len(d[q])<nvals:
            return False
        else:
            return True
        
def cq_in_keys(d,q, states=False, state=False):
    """
    Parameters
    ----------
    d : dict
        the dictionary to check
    q : str
        the quantity to check
    states : bool/int
        Whether it is a simple value or several. 

    Returns
    -------
    bool
        DESCRIPTION.

    """
    if states and state:
        print("You used both states and state. states trumps state.")
        state=False
    if not states and q not in d.keys():
        return False
    elif states:
        if q in d.keys():
            if type(d[q]) in [list, tuple, dict]:  # a data container 
                return len(d[q]) >= states
            else:  # float, str, etc...
                return False
        else:
            qs = ["{}_{}".format(q,s) for s in range(states)]
            return np.array([x in d.keys() for x in qs]).all()
    elif state:
        q = "{}_{}".format(q,state)
        return q in d.keys()
        
def setupLogger(to_console=True, to_log=False, logname="CCDatabase.log", printlevel=10):  
    """
    """
    if "ccdlog" not in globals(): 
        global ccdlog
        ccdlog = logging.getLogger("ccd")
    else:
        ccdlog = globals()["ccdlog"]
    handlernames = [i._name for i in ccdlog.handlers]
    handlerdict = {}
    if False in [type(to_console) == bool, type(to_log) == bool]:
        ccdlog.critical("\"to_console\" and \"to_log\" should be boolean!!")
        if type(to_log) == str and "." in to_log and logname == "CCDatabase.log":
            ccdlog.critical("""It appears you might have given your desired logname as \"to_log\" instead of changing \"logname\".+
                            I am using  {} as logname and setting to_log to True""".format(to_log))
            to_log, logname = True, to_log
    if to_console:
        handlerdict["console"] = logging.StreamHandler()
    if to_log:
        handlerdict["file"] = logging.FileHandler(logname)
    if handlerdict:  # there is at least a handler
        for k,handler in handlerdict.items():
            if k not in handlernames:
                handler._name = k
                handler.setLevel(printlevel)
                handler.setFormatter(logging.Formatter('%(levelname)s - %(name)s.%(funcName)s:  %(message)s'))
                ccdlog.addHandler(handler)
    else: # there is no handler
        ccdlog.setLevel(50)
    