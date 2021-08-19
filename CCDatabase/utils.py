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
import re
import cachetools
import pandas as pd


caches = {"npz": cachetools.LRUCache(maxsize=8),
          "json": cachetools.LRUCache(maxsize=10),
          "non_fdet_terms": cachetools.LRUCache(maxsize=1),
          "non_fdet_int": cachetools.LRUCache(maxsize=1),
          "fdet_terms": cachetools.LRUCache(maxsize=1),
          "fdet_grouped": cachetools.LRUCache(maxsize=1),
          "find_A": cachetools.LRUCache(maxsize=4),
          "find_3": cachetools.LRUCache(maxsize=4),
          "elst_ref": cachetools.LRUCache(maxsize=1),
          "elst_sum_iso": cachetools.LRUCache(maxsize=1),
          "elst_change_ref": cachetools.LRUCache(maxsize=1)}


def clear_caches(cachedict):
    for v in cachedict.values():
        v.clear()

def russdoll(obj):
    """
    Note
    ----
    not used for the moment
    
    Parameters
    ----------
    obj: obj
        the nested object
    
    Returns
    -------
    obj 
        the first object in the nested structure whose length is not 1
    """
    while True:
        try:
            if len(obj) == 1:
                n = obj[0]
            if n == obj:
                break
            else:
                obj = n
        except:
            break
    return obj

def str2arr(s, np_ = True):
    """
    Note
    ----
    Takes a list/array(vectors, matrices) written as a string by pandas, and turns it back into a list/array.
    If given a non-str will return it as it is.
    
    Parameters
    ----------
    s: str, any
        the string  or object to process
    np_: bool
        True if np.array is desired, False if list is desired
    
    Returns
    -------
    np.ndarr/list or s
        the array/list if possible
    """
    if type(s) != str or not re.match(r"\[?(?:\[(?: ?[0-9]+,? ?)+\](?:,|\n)? ?)+\]?",s):
        return s
    else:
        s = re.sub(r"(\d) ", r"\1, ", s)
        s = re.sub(r"\n ", r", ", s)
        r = js.loads(s)
        if np_:
            r = np.array(r)
        return r
            
def vec_df(df,np_):
    """
    Note
    ----
    Turns every array/list(vectors, matrices) stored as text into array/list
    
    Parameters
    ----------
    df: pd.DataFrame
        the dataframe to process
    np_: bool
        True if np.ndarray is desired, False if list is desired
    
    Does
    ----
    changes values within df
    """
    lmbd = lambda x: str2arr(x,np_=np_)
    for c in df.columns:
        df[c] = df[c].apply(lmbd)
        
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

def path_basename(path):
    """
    Note
    ----
    same as os.path.basename(path) but avoids issues with trailing slash.
    
    Parameters
    ----------
    path: str
        the path to get the basename of
    
    Returns
    -------
    list
        the basename
    """
    basename = os.path.basename(path)
    if basename == "":
        basename = os.path.basename(os.path.split(path)[0])
    return basename

@cachetools.cached(cache=caches["json"])
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
    with open(fname, "r") as f:
        try:
            return js.load(f)
        except:
            return {}  

def dump_js(obj, fname):
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
    with open(fname, "w") as f:
        js.dump(obj, f)

def dict_from_file(filepath):
    """
    Note
    ----
    Returns whatever is in the file in the form of a dictionary.
    For dataframe, strips the NaN to allow counting (rq_in_keys)
    If the file does not exist, returns empty dictionary
    
    Parameters
    ----------
    filepath: str
        the filepath
    
    Returns
    -------
    dict
        the dictionary of what is contained in the file
    """
    if not os.path.exists(filepath):
        return {} 
    type_ = filepath.split(".")[-1]
    if type_ == "json":
        return load_js(filepath)
    elif type_ == "xlsx":  # excel file
        df = pd.read_excel(filepath)
        return {k: v.dropna() for k,v in dict(df).items()}
    elif type_ == "csv":  # csv file
        df = pd.read_csv(filepath)
        return {k: v.dropna() for k,v in dict(df).items()}
    else:
        raise TypeError("This type of non-json file is not implemented yet")
                    
def deal_with_type(obj, condition=False, to=None):
    """
    Note
    ----
    Changes type to handle different inputs
    
    Parameters
    ----------
    obj: any
        the object whose type you want to change
    condition: list of types or type
        only change if of a type in the list. You can give both None and NoneType for NoneType
    to: type/func
        type to convert to, function to return. e.g. list, os.getcwd
    
    Returns
    -------
    obj
        the object with the desired type (or processed as desired)
    """
    if condition != False:
        condition = type(None) if condition is None else condition   # None => NoneType
        condition = [condition] if type(condition)==type else condition  # if only one type, make it list of types
        if type(obj) in condition:
            try:
                return to(obj)
            except:
                print("utils.deal_with_type: failed to return to(obj), returning obj")
                return obj
        else:
            return obj
    else:
        return to(obj)
    
def rq_in_keys(d, q, nvals=1):
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
        if len(d[q]) < nvals:
            return False
        else:
            return True
        
def cq_in_keys(d, q, atomstring=False, states=False, state=False, exc=False):
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
        Whether the complex quantity is present in d.
        If state, if it's present for that state.
        If states, if it's present up to that state.

    """
    if states and state:
        print("You used both states and state. states trumps state.")
        state = False
    shift = 1 if exc else 0
    if not states and q not in d.keys():
        return False
    elif states:
        if q in d.keys():
            if type(d[q]) in [list, tuple, dict]:  # a data container 
                if atomstring:
                    if len(d[q]) < states - shift:
                        return False
                    else:
                        for s in range(shift, states+1):
                            if np.array([atomkey not in d[q][s].keys() for atomkey in atomstring.split(",")]).any():
                                return False
                        return True
                        
                else:
                    return len(d[q]) >= states - shift
            else:  # float, str, etc...
                return False
        else:
            qs = ["{}_{}".format(q, s) for s in range(shift, states+1)]
            if atomstring:
                if np.array([x not in d.keys() for x in qs]).any():
                    return False
                for qn in qs:
                    if np.array([atomkey not in d[qn].keys() for atomkey in atomstring.split(",")]).any():
                        return False
                return True
            else:
                return np.array([x in d.keys() for x in qs]).all()
    elif state:
        q = "{}_{}".format(q, state)
        if atomstring:
            pass
        else:
            return q in d.keys()
        
def setupLogger(to_console=True, to_log=False, logname="CCDatabase.log", printlevel=10):  
    """
    Note
    ----
    Never add handlers directly, particularly within functions that get called several times, as this will result in handler multiplication.
    Use this function which prevents it and ensure 1 handler per type is present at most.
    
    Parameters
    ----------
    to_console: bool
        whether to print logging to the console
    to_log: bool
        whether to print logging to file
    logname: str
        the log filename. Default is CCDatabase.log
    printlevel: int
        logging level: 10,20,30,40,50. Accessible also as logging.DEBUG/INFO/WARNING/ERROR/CRITICAL (respectively)
        
    Does
    ----
    Sets up the logger avoiding handlers repetition
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
    