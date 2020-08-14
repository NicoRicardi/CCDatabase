#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some notes:
    1. For the time being, in looping functions (raw_quantities, complex_quantities), 
    ext(i.e. the wildcard for your output files, e.g. "*.out") and ignore (i.e. wildcard for files to ignore, e.g. "slurm*")
    need to be constant. Either launch these functions for different types of outputs,or implement different outputs.

"""

import json as js
import itertools as ittl
import datetime as dt
import CCParser as ccp
import glob as gl
import numpy as np
import pandas as pd
import os
import re

### Simple utilities
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

### Database cleaning

def move_to_trash(expr, startdir="", trash="", keep=True, logfile="operations.log", histfile="history.json"):
    """
    Parameters
    ----------
    expr: list/str
        wildcard expression (e.g. *.out, Densmat*.txt...) of files to keep
    startdir; str
        starting directory. default is current working directory
    trash: str
        trash location. Default is "./Trash", creates if necessary
    keep: bool
        whether to keep the matching ones (and move the others), or to move the matching ones(and keep the others).
        default is the former
    logfile: str
        operation logfile to add entry to
    histfile: str
        history file, a more detail file which allows to undo a move
    Does
    -------
    Either keeps only matching files or moves only matching files to "Trash" directory.
    This is added as an entry in the logfile (date, moved to trash files [not] matching {expr}).
    The exact move operations are printed to screen and written in the history json file.
    Files are never overwritten, and their modified names (e.g. file_2.ext) are stored, which allows to undo moving operations.
    """
    ### Handling type of expr
    if type(expr) == str:
        expr = [expr]
    if type(expr) != list:
        raise TypeError("Expr must be either list/tuple or string")
    ### Turning wildcards into regex
    expr = [i.replace("*",".+").replace("?",".") for i in expr]    
    ### Handling default startdir if none specified
    startdir = os.getcwd() if not startdir else startdir
    ### Handling default trash if none specified
    trash = os.path.join(os.getcwd(),"Trash") if not trash else trash
    mkdif(trash)
    try:
        histdict = load_js(histfile)  # read history file and look for highest ID
        ID = max([int(i) for i in histdict.keys()]) + 1  # get ID and add 1
    except FileNotFoundError:  #history file is not there yet
        histdict = {}
        ID = 1
    moved = []  # will be list of files moved
    for dirname,subdirlist,filelist in os.walk(startdir):   # walking through tree
        for fname in filelist:  # loop over files found
            for e in expr:  # loop over wildcard expressions
                if not(keep^bool(re.match(e,fname))):  #move or not based on "keep"
                    src = os.path.join(dirname,fname)
                    dest = os.path.join(trash,fname)
                    if os.path.exists(dest):
                        splt = dest.split(".")
                        bname = ".".join(splt[:-1])
                        ext = splt[-1]
                        count = 2
                        while os.path.exists(dest):  # avoid overwriting it
                            dest = "{}_{}.{}".format(bname,count,ext)
                            count += 1
                    os.move(src,dest)
                    print("moved {} to {}".format(os.path.join(dirname,fname), trash))
                    moved.append((os.path.abspath(src),os.path.abspath(dest)))
    histdict.update({ID:moved})  # {..., ID:[(src1,dest1),(src2,dest2),...]}      
    dump_js(histdict,histfile)          
    with open(logfile,"a") as log:
        neg = "not" if keep else ""
        date = dt.date.today().strftime("%d-%b-%Y")
        exprstr = ",".join(expr)
        log.write("""{}, ID: {}, moved to trash files {} matching with "{}"\n""".format(date,ID,neg,exprstr))

def keep_only(expr, startdir="", trash="", logfile="operations.log", histfile="history.json"):
    """
    Parameters
    ----------
    expr: list/str
        wildcard expression (e.g. *.out, Densmat*.txt...) of files to keep
    startdir; str
        starting directory. default is current working directory
    trash: str
        trash location. Default is "./Trash", creates if necessary
    Does
    -------
    Keeps only matching files and moves the others to "Trash" directory
    Writes to logfile and histfile see "move_to_trash"
    """
    move_to_trash(expr, startdir=startdir, trash=trash, keep=True, logfile=logfile, histfile=histfile)
    
def remove_matching(expr, startdir="", trash="", logfile="operations.log", histfile="history.json"):
    """
    Parameters
    ----------
    expr: list/str
        wildcard expression (e.g. *.out, Densmat*.txt...) of files to keep
    startdir; str
        starting directory. default is current working directory
    trash: str
        trash location. Default is "./Trash", creates if necessary
    Does
    -------
    Moves only matching files and keeps the others to "Trash" directory.
    Writes to logfile and histfile see "move_to_trash"
    """
    move_to_trash(expr, startdir=startdir, trash=trash, keep=False, logfile=logfile, histfile=histfile)

def undo_move(IDs, logfile="operations.log", histfile="history.json"):
    """
    Parameters
    ----------
    IDs : list/int/str
        the ID(s) of the move operation one wishes to undo. This is reported in the logfile
    histfile : str
        the history file which contains all filepaths to restore the files

    Does
    -------
    Restores all files in their original locations. 
    If a file with the same filepath has been added after the deletion, 
    overwriting is avoided and the file is copied as file_n.ext.
    """
    if type(IDs) == str:  # if only 1 ID
        IDs = [IDs]
    histdict = load_js(histfile)
    with open(logfile, "r") as f:
        lines = f.readlines()
    for ID in IDs:  # looping over all IDs
        to_restore = histdict[str(ID)]  # (str because from json) list of items to restore [[original_path1,trash_path1],[original_path2,trash_path2],...]
        for file in to_restore:
            dest = file[0]  # where we would like to restore it
            src = file[1]  # where it is in the trash
            splt = dest.split(".")  
            bname = ".".join(splt[:-1]) #***checkfix
            ext = splt[-1]
            count = 2
            while os.path.exists(dest):  # avoid overwriting it
                dest = "{}_{}.{}".format(bname,count,ext)
                count += 1
            os.move(src,dest)
            print("restored {}".format(dest))
        del histdict[str(ID)]  # removing entry from histfile
        for n,line in enumerate(lines):
            if "ID: {}".format(ID) in line:
                del lines[n]  # removing entry from logfile
    dump_js(histdict,histfile)  # updating histfile

    with open(logfile, "w") as f:  # updating logfile
        for line in lines:
                f.write(line)   

### Quantities-related
                
def q_in_keys(d,q, times=1):
    """
    Parameters
    ----------
    d : dict
        the dictionary to check
    q : str
        the quantity to check
    times : int, optional
        How many items should be available for d[q]. The default is 1.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    if q not in d.keys():
        return False
    else:
        if len(d[q])<times:
            return False
        else:
            return True
        
def find_and_parse(path, ext="*.out", ignore="slurm*", parser=None, parser_args=None, parser_kwargs=None): 
    """
    Parameters
    ----------
    path: str
        the folder to work in
    ext: str
        wildcard-like expression to match. Default is "*.out"
    ignore: str
        wildcard-like expression of files to ignore. Default "slurm*"
    parser: None/function
        parser to call. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs as parser_kwargs
        Your parser should take care of writing the data to file
    parser_args: list/tuple
        provide a list of the positional arguments for your parser
    parser_kwargs: dict
        dictionary of {kw1: arg1, kw2: arg2,...}
    Returns
    -------
    What your parser returns
        the parser container
    """
    path = os.path.abspath(path)
    files = [i for i in gl.glob(os.path.join(path,ext)) if i not in gl.glob(os.path.join(path,ignore))]  # e.g. all *.out which are not slurm*.out
    if len(files) == 1:
        file = files[0]
    else:
        lens = [len(f) for f in files]
        idx = lens.index(min(lens))        
        file = files[idx]
        print("Several matching files!! Using the shortest file ({})".format(file))
    ### Not sure if this next line can be broken
    parser_args = list(parser_args) if type(parser_args)==tuple else parser_args  
    parser, parser_args, parser_kwargs = (ccp.Parser,[file],dict(software="qchem",to_file=True, to_console=False, to_json=True))  if parser == None else (parser, [file]+parser_args, parser_kwargs)
    data = parser(*parser_args,**parser_kwargs)
    return data

def get_joblist(fname="variables.json"):
    """  

    Parameters
    ----------
    fname : str
        json file with the joblist (or specifications) to analyse

    Returns
    -------
    list(tuples)
        the list of "jobs" (e.g. [job1,job2,job3...] with job being ("A","B","basis",[...],"calc")) 
    """
    data = load_js(fname) 
    if "joblist" in data.keys():
        return(data["joblist"])
    else:
        ignore = data["ignore"] if "ignore" in data.keys() else []
        return [i for i in  ittl.product(*[data[i] for i in data["levels"]]) if i not in ignore]

def raw_quantities(path=None, qlist="variables.json", ext="*.out", ignore="slurm*", parser_file="CCParser.json", parser=None, parser_args=None, parser_kwargs=None):
    """
    Parameters
    ----------
    path: str
        the path to check 
    qlist : list/str
        raw quantities to check. Default checks from variables.json
        If there is no ",", searches for the quantity directly in the parser's json file in path.
        Otherwise, the string before the comma is considered as the folder name to look into.
        If this folder name is available as a subfolder (e.g. F&T, Macrocycles, Prepol)
    ext: str
        wildcard-like expression to match. Default is "*.out"
    ignore: str
        wildcard-like expression of files to ignore. Default "slurm*"
    parser_file: str
        name of the parser's json file. Default is "CCParser.json"
    parser: dict{q: func}, func, None
        parser to call (for each q). If not dictionary, turns into dict. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
        NB. If the file is human-generated and cannot be obtained from a parser, pass False
    parser_args: dict[parser: list/tuple]
        provide a list of the positional arguments for your parser, for each parser.
        If one parser only you can pass the list alone, it is turned into a dict
    parser_kwargs: dict[parser: dict]
        provide a dict of the keyword arguments for your parser, for each parser.
        If one parser only you can pass the dict alone, it is turned into dict[parser: dict]
    Does
    ----
    checks if the quantities are available in all "jobs". If necessary, it parses.
    Then prints what is missing
    
    Returns
    -------
    list
        the quantities missing in path
    """        
    ### Processing user input
    # qlist
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            qlist = load_js(qlist)["raw_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    # path
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    ### parser, parser_args, parser_kwargs
    # parser
    if type(parser) == list:
        if len(parser) == len(qlist):
            parser = {qlist[n]: parser[n] for n in range(len(qlist))}
        else:
            raise ValueError("""You gave a list as "parser" but the length does not match that of "qlist" """)
    elif type(parser) == dict:
        if np.array([q not in parser.keys() for q in qlist]).any():
            raise ValueError("""No parser for some of your quantities""")
    else:
        parser = {q:parser for q in qlist}
    # parser_args    
    if type(parser_args) == list:
        if len(parser_args) == len(qlist):
            parser_args = {qlist[n]: parser_args[n] for n in range(len(qlist))}
        else:
            raise ValueError("""You gave a list as "parser_args" but the length does not match that of "qlist" """)
    elif type(parser_args) == dict:
        if np.array([q not in parser_args.keys() for q in qlist]).any():
            raise ValueError("""No parser_args for some of your quantities """)
    else:
        parser_args = {q:parser_args for q in qlist}   
    # parser_kwargs    
    if type(parser_kwargs) == list:
        if len(parser_kwargs) == len(qlist):
            parser_kwargs = {qlist[n]: parser_kwargs[n] for n in range(len(qlist))}
        else:
            raise ValueError("""You gave a list as "parser_kwargs" but the length does not match that of "qlist" """)
    elif type(parser_kwargs) == dict:
        if np.array([q not in parser_kwargs.keys() for q in qlist]).any():
            raise ValueError("""No parser for some of your quantities""")
    else:
        parser_kwargs = {q:parser_kwargs for q in qlist}   
    ### Let's get started        
    missing = []
    reparsed = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    data = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    for n,q in enumerate(qlist):
        if type(q)=="str":
            times = 1
        elif type(q) in [list,tuple] and len(q)==2:  # multiple values of a function (e.g. 5 ex. en.)
            times,q = q[0],q[1]
        else:
            raise TypeError("quantity {} not understood".format(q))
        """
        q can be as follows:
            1 "quant": checks for "quant" in the current jobfolder
            2 "subfol,quant": checks for "quant" in subfol
            3 "paralfol,quant": checks for "quant" in paralfol
            4 "non_stdfile,quant": checks for "quant" in "non_stdfile" in current jobfolder
            5 "fol,non_stdfile,quant" or "paralfol,file,quant": loks for "quant" in "non-stdfile" in sub/paral-fol
        """
        stdfile = True  # q is supposed to be in a standard json file
        if "," in q:  # in another calc (not case 1)
            splt = q.split(",")
            if len(splt)==2 and "." in splt[0]:  # case 4
                fname = splt[0]
                type_ = fname.split("."[-1])
                path_tmp = path
                stdfile = True if type_ == "json" else False  # actually not in a standard parser_file
            else:  # cases 2,3,5
                fol = splt[0]
                if len(splt)==3:  # case 5
                    fname = splt[1]
                    type_ = fname.split("."[-1])
                    stdfile = True if type_ == "json" else False  # actually not in a standard parser_file
                else:  # cases 2,3
                    fname = parser_file
                subdirs = gl.glob(os.path.join(path,"*",""))  # subdirectories of path (e.g. MP2_A for F&T)
                paraldirs = gl.glob(os.path.join(split_path(path)[0],"*",""))  # directories in the parent folder (e.g. iso if path="emb")
                if fol in subdirs and fol not in paraldirs:  # it is a subdir (case 2)
                    path_tmp = os.path.join(path,fol)
                elif fol not in subdirs and fol in paraldirs:  # it is a parallel dir (case 3)
                    path_tmp = os.path.join(split_path(path)[0],fol)  
                else:
                    raise ValueError("""the location of your quantity {} cannot be understood.
                             Most likely it is either not present or double""".format(qlist[n]))
            q = splt[-1]
            filepath = os.path.join(path_tmp, fname)
            if not stdfile:
                if type_ == "xlsx":  # excel file
                    df = pd.read_excel(filepath)
                    if q not in df.columns:  # times not available yet
                        missing.append(qlist[n])  # original q, not split
                    continue  # no need to run all other ifs
                elif type_ == "csv":  # csv file
                    df = pd.read_csv(filepath)
                    if q not in df.columns:
                        missing.append(qlist[n])  # original q, not split
                    continue  # no need to run all other ifs
                else:
                    raise TypeError("This type of non-json file is not implemented yet")
        else:  # case 1
            path_tmp = path
            filepath = os.path.join(path_tmp, parser_file)
        ### case is determined. path_tmp and filepath are set
        if path_tmp not in data.keys():
            data[path_tmp] = load_js(filepath) if os.path.exists(filepath) else {}  # load if exists, else empty, will reparse
            if path_tmp not in reparsed.keys():
                reparsed[path_tmp] = {}
            reparsed[path_tmp][parser[qlist[n]]] = False  # full q, not split
        if not q_in_keys(data[path_tmp], q, times=times) and not reparsed[path_tmp][parser[qlist[n]]]:  # quantity not in json, not reparsed yet
            if parser[qlist[n]] != False:
                find_and_parse(path_tmp, ext=ext, ignore=ignore, parser=parser[qlist[n]], parser_args=parser_args[qlist[n]], parser_kwargs=parser_kwargs[qlist[n]])  # CCParser takes care of dumping json
                data[path_tmp] = load_js(filepath)  # let's read reparsed json
            reparsed[path_tmp][parser[qlist[n]]] = True  
        if not q_in_keys(data[path_tmp], q, times=times) and reparsed[path_tmp][parser[qlist[n]]]:
            missing.append(qlist[n])  # original q, not split
            print("{} missing {}".format(path_tmp,q))
    return missing   

def complex_quantities(path=None, qlist="variables.json", reqs=None, ext="*.out", ignore="slurm*", parser=None, parser_args=None, parser_kwargs=None):
    """
    Parameters
    ----------
    path: str
        the path to check 
    qlist : list/str
        coplex quantities to check. Default checks from variables.json
        If json, leave reqs as None (reads file)
        If single quantity(str), give a list as reqs
        If list of quantities, give a dict as reqs
    reqs: None/dict/list
        the requisites for the quantities (see above)
    ext: str
        wildcard-like expression to match. Default is "*.out"
    ignore: str
        wildcard-like expression of files to ignore. Default "slurm*"
    parser: None/function
        parser to call. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
    parser_args: list/tuple
        provide a list of the positional arguments for your parser
        
    Does
    ----
    checks if the quantities are available in all "jobs". If necessary, it parses.
    Then prints which folders do not have.
    
    Returns
    -------
    list
        [quantity1,quantity2, ...] of what is missing
    """   
    ### Processing user input
    # qlist
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            jsdict = load_js(qlist)
            qlist = jsdict["complex_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    else:
        raise TypeError("""qlist can be a list, a json filename to extract it from, or a single quantity(str),..],...}""")
    # reqs
    if type(reqs) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            if "jsdict" not in locals():
                jsdict = load_js(qlist)    
            reqs = jsdict["requisites"]
        else:
            reqs = [reqs]
            print("You only gave one raw quantity as requisite. That can happen but is rare. Sure about it?")
    elif type(reqs) != dict:
        raise TypeError("""reqs should be a dictionary {q1: [req1,req2,..], q2: [req1,req2,..],...} \n
                                                        otherwise a json to get it from (as jsdict["reqs"]) or a single quantity(str)""")
    all_raws = set(ittl.chain.from_iterable(reqs.values()))
    # parser
    if type(parser) == dict:
        if np.array([ i not in parser.keys() for i in all_raws]).any():
            raise ValueError("No parser for some of your requisite raw quantities")
        parserdict = parser
    else:
        parserdict = {q: parser for q in all_raws}
    # parser_args
    if type(parser_args) == dict:
        if np.array([ i not in parser_args.keys() for i in all_raws]).any():
            raise ValueError("No parser_args for some of your requisite raw quantities")
        parserargsdict = parser_args
    else:
        parserargsdict = {q: parser_args for q in all_raws}
    # parser_kwargs
    if type(parser_kwargs) == dict:
        if np.array([ i not in parser_kwargs.keys() for i in all_raws]).any():
            raise ValueError("No parser_kwargs for some of your requisite raw quantities")
        parserkwargsdict = parser_kwargs
    else:
        parserkwargsdict = {q: parser_kwargs for q in all_raws}
    # path
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    ### Let's get started
    missing = []
    datafp = os.path.join(path, "data.json")  # path of complex quantities json file
    if os.path.exists(datafp):
        data = load_js(datafp)
    else:
        data = {}
    for q in qlist:
        if q == str:
            times = 1
        if type(q) in [list,tuple] and len(q)==2:
            times,q = q[0],q[1]
        else:
            raise TypeError("quantity {} not understood".format(q))
        if not q_in_keys(data,q,times):
            miss = raw_quantities(path, qlist=reqs[q], ext=ext, ignore=ignore, parser=parserdict, parser_args=parserargsdict, parser_kwargs=parserkwargsdict)  #Check requisite raw quantities are there (parses if necessary)
            if len(miss) == 0:
                qval = quantfuncs[q](path)  # dictionary of {q: func} that calculate complex quantities (reads specific jsons to get raw, processes, returns)
                data[q] = qval
            else:
                print("cannot calculate {} because of missing raw quantity/ies".format(q))
                missing.append(q)
    dump_js(data,datafp)                        
    return missing

def loopthrough(funcdict,joblist):
    """
    Parameters
    ----------
    funcdict : dict
        {func: {"kw": arg,...}, func2: {"kw": arg,...}, ...}
    joblist : list/str
        list of job tuples. or file to obtain it from

    Returns
    -------
    Performs all the functions in funcdict for all jobs (e.g. complex_quantity, molden file)

    """    
    if type(joblist) == str:
        joblist = get_joblist(joblist)   
    for job in joblist:  # every job is a tuple
        for func,kwargs in funcdict.items():
            ckwargs = kwargs.copy()
            ckwargs.update({"path": os.path.join(*job)})  # we join the tuple elements to obtain the path
            func(**ckwargs)

def collect_data(joblist, levels=["A","B","basis","calc"], qlist="variables.json",reqs=None, ext="*.out", ignore="slurm*", parser=None, parser_args=None, parser_kwargs=None):
    """

    Returns
    -------
    None.

    """
    ### Processing user input
    # qlist
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            jsdict = load_js(qlist)
            qlist = jsdict["complex_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    else:
        raise TypeError("""qlist can be a list, a json filename to extract it from, or a single quantity(str),..],...}""")
    # reqs
    if type(reqs) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            if "jsdict" not in locals():
                jsdict = load_js(qlist)    
            reqs = jsdict["requisites"]
        else:
            reqs = [reqs]
            print("You only gave one raw quantity as requisite. That can happen but is rare. Sure about it?")
    elif type(reqs) != dict:
        raise TypeError("""reqs should be a dictionary {q1: [req1,req2,..], q2: [req1,req2,..],...} \n
                                                        otherwise a json to get it from (as jsdict["reqs"]) or a single quantity(str)""")
    all_raws = set(ittl.chain.from_iterable(reqs.values()))
    # parser
    if type(parser) == dict:
        if np.array([ i not in parser.keys() for i in all_raws]).any():
            raise ValueError("No parser for some of your requisite raw quantities")
        parserdict = parser
    else:
        parserdict = {q: parser for q in all_raws}
    # parser_args
    if type(parser_args) == dict:
        if np.array([ i not in parser_args.keys() for i in all_raws]).any():
            raise ValueError("No parser_args for some of your requisite raw quantities")
        parserargsdict = parser_args
    else:
        parserargsdict = {q: parser_args for q in all_raws}
    # parser_kwargs
    if type(parser_kwargs) == dict:
        if np.array([ i not in parser_kwargs.keys() for i in all_raws]).any():
            raise ValueError("No parser_kwargs for some of your requisite raw quantities")
        parserkwargsdict = parser_kwargs
    else:
        parserkwargsdict = {q: parser_kwargs for q in all_raws}
    # levels
    if len(levels) != len(joblist[0]):
        raise TypeError("Levels and joblist do not match!")
        
    ### Let's get started
    rows = levels + qlist
    columns = []
    for j in joblist:
        path = os.path.join(*j)
        jsfp = os.path.join(path, "data.json")
        cnt = 0
        while cnt <= 1:  # counter to reparse at most once
            data = load_js(jsfp) if os.path.exists(jsfp) else {}
            column = list(j)  
            missing = []
            for n,q in enumerate(qlist):
                to_add = data[q] if q in data.keys() else np.nan
                column.append(to_add)
                if q not in data.keys() and cnt == 0:
                    missing.append(q)
            if len(missing) != 0 and cnt == 0:
                complex_quantities(path=path, qlist=missing, reqs=reqs, ext=ext, 
                                   ignore=ignore, parser=parserdict, parser_args=parserargsdict, parser_kwargs=parserkwargsdict)  # calculated/reparsed
                cnt += 1
            else:
                cnt = 2  # either available or already tried calculating/reparsing
        columns.append(pd.Series(column))
    df = pd.concat(columns, axis=1)  # every column is a j 
    df.index = rows # index
    dt = df.T  # Transpose: now every column is either a level spec(A,B,basis,calc) or a property, every row is a job
    return df

                   
        
