#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some notes:
    1. For the time being, in looping functions (raw_quantities, complex_quantities), 
    ext(i.e. the wildcard for your output files, e.g. "*.out") and ignore (i.e. wildcard for files to ignore, e.g. "slurm*")
    need to be constant. Either launch these functions for different types of outputs,or implement different outputs.

"""

import itertools as ittl
import datetime as dt
import CCParser as ccp
import CCDatabase.QCEasyParser as qcep
import glob as gl
import numpy as np
import pandas as pd
import os
import re
import CCDatabase.utils as ut
import logging
import copy as cp
import shutil as sh
import traceback
from CCDatabase.quantity_functions import ccp_funcs, qcep_ccp_funcs, qcep_funcs

# set up logger. NB avoid homonimity with other module's loggers (e.g. ccp)
ccdlog = logging.getLogger("ccd")
del ccdlog.handlers[:]  # avoid issues with stubborn handlers
ccdlog.setLevel(10)
ccdlog.critical("beware that \"print\" is often delayed, while loggin isn't")

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
    expr = ut.deal_with_type(expr, orig=str, to=lambda x: [x])  # if str turn to list
    if type(expr) is not list:
        raise TypeError("Expr must be either list/tuple or string")
    ### Turning wildcards into regex
    expr = [i.replace(".", "\.").replace("*", ".+").replace("?", ".") for i in expr]    
    ### Handling default startdir if none specified
    startdir = os.getcwd() if not startdir else startdir
    ### Handling default trash if none specified
    trash = os.path.join(os.getcwd(), "Trash") if not trash else trash
    ut.mkdif(trash)
    try:
        histdict = ut.load_js(histfile, cached=False)  # read history file and look for highest ID
        ID = max([int(i) for i in histdict.keys()]) + 1  # get ID and add 1
    except FileNotFoundError:  #history file is not there yet
        histdict = {}
        ID = 1
    moved = []  # will be list of files moved
    for dirname,subdirlist,filelist in os.walk(startdir):   # walking through tree
        for fname in filelist:  # loop over files found
            for e in expr:  # loop over wildcard expressions
                if not(keep^bool(re.match(e, fname))):  # move or not based on "keep"
                    src = os.path.join(dirname, fname)
                    dest = os.path.join(trash, fname)
                    if os.path.exists(dest):
                        splt = dest.split(".")
                        bname = ".".join(splt[:-1])
                        ext = splt[-1]
                        count = 2
                        while os.path.exists(dest):  # avoid overwriting it
                            dest = "{}_{}.{}".format(bname, count,ext)
                            count += 1
                    sh.move(src, dest)
                    print("moved {} to {}".format(os.path.join(dirname,  fname), trash))
                    moved.append((os.path.abspath(src), os.path.abspath(dest)))
    histdict.update({ID:moved})  # {..., ID:[(src1,dest1),(src2,dest2),...]}      
    ut.dump_js(histdict, histfile)          
    with open(logfile, "a") as log:
        neg = "not" if keep else ""
        date = dt.date.today().strftime("%d-%b-%Y")
        exprstr = ",".join(expr)
        log.write("""{}, ID: {}, moved to trash files {} matching with "{}"\n""".format(date, ID, neg, exprstr))

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
    IDs = ut.deal_with_type(IDs, orig=str, to=lambda x: [x])  # if only 1 ID,turn to list
    histdict = ut.load_js(histfile, cached=False)
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
                dest = "{}_{}.{}".format(bname ,count, ext)
                count += 1
            sh.move(src, dest)
            print("restored {}".format(dest))
        del histdict[str(ID)]  # removing entry from histfile
        for n,line in enumerate(lines):
            if "ID: {}".format(ID) in line:
                del lines[n]  # removing entry from logfile
    ut.dump_js(histdict, histfile)  # updating histfile

    with open(logfile, "w") as f:  # updating logfile
        for line in lines:
                f.write(line)   

### Quantities-related
                
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
    data = ut.load_js(fname, cached=False) 
    if "joblist" in data.keys():
        return(data["joblist"])
    else:
        ignore = data["ignore"] if "ignore" in data.keys() else []
        return [i for i in  ittl.product(*[data[i] for i in data["levelnames"]]) if i not in ignore]
              
def find_and_parse(path, ext="*.out", ignore="slurm*", parser=None,
                   parser_args=None, parser_kwargs=None, check_input=True,
                   to_console=True, to_log=False, logname="CCDatabase.log", printlevel=20): 
    """
    Parameters
    ----------
    path: str
        the folder to work in. if None, uses os.getcwd()
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
    check_input: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other functions that check before for loops.
    to_console: bool
        whether to print logging to the console
    to_log: bool
        whether to print logging to file
    logname: str
        the log filename. Default is CCDatabase.log
    printlevel: int
        logging level: 10,20,30,40,50. Accessible also as logging.DEBUG/INFO/WARNING/ERROR/CRITICAL (respectively)
        
    Returns
    -------
    What your parser returns
        the parser container
    """
    ### Logger setup: file or console, level
    ut.setupLogger(to_console=to_console, to_log=to_log, logname=logname, printlevel=printlevel)

    if check_input:
        # path
        path = ut.deal_with_type(path, condition=None, to=os.getcwd)
        path = os.path.abspath(path)
        if not os.path.exists(path):
            ccdlog.critical("{} does not exist!".format(path))
            return None  # nothing to parse
        ccdlog.debug("path is {}".format(path))
        # parser
        if parser is not None:
            try:  # this changes with python versions, hence try/except
                is_func = callable(parser)  # whether it is a function
            except:
                is_func = hasattr(parser, "__call__")  # whether it is a function
            if not is_func:
                raise TypeError(""""Cannot process "parser", as it is none of: dict,list,func,None""")
        ccdlog.debug("parser passed inspection")
        # parser_args
        if parser_args is not None:
            if type(parser_args) is tuple:
                parser_args = list(parser_args)
            elif type(parser_args) is not list:
                parser_args = [parser_args]
                ccdlog.warning("Your \"parser_args\" is neither None nor list/tuple, turning it into a list (parser_args = [parser_args])")
        # parser_kwargs
        if parser_kwargs is not None and type(parser_kwargs) is not dict:
            raise TypeError(""""Cannot process "parser_kwargs", as it is none of: dict,None""") 
    # Finding file
    files = [i for i in gl.glob(os.path.join(path, ext)) if i not in gl.glob(os.path.join(path, ignore))]  # e.g. all *.out which are not slurm*.out
    if len(files) == 0:
        ccdlog.critical("No matching file in {}".format(path))
        return {}
    elif len(files) == 1:
        file = files[0]
        ccdlog.info("File found!")
    else:
        lens = [len(f) for f in files]
        idx = lens.index(min(lens))        
        file = files[idx]
        ccdlog.warning("Several matching files!! Using the shortest file ({})".format(file))
    #Let's get started and parse
    if parser in [None, "ccp", ccp.Parser]:
        parser, parser_args = ccp.Parser, [file]
        if not parser_kwargs: 
            parser_kwargs = dict(software="qchem", overwrite_file=False, 
                                 to_file=True, to_console=False,
                                 to_json=True, json_file="CCParser.json",
                                 large_fn="matrices.npz")
    else:
        parser_args = [file] + parser_args
    ccdlog.debug("parsing!")
    try:
        data = parser(*parser_args,**parser_kwargs)
    except Exception as e:
        ccdlog.critical("Exception while parsing! {}".format(e))
        ccdlog.critical("Traceback: {}".format(traceback.format_exc()))
    return data

def check_qlist(qlist, key, fp, jsdata):
    """
    Note
    ----
    Checks qlist and returns a processed one. Also returns fp and jsdata to avoid double-reading
    
    Parameters
    ----------
    qlist: list
        quantity list(whether raw or complex)
    key: str
        if item in a json file's dictionary (generally raw_quantities/complex_quantities)
    fp: None/str
        the filepath to json
    jsdata: None/obj (dict,list,..)
        what is read in the json file at filepath fp
        
    Returns
    -------
    tuple
        qlist, fp, jsdata
        fp and jsdata only if useful for later (jsonfile contains several arguments)
    """
    if type(qlist) is str:
        if re.match(".+\.json", qlist):  # actually a json file
            fp = cp.copy(qlist)
            jsdata = ut.load_js(fp, cached=False)
            if type(jsdata) is list:
                qlist = cp.copy(jsdata)
                fp, jsdata = None, None  # back to None because it only had qlist
            elif type(jsdata) is dict:
                if key in jsdata.keys():
                    qlist = jsdata[key]
                    ccdlog.info("obtained \"qlist\" from {}".format(fp))
                else:
                    raise KeyError("No {} in the json dictionary".format(key))
            else:  # no tuple option, they become lists in json
                raise TypeError("Your json file is neither dictionary nor list/tuple")
        else:  # actually a single quantity
            qlist = [qlist]
            ccdlog.info("obtained \"qlist\" as [{}]".format(qlist[0]))
    elif type(qlist) is list:
        if False in [type(i) == str or (type(i) == list and len(i) < 4) for i in qlist]:
            if len(qlist) > 4:
                raise TypeError(""""could not process your qlist.
                    Every item should be a string or:
                        - [maxstate, quantity]
                        - [atomicstring, quantity]
                        - [maxstate, atomicstring, quantity]
                        """)
            elif len(qlist) in [2,3]:
                qlist = [qlist]  # 1 quantity with atomiclist and/or stateslist. Will be checked later
            else:
                raise TypeError("""Qlist is a single non-string list. That cannot be processed.""")
    else:
        raise TypeError(""""could not process your qlist.
                        Please provide in one of the following ways:
                            - file.json (will read variable 'raw_quantities'
                            - "exc_energies"
                            - [2,"exc_energies"]
                            - ["exc_energies","osc_strength"]
                            - [[2,"exc_energies"],[2,"osc_strength"]]
                            - ["Na", "EFG_t"]
                            - [["Na", "EFG_t"], ["K", "EFG_t"]]
                            - [2, "Na", "EFG_t"]
                            - [[2, "Na", "EFG_t"], [2, "K", "EFG_t"]]
                            """)
    return qlist,fp,jsdata
 
def check_other_args(fp, jsdata, qlist, parserfuncs, parser, parser_args, parser_kwargs, ex_qs=None, reqs=None, raw=False):
    """
    Note
    ----
    only use: 
        a) within other functions
        b) after running check_qlist()
        c) when fp,jsdata are defined (at least as None)
                                
    Parameters
    ----------
    fp: None/str
        filepath to json file if useful
    jsdata: None/dict
        json dictionary if useful
    qlist: list
        pre-processed qlist
    parserfuncs: {},dict, callable
        {"pname": pfunc/None}
        None => ccp.
        if {pfunc: "pname"} it will be inverted
        callable works only if parser is only one.
    parser: None ,dict, str, func, list(if raw == True)                 
        {"q1": "pname1", "q2":"pname2"} or {"q1": pfunc1, "q2": pfunc2}
        if str/func, uses for all qs
        list allowed only if raw == True
    parser_args: None, dict, list
        {"pname1": args1, pname1, ..}
        list will be used for all parser
    parser_kwargs: None, dict:
        If dictionary of dictionary/None, will be used as is.
        Else, will be used for each parser
    ex_qs: list
        the excited state quantities in qlist
    raw: bool
        whethere it is raw quantities of complex quantities
        
    Returns
    -------
    tuple
        reqsdict, parserfuncs, parserdict, parserargsdict, parserkwargsdict
        reqsdict only for raw == False
        """
    if not raw and None in [ex_qs, reqs]:
        raise ValueError("You must specify ex_qs and reqs")
    if not raw:
        # ex_qs (tuple is fine as well)
        ex_qs = ut.deal_with_type(ex_qs, condition=str, to=lambda x: [x])  #deal with single q
        # reqs
        if reqs is None:
            if fp is not None:  # not declared as nonlocal, i.e. determined json file
                reqs = fp  # we have a json dictionary file
            else:
                raise ValueError("You must specify the required raw quantities for each complex quantity as reqs")
        elif type(reqs) is str:
            if re.match(".+\.json", reqs):  # actually a json file (qlist or its own)
                if jsdata is not None and reqs == fp:  # not already read
                    jsdata = ut.load_js(reqs, cached=False)
                if type(jsdata) is not dict:
                    raise ValueError("Your json file is not a dictionary!!")
                reqsdict = jsdata["requisites"] if "requisites" in jsdata.keys() else jsdata
                ccdlog.info("obtained \"reqsdict\" from {}".format(reqs))
            else:
                reqsdict = {q: [reqs] for q in qlist}
                ccdlog.warning("""You only gave one raw quantity as requisite. That can happen but is rare. Sure about it?
                            It is being converted as {q: [reqs] for q in qlist}""")
        elif type(reqs) in [list,tuple]: 
            reqs = ut.deal_with_type(reqs, condition=tuple, to=list)
            reqsdict = {q: reqs for q in qlist}
            ccdlog.info("You gave a list/tuple as \"reqs\". Obtained reqsdict as reqs = {q: reqs for q in qlist}")
        elif type(reqs) is dict:
            reqsdict = reqs
            for k,v in reqsdict.items():
                reqsdict[k] = ut.deal_with_type(v, condition=str, to=lambda x: [x])
            ccdlog.info("Turned any string value in your dictionary into a list")
        else:
            raise TypeError("""reqs should be a dictionary {q1: [req1,req2,..], q2: [req1,req2,..],...} \n
                                                            otherwise a json to get it from (as jsdict["reqs"]) or a single quantity(str)""")
        tmp = np.array([q in reqsdict.keys() for q in qlist])
        if not tmp.all():
            issues = np.where(tmp == False)[0]
            for i in issues:
                par_q = "_".join(qlist[i].split("_")[:-1])  # parent q, e.g. "ex_en_11" => "ex_en"
                if par_q in reqsdict.keys():
                    reqsdict[qlist[i]] = reqsdict[par_q]
                    ccdlog.info("You gave {} in qlist and {} in reqs. I suppose they are the same quantity, and will use them as such".format(qlist[i],par_q))
                else:
                    raise ValueError("Some quantity does not have its requisites specified")
        # any raw which is a req of a complex
        all_raws = [r[1] if type(r) in [tuple, list] else r for r in set(ittl.chain.from_iterable(reqsdict.values()))]  
    rawlist = qlist if raw else all_raws
    # parserfuncs
    """
    In the end it must be:
    parserfuncs = {"parsername1": parserfunc1,
                   "parsername2": parserfunc2, .. }
    None is interpreted as ccp.Parser
    """
    if not parserfuncs:  # default
        parserfuncs = {"ccp": None,
                       "qcep": qcep.parse_qchem}
    else:
        if type(parserfuncs) != dict:
            check_if_func = True
            if type(parserfuncs) in [tuple, list]:
                if len(parserfuncs) == 1:  # actually only one parser                                                      
                    parserfuncs = parserfuncs[0]
                    ccdlog.info("\"parserfuncs\" was a single element list. Will be used for all quantities")   
                if raw:
                    if type(parser) not in [tuple, list]:
                        raise ValueError("\"parserfuncs\" is a list but \"parser\" is not. Could not establish pname-pfunc correspondence")
                    if len(parserfuncs) == len(parser):                                                                       
                        parserfuncs = {parser[n]: parserfuncs[n] for n in range(len(qlist))}  
                        check_if_func = False                                   
                        ccdlog.info("obtained \"parserdict\"")                                                            
                    else:                                                                                                 
                        raise ValueError("""You gave a list as "parser" but the length does not match that of "qlist" """)
                else:
                    raise TypeError("\"parser\" cannot be a list. it can be dict/str/None")
            if check_if_func:
                try:  # this changes with python versions, hence try/except
                    is_func = callable(parserfuncs)  # whether it is a function
                except:
                    is_func = hasattr(parserfuncs,"__call__")  # whether it is a function
                if is_func: 
                    if type(parser) == str:
                        if parser == "ccp" and parserfuncs.__name__ != "Parser":
                            ccdlog.warning("parser=\"ccp\" would call CCParser.Parser, \
                                           but you provided another function. Your function will be used!")
                        parsername = parser  
                    elif parser is not None:
                        if type(parser) in [list, tuple]:
                            all_same = np.array([i == parser[0] for i in parser]).all() 
                        if type(parser) == dict:
                            all_same = np.array([i == list(parser.keys())[0] for i in parser.vals()]).all() 
                        if not all_same:
                            raise ValueError("You gave only one parserfunc but called for several parsernames!!")
                    else:
                        ccdlog.warning("parser=None would call CCParser.Parser, \
                                           but you provided another function. Your function will be used!")
                        parsername = parserfuncs.__name__
                    parserfuncs = {parsername: parserfuncs}
        try:  # by now it has to be a dictionary
            all_vals_funcs = np.array([callable(p) or (p==None) for p in parserfuncs.values()]).all()
            all_keys_funcs = np.array([callable(p) or (p==None) for p in parserfuncs.keys()]).all()
        except:
            all_vals_funcs = np.array([hasattr(p,"__call__") or (p==None) for p in parserfuncs.values()]).all()
            all_keys_funcs = np.array([hasattr(p,"__call__") or (p==None) for p in parserfuncs.keys()]).all()
        assert all_vals_funcs^all_keys_funcs, "If you wish to use other parsers than ccp&qcep, please provide parserfuncs={pname1: pfunc1, pname2: pfunc2,..}"
        if all_keys_funcs:
            parserfuncs = {parserfuncs[k]: k for k in parserfuncs.keys()}  # inverting if necessary
    # parser
    """
    In the end it must be:
        parserdict = {"q1": "parsername3",
                      "q2": "parsername4", ..}
    """
    if parser is None:  
        parserdict = {q: "ccp" for q in rawlist}
        ccdlog.info("obtained \"parserdict\" as \"ccp\" for all raw quantities")
    else:  #parser != None
        if type(parser) is str:
            toprint = parser
            parser = {q: parser for q in rawlist}  # not parserdict!!
            ccdlog.info("obtained \"parserdict\" as \"{}\" for all raw quantities".format(toprint))
        if type(parser) is dict:
            tmp = np.array([ i in parser.keys() for i in rawlist])
            if not tmp.all():
                issues = np.where(tmp == False)[0]
                for i in issues:
                    par_q = "_".join(rawlist[i].split("_")[:-1])  # parent q, e.g. "ex_en_11" => "ex_en"
                    if par_q in parser.keys():
                        parser[rawlist[i]] = parser[par_q]
                        ccdlog.info("You gave {} in {} and {} in parser. I suppose they are the same quantity, and will use them as such".format(rawlist[i],
                                    "qlist" if raw else "requisites",par_q))
                    else:
                        raise ValueError("No parser for some of your{} raw quantities".format("" if raw else " requisite"))
            parserdict = parser
            ccdlog.info("obtained \"parserdict\" with all necessary keys")
        elif type(parser) in [list,tuple]:
            if len(parser) == 1:  # actually only one parser                                                      
                    parserdict = {q: parser[0] for q in qlist} 
                    ccdlog.info("\"parserdict\" was a single element list. Will be used for all quantities")   
            if raw:
                if len(parser) == len(qlist):                                                                       
                    parserdict = {qlist[n]: parser[n] for n in range(len(qlist))}                                     
                    ccdlog.info("obtained \"parserdict\"")                                                            
                else:                                                                                                 
                    raise ValueError("""You gave a list as "parser" but the length does not match that of "qlist" """)
            else:
                raise TypeError("\"parser\" cannot be a list. it can be dict/str/None")
        else:
            try:  # this changes with python versions, hence try/except
                is_func = callable(parser)  # whether it is a function
            except:
                is_func = hasattr(parser,"__call__")  # whether it is a function
            if is_func: 
                name = parser.__name__
                parserfuncs[name] = parser
                parserdict = {q: name for q in rawlist}
                ccdlog.warning("You gave a function as \"parser\". Will be named {} and used for all quantities.".format(name))
                ccdlog.warning("Please consider giving parserfuncs={pname: pfunc} and parser=pname")
            else:
                raise TypeError("Cannot process \"parser\", as it is none of: {}dict,func,None".format("list/tuple," if raw else ""))
    for k,v in parserdict.items():
        if type(v) is str:
            if v not in parserfuncs.keys():
                raise ValueError("You asked at least once for an unrecognised parser.")
            continue
        elif v is None:
                parserdict[k] = "ccp"  # changing None to ccp
                continue
        elif v is False:
            continue
        else: 
            try:
                is_func = callable(v)
            except:
                is_func = hasattr(v, "__call__")
            if is_func:
                name = v.__name__
                parserfuncs[name] = v
                parserdict[k] = name
                ccdlog.warning("You gave q: func instead of q: parsername. We changed it to {q}:{name} and added {name} to parserfuncs".format(q=k,name=name))
            else:
                raise ValueError("some value in your parserdict is not valid. They should all be parsernames or False.")
    parsers_used = list(set(parserdict.values()))
    # parser_args
    """
    In the end it must be:
        parserargsdict = {"pname1": [arg1, arg2, ..],
                          "pname2": [arg10, arg11, ..], ..}
    None is accepted as *args for ccp.Parser
    """
    if parser_args is None:
        parserargsdict = {p: [] for p in parsers_used}
        ccdlog.info("obtained \"parserargsdict\" as empty list for all parsers")
    elif type(parser_args) is dict:
        if np.array([type(k) != str for k in parser_args.keys()]):
            for k,v in parser_args.items():
                if type(k) != str:
                    try:
                        is_func = callable(k)
                    except:
                        is_func = hasattr(k, "__call__")
                    if is_func:
                        name = k.__name__
                        parser_args[name] = v
                        parserfuncs[name] = k
                        ccdlog.warning("You gave pfunc: args instead of \"pname\": args. We changed it to \"pname\": args and added \"pname\": pfunc to parserfuncs")
                    else:
                        raise ValueError("At least one key in parser_args is neither a string nor a function")
        if np.array([p not in parser_args.keys() for p in parsers_used]).any():  # not elif!
            raise ValueError("No parser_args for at least one of your parsers")
        typelist = [type(parser_args[p]) for p in parsers_used]
        if np.array([t not in [list, tuple, type(None)] for t in typelist]).any():
            raise TypeError("At least one of the elements in parser_args is none of list,tuple,None")
        elif tuple in typelist:
            for k,v in parser_args.items():
                if type(v) is tuple:  # turn to list
                    parser_args[k] = list(v)
                    ccdlog.debug("turned parser_args element from tuple to list")
        parserargsdict = parser_args
        ccdlog.info("obtained \"parserargsdict\" from the provided dictionary")
    elif type(parser_args) in [list,tuple]:
        parser = ut.deal_with_type(parser, condition=tuple, to=list)  
        parserargsdict = {p: parser_args for p in parsers_used}
        ccdlog.info("obtained \"parserargsdict\" as your \"parser_args\" for all parsers")
    else:
        raise ValueError("""Cannot process your "parser_args" """)
    # parser_kwargs    
    """
    In the end it must be:
        parserkwargsdict = {"pname1": {kw1: arg1, kw2: arg2, ..},
                          "pname2": {kw10: arg10, kw11: arg11, ..}, ..}
    None is accepted as **kwargs for ccp.Parser
    """
    if parser_kwargs is None:
        parserkwargsdict = {p: {} for p in parsers_used}
        ccdlog.info("obtained \"parserkwargsdict\" as empty dictionary for all parsers")
    else:  # parser_kwargs != None
        if type(parser_kwargs) in [tuple, list]:
            if len(parser_kwargs) == 1 and type(parser_kwargs[0]) is dict:
                parser_kwargs = parser_kwargs[0]
                ccdlog.debug("You gave \"parser_kwargs\" as [dict]. The dictionary will be processed.")
            else:
                raise TypeError("Please provide \"parser_kwargs\" either as a dictionary (to use for all parsers) or as a dictionary of dictionaries")
        if type(parser_kwargs) is dict:  # not elif!!
            if np.array([type(i) not in [dict,None] for i in parser_kwargs.values()]).any():  # It's not a dict of dicts/None
                parserkwargsdict = {p: parser_kwargs for p in parsers_used}
                ccdlog.warning("You gave a dictionary, instead of a dictionary of dictionaries, as  \"parser_kwargs\". It will be used for all parsers.")
            else:
                if np.array([type(k) != str for k in parser_kwargs.keys()]).any():
                    for k,v in parser_args.items():
                        if type(k) != str:
                            try:
                                is_func = callable(k)
                            except:
                                is_func = hasattr(k, "__call__")
                            if is_func:
                                name = k.__name__
                                parser_kwargs[name] = v
                                parserfuncs[name] = k
                                ccdlog.warning("You gave pfunc: args instead of \"pname\": kwargs. We changed it to \"pname\": kwargs and added \"pname\": pfunc to parserfuncs")
                            else:
                                raise ValueError("At least one key in parser_args is neither a string nor a function")
                if np.array([p not in parser_kwargs.keys() for p in parsers_used]).any():  # allows unnecessary keys. Not elif!
                    raise ValueError("""No parser_kwargs for at least one of your parsers""")
                else:
                    parserkwargsdict = parser_kwargs
                    ccdlog.info("obtained \"parserkwargsdict\" as your \"parser_kwargs\"")
        else:
            raise ValueError(""""parser_kwargs" can be None, dict (to use for all quantities), or dict of dict. Your type is not recognised.""")
    to_return = (parserfuncs, parserdict, parserargsdict, parserkwargsdict) if raw else (reqsdict, parserfuncs, parserdict, parserargsdict, parserkwargsdict)
    return to_return

def get_stateslist_atomiclist(qlist):
    """
    """
    stateslist, atomiclist, nqlist = [], [], []
    for q in qlist:
        if type(q) == str:
            stateslist.append(False)
            atomiclist.append(False)
            nqlist.append(q)
        elif type(q) in [list, tuple]:
            if type(q[-1]) != str:
                raise ValueError("Provide stateslist/atomiclist first, quantity last!")
            if len(q) == 2:
                if type(q[0]) == str:
                    stateslist.append(False)
                    atomiclist.append(q[0])
                    nqlist.append(q[1])
                elif type(q[0]) == int:
                    stateslist.append(q[0])
                    atomiclist.append(False)
                    nqlist.append(q[1])
            elif len(q) == 3:
                if type(q[0]) == str and type(q[1]) == int:
                    stateslist.append(q[1])
                    atomiclist.append(q[0])
                    nqlist.append(q[2])
                elif type(q[1]) == str and type(q[0]) == int:
                    stateslist.append(q[0])
                    atomiclist.append(q[1])
                    nqlist.append(q[2])
                else:
                    raise ValueError("Something wrong with one element of qlist")
        else:
            raise ValueError("Something wrong with one element of qlist")
    return stateslist, atomiclist, nqlist
                
def raw_quantities(path=None, qlist="variables.json", ext="*.out", ignore="slurm*",
                   parser_file="CCParser.json", parserfuncs={},parser=None, 
                   parser_args=None, parser_kwargs=None, check_input=True,
                   to_console=True, to_log=False, logname="CCDatabase.log", 
                   printlevel=20):
    """
    Parameters
    ----------
    path: str
        the path to check 
    qlist : list/str
        raw quantities to check. Default checks from variables.json
        If there is no ",", searches for the quantity directly in the parser's json file in path.
        Otherwise, the string before the comma is considered as the folder name to look into.
        This folder can be a subfolder (e.g. F&T, Macrocycles, Prepol) or a parallel folder (e.g. "../sup/").
        One can also give q as "folder, file, quantity"
    ext: str
        wildcard-like expression to match. Default is "*.out"
    ignore: str
        wildcard-like expression of files to ignore. Default "slurm*"
    parser_file: str
        name of the parser's json file. Default is "CCParser.json"
    parserfuncs: {},dict
        {"pname": pfunc/None}
        None => ccp.
        if {pfunc: "pname"} it will be inverted
    parser: None ,dict, str, func, list(if raw == True)                 
        {"q1": "pname1", "q2":"pname2"} or {"q1": pfunc1, "q2": pfunc2}
        if str/func, uses for all qs
        list allowed
    parser_args: None, dict, list
        {"pname1": args1, "pname2": args2, ..}
        list will be used for all parser
    parser_kwargs: None, dict:
        If dictionary of dictionary/None, will be used as is.
        Else, will be used for each parser
    check_input: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other functions that check before for loops.
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
    checks if the quantities are available in all "jobs". If necessary, it parses.
    Then prints what is missing
    
    Returns
    -------
    list
        the quantities missing in path
    """        
    ### Logger setup: file or console, level
    ut.setupLogger(to_console=to_console, to_log=to_log, logname=logname, printlevel=printlevel)
    fp, jsdata = None, None  # to allow nonlocal declaration
    ### Processing user input
    if check_input:
        qlist, fp, jsdata = check_qlist(qlist, "raw_quantities", fp, jsdata)
    nvalslist = [q[0] if type(q) in [tuple,list] else 1 for q in qlist]
    qlist = [q[1] if type(q) in [tuple,list] else q for q in qlist]
    ccdlog.debug("divided qlist into qlist and nvalslist")
 
    if check_input:
        # path
        if path is None:
            path = os.getcwd()
        else:
            path = os.path.normpath(path)
        if not os.path.exists(path):
            ccdlog.critical("{} does not exist!".format(path))
            return qlist  # all quantities are missing!
        ccdlog.debug("path is {}".format(path))
        parserfuncs, parserdict, parserargsdict, parserkwargsdict = check_other_args(fp,
                                                                                     jsdata, 
                                                                                     qlist, 
                                                                                     parserfuncs,
                                                                                     parser, 
                                                                                     parser_args, 
                                                                                     parser_kwargs,
                                                                                     raw=True)
    else:
        parserdict = parser
        parserargsdict = parser_args
        parserkwargsdict = parser_kwargs
        ccdlog.debug("Assigned parser-related arguments without checking")
    ### Let's get started      
    missing = []
    reparsed = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    data = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    for n,q in enumerate(qlist):
        parsername = parserdict[qlist[n]]
        ccdlog.debug("quantity: {}".format(q))
        nvals = nvalslist[n]
        """
        q can be as follows:
            1 "quant": checks for "quant" in the current jobfolder
            2 "subfol,quant": checks for "quant" in subfol
            3 "paralfol,quant": checks for "quant" in paralfol
            4 "non_stdfile,quant": checks for "quant" in "non_stdfile" in current jobfolder
            5 "fol,non_stdfile,quant" or "paralfol,file,quant": loks for "quant" in "non-stdfile" in sub/paral-fol
        """
        if "," in q:  # in another calc (not case 1)
            splt = q.split(",")
            if len(splt) == 2 and "." in splt[0]:  # case 4
                ccdlog.debug("case 4")
                fname = splt[0]
                path_tmp = path
            else:  # cases 2,3,5
                fol = splt[0]
                if len(splt) == 3:  # case 5
                    ccdlog.debug("case 5")
                    fname = splt[1]
                else:  # cases 2,3
                    fname = parser_file
                as_subfol = os.path.join(path, fol)
                is_subfol = os.path.isdir(as_subfol)
                as_paralfol = os.path.join(ut.split_path(path)[0], fol)
                is_paralfol = os.path.isdir(as_paralfol)
                if is_subfol and not is_paralfol:  # it is a subdir (case 2)
                    ccdlog.debug("case 2")
                    path_tmp = as_subfol
                elif not is_subfol and is_paralfol:  # it is a parallel dir (case 3)
                    ccdlog.debug("case 3")
                    path_tmp = as_paralfol
                else:
                    ccdlog.critical("the location of your quantity {} cannot be understood.\
                             Most likely it is either not present or double. Continuing with other quantities".format(qlist[n]))
                    continue
            q = splt[-1]
            filepath = os.path.join(path_tmp, fname)
            ccdlog.debug("filepath is {}".format(filepath))
        else:  # case 1
            ccdlog.debug("case 1")
            path_tmp = path
            filepath = os.path.join(path_tmp, parser_file)
        ### case is determined. path_tmp and filepath are set
        if path_tmp not in data.keys():
            ccdlog.debug("No data from {} yet".format(path_tmp))
            data[path_tmp] = ut.dict_from_file(filepath, cached=True)
            if path_tmp not in reparsed.keys():
                reparsed[path_tmp] = []
        if not ut.rq_in_keys(data[path_tmp], q, nvals=nvals) and parsername not in reparsed[path_tmp]:  # quantity not in json, not reparsed yet
            if parsername != False:  # no parsing for human-generated quantities (e.g. correspondance)
                # Here parserdict should be {q1: pname1,..}, parserargsdict {q1:[arg1,..],..}, parserkwargsdict {q:{kw1:arg1,..},..}
                ccdlog.info("reparsing in folder {}".format(path_tmp))
                parser = parserfuncs[parsername]
                try:
                    find_and_parse(path_tmp, ext=ext, ignore=ignore,
                                   parser=parser, 
                                   parser_args=parserargsdict[parsername], 
                                   parser_kwargs=parserkwargsdict[parsername],
                                   check_input=False, to_console=to_console,
                                   to_log=to_log, logname=logname, printlevel=printlevel)  # parser takes care of dumping json
                except:
                    ccdlog.critical("could not parse in {}".format(path_tmp))
                if os.path.exists(filepath):
                    data[path_tmp].update(ut.load_js(filepath, cached=False))  # let's read reparsed json/xlsx/df
                else:
                    ccdlog.critical("{} does not exist. Either did not parse or parsed and saved elsewhere".format(filepath))
            reparsed[path_tmp].append(parsername) 
        if not ut.rq_in_keys(data[path_tmp], q, nvals=nvals) and parsername in reparsed[path_tmp]:
            missing.append(qlist[n])  # original q, not split
            ccdlog.error("{} missing {}".format(path_tmp,q))
    return missing   

def complex_quantities(path=None, qlist="variables.json", ex_qs=[], reqs=None, ext="*.out",
                       ignore="slurm*", parser_file="CCParser.json",parserfuncs={},
                       parser=None, parser_args=None, parser_kwargs=None, check_input=True,
                       funcdict="ccp",to_console=True, to_log=False,
                       logname="CCDatabase.log", printlevel=20):
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
    parserfuncs: {},dict
        {"pname": pfunc/None}
        None => ccp.
        if {pfunc: "pname"} it will be inverted
    parser: None ,dict, str, func, list(if raw == True)                 
        {"q1": "pname1", "q2":"pname2"} or {"q1": pfunc1, "q2": pfunc2}
        if str/func, uses for all qs
        list allowed only if raw == True
    parser_args: None, dict, list
        {"pname1": args1, pname1, ..}
        list will be used for all parser
    parser_kwargs: None, dict:
        If dictionary of dictionary/None, will be used as is.
        Else, will be used for each parser
    check_input: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other functions that check before for loops.
    funcdict: dict/str
        dictionary {"complex_quantity": function, ...}
        or one of these strings: "ccp","qcep-ccp","qcep"
        strings recognition is case insensitive and removes non-alphabetic characters
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
    checks if the quantities are available in all "jobs". If necessary, it parses.
    Then prints which folders do not have.
    
    Returns
    -------
    list
        [quantity1,quantity2, ...] of what is missing
    """
    ### Logger setup: file or console, level
    ut.setupLogger(to_console=to_console, to_log=to_log, logname=logname, printlevel=printlevel)
    fp, jsdata = None, None  # to allow nonlocal declaration
    ### Processing user input
    #funcdict
    if type(funcdict) is str:
        funcdict = (''.join([i for i in funcdict if i.isalpha()])).lower()  # removes non-aplhabetic
        if funcdict == "ccp":
            funcdict = ccp_funcs  
            ccdlog.info("imported ccp_funcs as funcdict")
        elif funcdict == "qcepccp":
            funcdict = qcep_ccp_funcs
            ccdlog.info("imported qcep_ccp_funcs as funcdict")
        elif funcdict == "qcep":
            funcdict = qcep_funcs
            ccdlog.info("imported qcep_funcs as funcdict")
    elif funcdict != dict:
        raise TypeError(""""funcdict" should be a dictionary or a string""")
        
    if check_input:
        qlist, fp, jsdata = check_qlist(qlist, "complex_quantities", fp, jsdata)
    stateslist, atomiclist, qlist = get_stateslist_atomiclist(qlist)
    ccdlog.debug("divided qlist into qlist, stateslist, atomiclist")

    if check_input:
        # path
        if path is None:
            path = os.getcwd()
        else:
            path = os.path.normpath(path)
        if not os.path.exists(path):
            ccdlog.critical("{} does not exist!".format(path))
            return qlist  # all quantities are missing!
        ccdlog.debug("path is {}".format(path))
        reqsdict, parserfuncs, parserdict, parserargsdict, parserkwargsdict = check_other_args(fp,
                                                                                               jsdata, 
                                                                                               qlist, 
                                                                                               parserfuncs,
                                                                                               parser, 
                                                                                               parser_args, 
                                                                                               parser_kwargs,
                                                                                               ex_qs=ex_qs, 
                                                                                               reqs=reqs,
                                                                                               raw=False)
    else:  # check_input==False
        reqsdict = reqs
        parserdict = parser
        parserargsdict = parser_args
        parserkwargsdict = parser_kwargs
        ccdlog.debug("Assigned parser-related arguments without checking")
    ### Let's get started
    missing = []
    datafp = os.path.join(path, "data.json")  # path of complex quantities json file
    if os.path.exists(datafp):
        try:
            data = ut.load_js(datafp, cached=True)
            ccdlog.info("Loaded {}".format(datafp))
        except:
            data = {}
            ccdlog.warning("issues in reading  {}!".format(datafp))
    else:
        data = {}
        ccdlog.debug("{} doesn't exist yet".format(datafp))
    for n,q in enumerate(qlist):
        ccdlog.debug("quantity {}".format(q))
        if q not in funcdict.keys(): # dictionary of {q: func} that calculate complex quantities (reads specific jsons to get raw, processes, returns)
            par_q = "_".join(q.split("_")[:-1])  # parent q, e.g. "ex_en_11" => "ex_en"
            if par_q in funcdict.keys():
                func = funcdict[par_q]
                state_num = int(q.split("_")[-1])  # the state(only one) we want to look at
                ccdlog.info("{}: will use function for {} on state {}".format(q, par_q, state_num))
            else:
                raise ValueError("There is no function for complex quantity {}".format(q))
        else:
            func = funcdict[q]
            state_num = False
        states = stateslist[n]
        if not ut.cq_in_keys(data, q,  states=states, exc=(q in ex_qs)):
            ccdlog.debug("{} not in data{}. Trying to obtain it".format(q, " {} times".format(states) if states else ""))
            miss = raw_quantities(path, qlist=reqsdict[q], ext=ext, ignore=ignore,
                                  parser_file=parser_file,
                                  parserfuncs=parserfuncs,
                                  parser=parserdict,
                                  parser_args=parserargsdict,
                                  parser_kwargs=parserkwargsdict,
                                  check_input=False, to_console=to_console,
                                  to_log=to_log, logname=logname, printlevel=printlevel)  #Check requisite raw quantities are there (parses if necessary)
            if len(miss) == 0:
                ccdlog.debug("no requisite raw quantities missing")
                shift = 1 if q in ex_qs else 0  # 1 if only excited state
                if states:
                    vals={}
                    for s in range(shift, states+1):  # if ES, python counting => fortran/human counting
                        try:
                            if atomiclist[n]:
                                qval = func(path=path, atomstring=atomiclist[n], n=s)
                            else:
                                qval = func(path=path, n=s)
                            vals[s] = qval  # nb will become a str when dumped
                            data["{}_{}".format(q, s)] = qval  # e.g. ex_en_1 : val1
                        except Exception as e:
                            ccdlog.error("Errors while calculating {}, state {}, in {}".format(q, s, path))
                            ccdlog.error("Error message below \n {}".format(e))
                            missing.append(q)
                    data[q] = vals  # e.g ex_en:{1:val1,2:val2}
                    ccdlog.debug("saved {} as both dictionary and individual values".format(q))
                elif state_num:
                    try:
                        if atomiclist[n]:
                            qval = func(path=path, atomstring=atomiclist[n], n=s)
                        else:
                            qval = func(path=path, n=state_num)
                        data[q] = qval
                    except Exception as e:
                        ccdlog.error("Errors while calculating {} in {}".format(q, path))
                        ccdlog.error("Error message below \n {}".format(e))
                        missing.append(q)
                else:  # no state specified, could be 1 or "as many as possible"
                    failed = False
                    s = shift
                    vals = {}
                    while not failed:  # trying to obtain as many vals as possible
                        try:
                            if atomiclist[n]:
                                qval = func(path=path, atomstring=atomiclist[n], n=s)
                            else:
                                qval = func(path=path, n=s)
                            vals[s] = qval  # nb will become a str when dumped
                            ccdlog.debug("q:{}, s:{} = {}".format(q, s, qval))
                            s += 1
                        except Exception as e:  # max num of vals
                            if s == 0:  # ex_q not declared, e.g. ex_en_0
                                failed = None
                                s += 1
                                ccdlog.info("""Probably q had to be added to ex_qs but was not.
                                             Error while calculating {}_{} in {}. Error message below.\
                                             Traceback if printlevel=10 \n {}""".format(q, s-1, path, e))
                                ccdlog.debug("Traceback:\n {}".format(traceback.format_exc()))
                                
                            else:
                                failed = True 
                                ccdlog.info("Error while calculating {}_{} in {}. Error message below.\
                                            Traceback if printlevel=10 \n {}""".format(q, s, path, e))
                                ccdlog.debug("Traceback:\n {}".format(traceback.format_exc()))
                    if len(vals) == 0:
                        missing.append(q)
                    elif len(vals) == 1:
                        data[q] = vals[0]
                        ccdlog.debug("{} is a single value".format(q))
                    else:
                        data[q] = vals
                        for s,v in vals.items():
                            data["{}_{}".format(q, s)] = v  # e.g. ex_en_1 : val1
                        ccdlog.debug("{} has many values. Added both dict and individual values".format(q))
            else:
                ccdlog.error("cannot calculate {} because of missing raw quantity/ies".format(q))
                ccdlog.debug("missing: {}".format(miss))
                missing.append(q)
    ut.dump_js(data, datafp)
    ccdlog.info("dumped into {}".format(datafp))                        
    return missing

def collect_data(joblist, levelnames=["A","B","basis","calc"], qlist="variables.json",
                 ex_qs=[], reqs=None, ext="*.out", ignore="slurm*", parser=None,
                 parserfuncs={},parser_file="CCParser.json",parser_args=None, 
                 parser_kwargs=None, check_input=True, funcdict="ccp",
                 look_for_more_states=True, to_console=True, to_log=False, 
                 logname="CCDatabase.log", printlevel=20):
    """
    Parameters
    ----------
    joblist: list[tuples/list]/str
        Either list of jobs, or json file to retrieve
    levelnames: list[str]
        list of level names.Must match len of joblist items. Default is ["A","B","basis","calc"]
    qlist: list[str]/str
        list of complex quantities to get
    reqs: dict
        dictionary of requisite raw quantities for each complex quantity
    ext: str
        wildcard expression to match for find_and_parse 
    ignore: str
        wildcard expression to ignore for find_and_parse 
    parserfuncs: {},dict
        {"pname": pfunc/None}
        None => ccp.
        if {pfunc: "pname"} it will be inverted
    parser: None ,dict, str, func,# list(if raw == True)                 
        {"q1": "pname1", "q2":"pname2"} or {"q1": pfunc1, "q2": pfunc2}
        if str/func, uses for all qs
        # list allowed only if raw == True
        # TODO: check if list works (possible only if reqs is also a list)
    parser_args: None, dict, list
        {"pname1": args1, pname1, ..}
        list will be used for all parser
    parser_kwargs: None, dict:
        If dictionary of dictionary/None, will be used as is.
        Else, will be used for each parser
    funcdict: dict/str
        dictionary {"complex_quantity": function, ...}
        or one of these strings: "ccp","qcep-ccp","qcep"
        strings recognition is case insensitive and removes non-alphabetic characters
    to_console: bool
        whether to print logging to the console
    to_log: bool
        whether to print logging to file
    logname: str
        the log filename. Default is CCDatabase.log
    printlevel: int
        logging level: 10,20,30,40,50. Accessible also as logging.DEBUG/INFO/WARNING/ERROR/CRITICAL (respectively)
        
   Returns
    -------
    pd.DataFrame
        a DataFrame with the desired quantities for each calculation
    """
    ### Logger setup: file or console, printlevel
    ut.setupLogger(to_console=to_console, to_log=to_log, logname=logname, printlevel=printlevel)
    fp, jsdata = None, None
    ### Processing user input
    # joblist
    joblist = ut.deal_with_type(joblist,condition=str, to=get_joblist)
    ccdlog.info("obtained \"joblist\"")
    if check_input:
        qlist, fp, jsdata = check_qlist(qlist, "complex_quantities", fp, jsdata)
    oldqlist = qlist.copy()
    stateslist, atomiclist, qlist = get_stateslist_atomiclist(qlist)
    ccdlog.debug("divided qlist into qlist, stateslist, atomiclist")

    if check_input:
        # funcdict
        if type(funcdict) is str:
            funcdict = (''.join([i for i in funcdict if i.isalpha()])).lower()  # removes non-aplhabetic
            if funcdict not in ["ccp","qcepccp","qcep"]:
                raise ValueError("Unrecognised string value for \"funcdict\"")
        elif type(funcdict) is not dict:
            raise TypeError("\"funcdict\" can be a dictionary or a recognisable string (\"ccp\",\"qcepccp\",\"qcep\")")            
        reqsdict, parserfuncs, parserdict, parserargsdict, parserkwargsdict = check_other_args(fp,
                                                                                               jsdata, 
                                                                                               qlist, 
                                                                                               parserfuncs,
                                                                                               parser, 
                                                                                               parser_args, 
                                                                                               parser_kwargs,
                                                                                               ex_qs=ex_qs, 
                                                                                               reqs=reqs,
                                                                                               raw=False)
    else:  # check_input == False
        reqsdict = reqs
        parserdict = parser
        parserargsdict = parser_args
        parserkwargsdict = parser_kwargs
        ccdlog.debug("Assigned parser-related arguments without checking")
    # levelnames
    if len(levelnames) != len(joblist[0]):
        raise TypeError("levelnames and joblist do not match!")
    
    ### Let's get started
    ccdlog.debug("done checking input. Starting collection")
    nsl = stateslist.copy()  # new stateslist
    natmklst = [[1 for i in s.split(",")] if s else [1] for s in atomiclist]  # max natoms per atomkey
    columns = []
    for j in joblist:
#        ut.clear_caches(ut.caches)
        path = os.path.join(*j)
        if not os.path.exists(path):
            ccdlog.critical("{} does not exist! skipping it".format(path))
            continue
        ccdlog.info("working on {}".format(path))
        jsfp = os.path.join(path, "data.json")
        cnt = 0
        while cnt <= 1:  # counter to reparse at most once
            data = ut.load_js(jsfp, cached=True) if os.path.exists(jsfp) else {}
            column = list(j)  
            missing = []
            if data == {}:  # no need to loop over qlist
                if cnt == 0:  # let's obtain our complex quantities
                    ccdlog.debug("data.json was not there or empty")
                    missing = oldqlist
                else:
                    ccdlog.critical("Your datafile is empty after parsing. Probably filename mismatch or nothing written to file")
            else:
                for n,q in enumerate(qlist): 
                    tmp = natmklst[n].copy()  # natoms for each atomkey before this iteration
                    ccdlog.debug("quantity: {}".format(q))
                    shift = 1 if q in ex_qs else 0  # 1 if only excited state
                    item = data[q] if q in data.keys() else False  # can be value, or many values!
                    if stateslist[n]:  # NB not "is not False"
                        try:
                            if atomiclist[n]:
                                values = []
                                for s in range(shift, stateslist[n]+1):
                                    for natmk,atomkey in enumerate(atomiclist[n].split(",")):
                                        to_ext = item[str(s)][atomkey]
                                        to_ext.extend((natmklst[n][natmk] - len(to_ext))*[np.nan])  # if fewer atoms than other job
                                        values.extend(to_ext)
                                        if len(to_ext) > natmklst[n][natmk]:
                                                natmklst[n][natmk] = len(to_ext)
                            else:
                                values = [item[str(v)] for v in range(shift, stateslist[n]+1)]
                            # if q or str[v] gives key error goes to except
                            column.extend(values)
                            ccdlog.debug("gotten values for {} from its dictionary".format(q))
                        except:
                            missing_states, partial_states, missing_atmks = [], [], []
                            for s in range(shift,stateslist[n]+1):
                                qn = "{}_{}".format(q,s)
                                if qn in data.keys():
                                    if atomiclist[n]:
                                        values = []
                                        for natmk,atom_key in enumerate(atomiclist[n].split(",")):
                                            if atom_key in data[qn].keys():
                                                to_ext = data[qn][atomkey]
                                                to_ext.extend((natmklst[n][natmk] - len(to_ext))*[np.nan])  # if fewer atoms than other job
                                                values.extend(to_ext)
                                                if len(to_ext) > natmklst[n][natmk]:
                                                    natmklst[n][natmk] = len(to_ext)
                                            else:
                                                values.append(np.nan)
                                                if qn not in partial_states:
                                                    partial_states.append(qn)  # we only add it the first time it is not there
                                                if atom_key not in missing_atmks:  # we only add it the first time it is not there
                                                    missing_atmks.append(atom_key)
                                    else:
                                        values = [data[qn]]
                                else:  # qn not in data.keys()
                                    missing_states.append(qn)
                                    ccdlog.debug("added {} to \"missing\"".format(qn))
                                    values = [np.nan]
                                column.extend(values)
                            if True in [bool(i) for i in [missing_states, partial_states, missing_atmks]]:  # something missing
                                if cnt != 0:
                                    ccdlog.debug("{} could not be obtained".format([stateslist[n], q]))
                                else:
                                    if len(missing_states) == stateslist[n]:  # all states missing
                                        if atomiclist[n]:
                                            missing.append([stateslist[n], atomiclist[n], q])  # e.g. [5,"C1-C3",q]
                                        else:
                                            missing.append([stateslist[n], q])  # e.g. [5,q]
                                    elif len(partial_states) == stateslist[n]:  # all states partial
                                            atmstr = ",".join(missing_atmks)
                                            missing.append([stateslist[n], atmstr, q])  # e.g. [5,"C1-C3",q]
                                    else:
                                        if atomiclist[n]:
                                            missing.extend([[k,atomiclist[n]] for k in missing_states])  # e.g. [["q_3","C1-C3"], ["q_5","C1-C3"]]
                                            missing.extend([[k,",".join(missing_atmks)] for k in partial_states])  # e.g. ["q_7","C2"], ["q_8","C2"]]
                                        else:
                                            missing.extend(missing_states)  # e.g. ["q_3","q_5"]
                                    ccdlog.warning("had to get values for {} from individual values, check dictionary".format(q))  # here
                        
                        #
                        start = len(column)  # starting point for this set of vals
                        if atomiclist[n]:
                            if natmklst[n] != tmp:
                                for s in range(nsl[n]):
                                    for nk, v in enumerate(natmklst[n]):
                                        if tmp[nk] != v:
                                            for col in columns[:-1]:  # all the previous ones
                                                col = pd.concat(
                                                        [col[:start+sum([natmklst[:nk]])+tmp[nk]],
                                                             pd.Series((v-tmp[nk])*[np.nan]),
                                                             col[start+sum([natmklst[:nk]])+tmp[nk]:]],
                                                        ignore_index=True)  #  padding other columns
                        #
                    else:  # stateslist[n] == False
#                        to_add = data[q] if q in data.keys() else np.nan  # del?
                        if q not in data.keys():
                            if cnt == 0:
                                missing.append(oldqlist[n])
                                column.append(item if item else np.nan)
                                ccdlog.debug("added {} to \"missing\"".format(oldqlist[n]))
                            else:
                                column.append(item if item else np.nan)
                                ccdlog.debug("{} could not be obtained".format(oldqlist[n]))
                        elif type(item) is dict:
                            if not list(item.keys())[0].isnumeric():  # only one state, atomic values
                                item = {str(shift): item}  # make it {s: atomicvalsdict}
                            ccdlog.warning("You did not specify states for {q}. This makes collection somewhat slower. Please consider using [max_state,{q}] or \"{q}_n\" with n being the desired state".format(q=q))
                            values = []
                            try:
                                if atomiclist[n]:
                                    for v in range(shift, len(item)+shift):
                                        for natmk,atomkey in enumerate(atomiclist[n].split(",")):
                                            to_ext = item[str(v)][atomkey]
                                            to_ext.extend((natmklst[n][natmk] - len(to_ext))*[np.nan])  # if fewer atoms than other jobs
                                            values.extend(to_ext)
                                            if len(to_ext) > natmklst[n][natmk]:
                                                    natmklst[n][natmk] = len(to_ext)
                                            if look_for_more_states:
                                                missing.append(["{}_{}".format(q,v+1),atomiclist[n]])  # try to parse one more state
                                else:
                                    values = [item[str(v)] for v in range(shift, len(item)+1)]  # all in order
                                    if look_for_more_states:
                                                missing.append("{}_{}".format(q,len(values)+1))  # try to parse one more state
                                ccdlog.info("states {shift}-{n} in {q}".format(shift=shift,n=len(item) - 1 + shift, q=q))
                            except KeyError:  # not ordered (1,3,4,..) or atmk missing
                                max_state = max([int(i) for i in  item.keys()])
                                missing_states, partial_states, missing_atms = [], [], []
                                if atomiclist[n]:
                                    for v in range(shift, max_state+1):
                                        if str(v) in item.keys():
                                            for atomkey in atomiclist[n].split(","):
                                                if atomkey in item[str(v)].keys():
                                                    to_ext = item[str(v)][atomkey]
                                                    to_ext.extend((len(to_ext) - len(to_ext))*[np.nan])  # if fewer atoms than other job
                                                    values.extend(to_ext)
                                                    if len(to_ext) > natmklst[n][natmk]:  # if more atoms than other jobs
                                                        natmklst[n][natmk] = len(to_ext)
                                            else:
                                                values.append(natmklst[n][natmk]*[np.nan])
                                                if v not in partial_states:
                                                    partial_states.append(v)
                                                if atomkey not in missing_atms:
                                                    missing_atms.append(atomkey)
                                        else:
                                            missing_states.append(v)
                                    if look_for_more_states:                                        
                                        missing.extend([["{}_{}".format(q,k),atomiclist[n]] for k in missing_states])  # e.g. [["q_3","C1-C3"], ["q_5","C1-C3"]]
                                        missing.extend([["{}_{}".format(q,k),",".join(missing_atmks)] for k in partial_states])  # e.g. ["q_7","C2"], ["q_8","C2"]]
                                else:
                                    values = [item[str(v)] if str(v) in item.keys() else np.nan for v in range(shift, max_state+1)]
                                    if np.nan in values and look_for_more_states:
                                        missing.extend([n+shift for n,k in enumerate(values) if k == np.nan])
                                ccdlog.warning("some state or atomkey missing in {q}, using np.nan for it".format(q=q))
                            start = len(column)  # starting point for this set of vals
                            if atomiclist[n]:
                                if natmklst[n] != tmp:
                                    for s in range(nsl[n]):
                                        for nk, v in enumerate(natmklst[n]):
                                            if tmp[nk] != v:
                                                for col in columns[:-1]:  # all the previous ones
                                                    col = pd.concat(
                                                            [col[:start+sum([natmklst[:nk]])+tmp[nk]],
                                                                 pd.Series((v-tmp[nk])*[np.nan]),
                                                                 col[start+sum([natmklst[:nk]])+tmp[nk]:]],
                                                            ignore_index=True)  #  padding other columns
                            else:
                                values.extend((nsl[n] - len(values))*[np.nan])  # padding to match other columns
                            column.extend(values)
                            
                            if len(item) > nsl[n]:  # more states than before
                                for col in columns[:-1]:  # all the previous ones
                                    col = pd.concat(
                                            [col[:start+nsl[n]*sum(natmklst[n])],
                                                 pd.Series((len(values)-nsl[n]*sum(natmklst[n]))*[np.nan]),
                                                 col[start+nsl[n]*sum(natmklst[n]):]],
                                            ignore_index=True)  #  padding other columns
                                nsl[n] = len(item) - 1 + shift
                                ccdlog.debug("Padded previous columns with np.nan and adjusted length of next ones")
                        else:  # not a dict
                            column.append(item if item else np.nan)
                            ccdlog.debug("added column")
            if len(missing) != 0 and cnt == 0:
                ccdlog.info("trying to recalculate missing quantities")
                complex_quantities(path=path, qlist=missing, ex_qs = ex_qs,
                                   reqs=reqsdict, ext=ext, ignore=ignore, 
                                   parser_file=parser_file, parser=parserdict,
                                   parserfuncs=parserfuncs,
                                   parser_args=parserargsdict, 
                                   parser_kwargs=parserkwargsdict, check_input=False,
                                   funcdict=funcdict, to_console=to_console,
                                   to_log=to_log, logname=logname, printlevel=printlevel)  # calculated/reparsed
                cnt += 1
            else:
                cnt = 2  # either available or already tried calculating/reparsing
        columns.append(pd.Series(column))
    xp_qlist = []  # expanded qlist
    for n,q in enumerate(qlist):  # Could not manage with list comprehension
        if nsl[n] is not False:
            shift = 1 if q in ex_qs else 0  # 1 if only excited state
            if atomiclist[n]:
                myint = lambda x: int(x) if x else 1
                if nsl[n] - shift == 0:  # only one state
                    xp_qlist += ["{}_{}{}".format(q,
                                 "".join([i for i in atmk.split("-")[0] if i.isalpha()]),
                                 vl + myint("".join([i for i in atmk.split("-")[0] if i.isnumeric()]))) \
                                    for natmk,atmk in enumerate(atomiclist[n].split(",")) \
                                    for vl in range(natmklst[n][natmk])]
                else:
                    xp_qlist += ["{}_{}_{}{}".format(q, s+shift,
                                 "".join([i for i in atmk.split("-")[0] if i.isalpha()]),
                                 vl + myint("".join([i for i in atmk.split("-")[0] if i.isnumeric()]))) \
                                    for s in range(nsl[n] + 1 - shift) \
                                    for natmk,atmk in enumerate(atomiclist[n].split(",")) \
                                    for vl in range(natmklst[n][natmk])]
            else:
                if nsl[n] - shift == 0:  # only one state
                    xp_qlist.append(q)
                else:
                    xp_qlist += ["{}_{}".format(q, s+shift) for s in range(nsl[n] + 1 - shift)]
        else:
            xp_qlist.append(q)
    rows = levelnames + xp_qlist 
    ccdlog.debug("obtained rows")
    df = pd.concat(columns, axis=1)  # every column is a j 
    ccdlog.debug("concatenated")
    try:
        df.index = rows # index
    except ValueError:  # length mismatch (e.g. "ex_en_0")
        ccdlog.critical("Some issue in naming your columns!\n columns are {}".format(rows))
    ccdlog.info("done, return DataFrame")
    return df.T  # Transpose: now every column is either a level spec(A,B,basis,calc) or a property, every row is a job


def loopthrough(funcdict,joblist):
    """
    Parameters
    ----------
    funcdict : dict
        {func: {"args":[arglist],"kw": kwarg,...}, func2: {"args":[arglist],"kw": kwarg,...}, ...}
    joblist : list/str
        list of job tuples. or file to obtain it from

    Does
    -------
    Performs all the functions in funcdict for all jobs (e.g. complex_quantity, molden file)

    """    
    joblist = ut.deal_with_type(joblist, condition=str, to=get_joblist)
    for job in joblist:  # every job is a tuple
        for func,kwargs in funcdict.items():
            ckwargs = kwargs.copy()
            if "args" in ckwargs.keys():
                args = ckwargs["args"]
                del ckwargs["args"]
            else:
                args = []
            ckwargs.update({"path": os.path.join(*job)})  # we join the tuple elements to obtain the path
            func(*args,**ckwargs)

def correct_values(fplist, todel=[], toprocess={}):
    """
    """
    for fp in fplist:
        if not os.path.isfile(fp):
            ccdlog.warning("{} does not exist or is not a file".format(fp))
            continue
        d = ut.load_js(fp, cached=False)
        for td in todel:
            if td in d.keys():
                del d[td]
        for key, func in toprocess.items():
            if key in d.keys():
                try:
                    d[key] = func(d[key])
                except:
                    ccdlog.error("Could not use function {} on {}, {}".format(func.__name__, fp, key))
        ut.dump_js(d, fp)