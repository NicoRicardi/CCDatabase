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
import glob as gl
import numpy as np
import pandas as pd
import os
import re
import CCDatabase.utils as ut
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
    expr = ut.deal_with_type(expr,orig=str,to=lambda x: [x])  # if list turn to str
    if type(expr) != list:
        raise TypeError("Expr must be either list/tuple or string")
    ### Turning wildcards into regex
    expr = [i.replace("*",".+").replace("?",".") for i in expr]    
    ### Handling default startdir if none specified
    startdir = os.getcwd() if not startdir else startdir
    ### Handling default trash if none specified
    trash = os.path.join(os.getcwd(),"Trash") if not trash else trash
    ut.mkdif(trash)
    try:
        histdict = ut.load_js(histfile)  # read history file and look for highest ID
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
    ut.dump_js(histdict,histfile)          
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
    IDs = ut.deal_with_type(IDs,orig=str,to=lambda x: [x])  # if only 1 ID,turn to list
    histdict = ut.load_js(histfile)
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
    ut.dump_js(histdict,histfile)  # updating histfile

    with open(logfile, "w") as f:  # updating logfile
        for line in lines:
                f.write(line)   

### Quantities-related
                
def find_and_parse(path, ext="*.out", ignore="slurm*", parser=None,
                   parser_args=None, parser_kwargs=None, check_parsing=True): 
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
    check_parsing: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other quantities that check before for loops.
    Returns
    -------
    What your parser returns
        the parser container
    """
    ### Processing input
    # path
    path = ut.deal_with_type(path,condition=None,to=os.getcwd)
    path = os.path.abspath(path)
    if check_parsing:
        # parser
        try:  # this changes with python versions, hence try/except
            is_func = callable(parser)  # whether it is a function
        except:
            is_func = hasattr(parser,"__call__")  # whether it is a function
        if parser != None and not is_func:
            raise TypeError(""""Cannot process "parser", as it is none of: dict,list,func,None""")
        # parser_args
        if parser_args != None:
            if type(parser_args) == tuple:
                parser_args = list(parser_args)
            elif type(parser_args) != list:
                parser_args = [parser_args]
                print("""Your "parser_args" is neither None nor list/tuple, turning it into a list (parser_args = [parser_args])""")
        # parser_kwargs
        if parser_kwargs != None and type(parser_kwargs)!=dict:
            raise TypeError(""""Cannot process "parser_kwargs", as it is none of: dict,None""") 
    # Finding file
    files = [i for i in gl.glob(os.path.join(path,ext)) if i not in gl.glob(os.path.join(path,ignore))]  # e.g. all *.out which are not slurm*.out
    if len(files) == 0:
        raise FileNotFoundError("Either no matching file or inexistent path")
    elif len(files) == 1:
        file = files[0]
    else:
        lens = [len(f) for f in files]
        idx = lens.index(min(lens))        
        file = files[idx]
        print("Several matching files!! Using the shortest file ({})".format(file))
    #Let's get started and parse
    parser, parser_args, parser_kwargs = (ccp.Parser,[file],dict(software="qchem",to_file=True, to_console=False, to_json=True))  if parser == None else (parser, [file]+parser_args, parser_kwargs)
#    print("parser_args",parser_args)
#    print("parser_kwargs",parser_kwargs)
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
    data = ut.load_js(fname) 
    if "joblist" in data.keys():
        return(data["joblist"])
    else:
        ignore = data["ignore"] if "ignore" in data.keys() else []
        return [i for i in  ittl.product(*[data[i] for i in data["levels"]]) if i not in ignore]

def raw_quantities(path=None, qlist="variables.json", ext="*.out", ignore="slurm*",
                   parser_file="CCParser.json", parser=None, parser_args=None,
                   parser_kwargs=None, check_parsing=True):
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
    parser: dict{q: func}, func, None, False
        parser to call (for each q). If not dictionary, turns into dict. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
        NB. If the file is human-generated and cannot be obtained from a parser, pass False
    parser_args: dict[parser: list/tuple]/None
        provide a list of the positional arguments for your parser, for each parser.
        If one parser only you can pass the list alone, it is turned into a dict
    parser_kwargs: dict[parser: dict]/None
        provide a dict of the keyword arguments for your parser, for each parser.
        If one parser only you can pass the dict alone, it is turned into dict[parser: dict]
    check_parsing: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other quantities that check before for loops.
        
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
            qlist = ut.load_js(qlist)["raw_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    elif type(qlist) == list:
        if len(qlist) == 2 and type(qlist[0]) == int:
            qlist = [qlist]
    else:
        raise TypeError(""""could not process your qlist.
                        Please provide in one of the following ways:
                            - file.json (will read variable 'raw_quantities'
                            - "exc_energies"
                            - [2,"exc_energies"]
                            - ["exc_energies","osc_strength"]
                            - [[2,"exc_energies"],[2,"osc_strength"]]
                            """)
    nvalslist = [q[0] if type(q) in [tuple,list] else 1 for q in qlist]
    qlist = [q[1] if type(q) in [tuple,list] else q for q in qlist]
    # path
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    if check_parsing:
        ### parser, parser_args, parser_kwargs
        # parser
        if parser == None:
            parser = {q: parser for q in qlist}
        else:
            if type(parser) == tuple:
                parser = list(parser)
            if type(parser) == list:
                if len(parser) == 1:  # actually only one parser
                    try:
                        is_func = callable(parser[0])
                    except:
                        is_func = hasattr(parser[0],"__call__")
                    if not is_func:
                        raise ValueError("You gave parser as a list of 1 element, but this element is not a function")
                    parser = {q: parser[0] for q in qlist}
                elif len(parser) == len(qlist):
                    try:
                        are_funcs = np.array([callable(i)^(i==None) for i in parser]).all()
                    except:
                        are_funcs = np.array([hasattr(i,"__call__")^(i==None) for i in parser]).all()
                    if not are_funcs:
                        raise ValueError("You gave a list of parsers, at least one of those is neither a function nor None")
                    parser = {qlist[n]: parser[n] for n in range(len(qlist))}
                else:
                    raise ValueError("""You gave a list as "parser" but the length does not match that of "qlist" """)
            elif type(parser) == dict:
                if np.array([q not in parser.keys() for q in qlist]).any():
                    raise ValueError("""No parser for some of your quantities""")
                try:
                    are_funcs = np.array([callable(parser[q])^(parser[q]==None) for q in qlist]).all()
                except:
                    are_funcs = np.array([hasattr(parser[q],"__call__")^(parser[q]==None) for q in qlist]).all()
                if not are_funcs:
                    raise ValueError("You gave a dict of parsers, at least one of those is neither a function nor None")
                parser = parser
            else:
                try:  # this changes with python versions, hence try/except
                    is_func = callable(parser)  # whether it is a function
                except:
                    is_func = hasattr(parser,"__call__")  # whether it is a function
                if is_func:
                    parser = {q: parser for q in qlist}
                else:
                    raise TypeError(""""Cannot process "parser", as it is none of: dict,list,func,None""")
        # parser_args   
        if parser_args == None:
            parser_args = {q: parser_args for q in qlist}
        else:
            if type(parser_args) == tuple:  # turn to list
                parser_args = list(parser_args)
            if type(parser_args) == list:  # it was either list or tuple
                typelist = [type(i) for i in parser_args]
                if  not np.array([t in [tuple,list] for t in typelist]).all():  #  it's a single list of args
                    parser_args = {q:parser_args for q in qlist}  
                    print(""""You gave a single list "parser_args", it will be used for all quantities""")
                else:
                    if not np.array([t in [tuple,list,type(None)] for t in typelist]).all():
                        raise TypeError("""at least one element in "parser_args" has an invalid type""") 
                    if tuple in typelist:
                        for n,i in enumerate(parser_args):
                            if type(i) == tuple:  # turn to list
                                parser_args[n] = list(i)
                    if len(parser_args) == 1:  # actually only one list of args
                        parser_args = {q:parser_args[0] for q in qlist}  # I know it is in [list,None] (tuples have been turned to lists)
                    elif len(parser_args) == len(qlist):
                        parser_args = {qlist[n]: parser_args[n] for n in range(len(qlist))}  # I know all are in [list,None] (tuples have been turned to lists)
                    else:
                        raise ValueError("""You gave a list/tuple of lists/tuples as "parser_args" but the length does not match that of "qlist" """)
            elif type(parser_args) == dict:
                if np.array([q not in parser_args.keys() for q in qlist]).any():
                    raise ValueError("""No parser_args for some of your quantities """)
                typelist = [type(parser_args[q]) for q in qlist]  # allows unnecessary keys
                if tuple in typelist:
                        for k,v in parser_args.items():
                            if type(v) == tuple:  # turn to list
                                parser_args[k] = list(v)
                if np.array([t not in [list,type(None)] for t in typelist]).any():
                    raise TypeError("""at least one element in "parser_args" has an invalid type""") 
        # parser_kwargs    
        if parser_kwargs == None:
            parser_kwargs = {qlist[n]: parser_kwargs for n in range(len(qlist))}
        else:
            if type(parser_kwargs) == tuple:
                parser_kwargs = list(parser_kwargs)
            elif type(parser_kwargs) == list:
                if np.array([type(i) not in [dict,type(None)] for i in parser_kwargs]).any():
                    raise ValueError("""You can give a list of dictionaries with the same
                                     length of qlist, but your list's element are not only dictionaries or None""")
                if len(parser_kwargs) == len(qlist):
                    parser_kwargs = {qlist[n]: parser_kwargs[n] for n in range(len(qlist))}
                else:
                    raise ValueError("""You gave a list of dictionaries as "parser_kwargs" but the length does not match that of "qlist" """)
            elif type(parser_kwargs) == dict:
                if np.array([type(i) not in [dict,None] for i in parser_kwargs.values()]).any():  # It's not a dict of dicts/None
                    parser_kwargs = {qlist[n]: parser_kwargs for n in range(len(qlist))}
                elif np.array([q not in parser_kwargs.keys() for q in qlist]).any():
                    raise ValueError("""No parser_kwargs for some of your quantities""")
            else:
                raise ValueError(""""parser_kwargs" can be None, dict (to use for all quantities), or dict of dict. Your type is not recognised.""")

    ### Let's get started        
    missing = []
    reparsed = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    data = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    for n,q in enumerate(qlist):
        nvals = nvalslist[n]
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
                type_ = fname.split(".")[-1]
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
                paraldirs = gl.glob(os.path.join(ut.split_path(path)[0],"*",""))  # directories in the parent folder (e.g. iso if path="emb")
                if fol in subdirs and fol not in paraldirs:  # it is a subdir (case 2)
                    path_tmp = os.path.join(path,fol)
                elif fol not in subdirs and fol in paraldirs:  # it is a parallel dir (case 3)
                    path_tmp = os.path.join(ut.split_path(path)[0],fol)  
                else:
                    raise ValueError("""the location of your quantity {} cannot be understood.
                             Most likely it is either not present or double""".format(qlist[n]))
            q = splt[-1]
            filepath = os.path.join(path_tmp, fname)
            if not stdfile:
                if type_ == "xlsx":  # excel file
                    df = pd.read_excel(filepath)
                    if q not in df.columns:  # nvals not available yet
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
            data[path_tmp] = ut.load_js(filepath) if os.path.exists(filepath) else {}  # load if exists, else empty, will reparse
            if path_tmp not in reparsed.keys():
                reparsed[path_tmp] = {}
            reparsed[path_tmp][parser[qlist[n]]] = False  # full q, not split
        if not ut.rq_in_keys(data[path_tmp], q, nvals=nvals) and not reparsed[path_tmp][parser[qlist[n]]]:  # quantity not in json, not reparsed yet
            if parser[qlist[n]] != False:  # no parsing for human-generated quantities (e.g. correspondance)
                # Here parser should be {q1: parser1,..}, parser_args {q1:[arg1,..],..}, parser_kwargs {q:{kw1:arg1,..},..}
                find_and_parse(path_tmp, ext=ext, ignore=ignore, parser=parser[qlist[n]], parser_args=parser_args[qlist[n]], parser_kwargs=parser_kwargs[qlist[n]], check_parsing=False)  # CCParser takes care of dumping json
                data[path_tmp] = ut.load_js(filepath)  # let's read reparsed json
            reparsed[path_tmp][parser[qlist[n]]] = True  
        if not ut.rq_in_keys(data[path_tmp], q, nvals=nvals) and reparsed[path_tmp][parser[qlist[n]]]:
            missing.append(qlist[n])  # original q, not split
            print("{} missing {}".format(path_tmp,q))
    return missing   

def complex_quantities(path=None, qlist="variables.json", reqs=None, ext="*.out",
                       ignore="slurm*", parser_file="CCParser.json",parser=None,
                       parser_args=None, parser_kwargs=None, check_parsing=True,funcdict="ccp"):
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
    parser_file: str
        standard parser_file to pass to raw_quantities. Default is ccp
    parser: None/dict{functions}
        parser to call. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
    parser_args: None/dict{list/tuple}
        provide a list of the positional arguments for your parser
    parser_kwargs: None/dict{list/tuple}
        provide a list of the positional arguments for your parser
    check_parsing: bool
        whether parser,parser_args, parser_kwargs should be checked. Default is True.
        False is used by other quantities that check before for loops.
    funcdict: dict/str
        dictionary {"complex_quantity": function, ...}
        or one of these strings: "ccp","qcep-ccp","qcep"
        strings recognition is case insensitive and removes non-alphabetic characters
        
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
            jsdict = ut.load_js(qlist)
            qlist = jsdict["complex_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    elif type(qlist) == list:
        if len(qlist) == 2 and type(qlist[0]) == int:
            qlist = [qlist]
    else:
        raise TypeError(""""could not process your qlist.
                        Please provide in one of the following ways:
                            - file.json (will read variable 'raw_quantities'
                            - "exc_energies"
                            - [2,"exc_energies"]
                            - ["exc_energies","osc_strength"]
                            - [[2,"exc_energies"],[2,"osc_strength"]]
                            """)
                            
    stateslist = [q[0] if type(q) in [tuple,list] else False for q in qlist]
    qlist = [q[1] if type(q) in [tuple,list] else q for q in qlist]
    #funcdict
    if type(funcdict) == str:
        funcdict = (''.join([i for i in funcdict if i.isalpha()])).lower()  # removes non-aplhabetic
        if funcdict == "ccp":
            from CCDatabase.quantity_functions import ccp_funcs
            funcdict = ccp_funcs
        elif funcdict == "qcepccp":
            from CCDatabase.quantity_functions import qcep_ccp_funcs
            funcdict = qcep_ccp_funcs
        elif funcdict == "qcep":
            from CCDatabase.quantity_functions import qcep_funcs
            funcdict = qcep_funcs
    elif funcdict != dict:
        raise TypeError(""""funcdict" should be a dictionary or a string""")
    
    # reqs
    if reqs == None:
        raise ValueError("You must specify the required raw quantities for each complex quantity as reqs")
    elif type(reqs) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            if "jsdict" not in locals():
                jsdict = ut.load_js(qlist)    
            reqs = jsdict["requisites"]
        else:
            reqs = {q: [reqs] for q in qlist}
            print("You only gave one raw quantity as requisite. That can happen but is rare. Sure about it?")
    elif type(reqs) in [list,tuple]:
        reqs = {q: reqs for q in qlist}
    elif type(reqs) == dict:
        for k,v in reqs.items():
            reqs[k] = ut.deal_with_type(v,condition=str,to=lambda x: [x])
    else:
        raise TypeError("""reqs should be a dictionary {q1: [req1,req2,..], q2: [req1,req2,..],...} \n
                                                        otherwise a json to get it from (as jsdict["reqs"]) or a single quantity(str)""")
    assert (np.array([q in reqs.keys() for q in qlist]).all()), "Some quantity does not have its requisites specified" 
    # any raw which is a req of a complex
    all_raws = [r[1] if type(r) in [tuple,list] else r for r in set(ittl.chain.from_iterable(reqs.values()))]  
    if check_parsing:
        # parser
        if parser == None:  
            parserdict = {q: parser for q in all_raws}
        else:
            parser = ut.deal_with_type(parser,condition=tuple,to=list)
            if type(parser) == dict:
                if np.array([ i not in parser.keys() for i in all_raws]).any():
                    raise ValueError("No parser for some of your requisite raw quantities")
                try:
                    are_funcs = np.array([callable(parser[i])^(i==None) for i in all_raws]).all()
                except:
                    are_funcs = np.array([hasattr(parser[i],"__call__")^(i==None) for i in all_raws]).all()
                if not are_funcs:
                    raise TypeError("At least one of your parsers is neither a function nor None")  
                parserdict = parser
            elif type(parser) == list:
                raise TypeError(""""parser" cannot be a list. it can be None/func/dict """)
            else:
                try:  # this changes with python versions, hence try/except
                    is_func = callable(parser)  # whether it is a function
                except:
                    is_func = hasattr(parser,"__call__")  # whether it is a function
                if is_func:
                    parserdict = {q: parser for q in all_raws}
                else:
                    raise TypeError(""""Cannot process "parser", as it is none of: dict,list,func,None""")
        # parser_args
        if parser_args == None:
            parserargsdict = {q: parser_args for q in all_raws}
        elif type(parser_args) == dict:
            if np.array([i not in parser_args.keys() for i in all_raws]).any():
                raise ValueError("No parser_args for some of your requisite raw quantities")
            typelist = [type(parser_args[i]) for i in all_raws]
            if np.array([t not in [list,tuple,type(None)] for t in typelist]).any():
                raise TypeError("At least one of the elements in parser_args is none of list,tuple,None")
            parserargsdict = parser_args
        elif type(parser_args) in [list,tuple]:
            parser = ut.deal_with_type(parser,condition=tuple,to=list)
            parserargsdict = {q: parser_args for q in all_raws}
        else:
            raise ValueError("""Cannot process your "parser_args" """)
        # parser_kwargs
        if parser_kwargs == None:
            parserkwargsdict = {q: parser_kwargs for q in all_raws}
        elif type(parser_kwargs) == dict:
            if np.array([type(i) not in [dict,type(None)] for i in parser_kwargs.values()]).any():  # not a dict of dicts, using for all
                parserkwargsdict = {q: parser_kwargs for q in all_raws}
            else:
                if np.array([i not in parser_kwargs.keys() for i in all_raws]).any():
                    raise ValueError("No parser_kwargs for some of your requisite raw quantities")
                parserkwargsdict = parser_kwargs
        else:
            raise ValueError(""""parser_kwargs" can be None, dict (to use for all quantities), or dict of dict. Your type is not recognised.""")
    else:
        parserdict = parser
        parserargsdict = parser_args
        parserkwargsdict = parser_kwargs
    # path
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    ### Let's get started
    missing = []
    datafp = os.path.join(path, "data.json")  # path of complex quantities json file
    if os.path.exists(datafp):
        data = ut.load_js(datafp)
    else:
        data = {}
    for n,q in enumerate(qlist):
        if q not in funcdict.keys(): # dictionary of {q: func} that calculate complex quantities (reads specific jsons to get raw, processes, returns)
            par_q = "_".join(q.split("_")[:-1])  # parent q, e.g. "ex_en_11" => "ex_en"
            if par_q in funcdict.keys():
                func = funcdict[par_q]
                state_num = int(q.split("_")[-1])  # the state(only one) we want to look at
            else:
                raise ValueError("There is no function for complex quantity {}".format(q))
        else:
            func = funcdict[q]
            state_num = False
        states = stateslist[n]
        if not ut.cq_in_keys(data,q,states=states):
            miss = raw_quantities(path, qlist=reqs[q], ext=ext, ignore=ignore,
                                      parser_file=parser_file,
                                      parser={r:parserdict[r] for r in reqs[q]},
                                      parser_args={r:parserargsdict[r] for r in reqs[q]},
                                      parser_kwargs={r:parserkwargsdict[r] for r in reqs[q]},check_parsing=False)  #Check requisite raw quantities are there (parses if necessary)
            if len(miss) == 0:
                if states:
                    vals={}
                    for s in range(1,states+1):  # python counting => fortran/human counting
                        try:
                            qval = func(path=path,n=s)  
                            vals[s] = qval  # nb will become a str when dumped
                            data["{}_{}".format(q,s)] = qval  # e.g. ex_en_1 : val1
                        except:
                            print("Errors while calculating {}, state {}, in {}".format(q,s,path))
                            missing.append(q)
                    data[q] = vals  # e.g ex_en:{1:val1,2:val2}
                elif state_num:
                    try:
                        qval = func(path=path,n=state_num)  
                        data[q] = qval
                    except:
                        print("Errors while calculating {} in {}".format(q,path))
                        missing.append(q)
                else:  # no state specified, could be 1 or "as many as possible"
                    failed = False
                    s = 1
                    vals = {}
                    while not failed:  # trying to obtain as many vals as possible
                        try:
                            qval = func(path=path,n=s)  
                            vals[s] = qval  # nb will become a str when dumped
                            s += 1
                        except:  # max num of vals
                            failed = True 
                    if len(vals) == 0:
                        print("Errors while calculating {} in {}".format(q,path))
                        missing.append(q)
                    elif len(vals) == 1:
                        data[q] = vals[1]
                    else:
                        data[q] = vals
                        for s,v in vals.items():
                            data["{}_{}".format(q,s)] = v  # e.g. ex_en_1 : val1
            else:
                print("cannot calculate {} because of missing raw quantity/ies".format(q))
                missing.append(q)
    ut.dump_js(data,datafp)                        
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
    joblist = ut.deal_with_type(joblist,condition=str,to=get_joblist)
    for job in joblist:  # every job is a tuple
        for func,kwargs in funcdict.items():
            ckwargs = kwargs.copy()
            ckwargs.update({"path": os.path.join(*job)})  # we join the tuple elements to obtain the path
            func(**ckwargs)

def collect_data(joblist, levels=["A","B","basis","calc"], qlist="variables.json",
                 reqs=None, ext="*.out", ignore="slurm*", parser=None,
                 parser_file="CCParser.json",parser_args=None, 
                 parser_kwargs=None, check_parsing=True, funcdict="ccp"):
    """
    Note
    ----
    
    Parameters
    ----------
    joblist: list[tuples/list]/str
        Either list of jobs, or json file to retrieve
    levels: list[str]
        list of levels.Must match len of joblist items. Default is ["A","B","basis","calc"]
    qlist: list[str]/str
        list of complex quantities to get
    reqs: dict
        dictionary of requisite raw quantities for each complex quantity
    ext: str
        wildcard expression to match for fin_and_parse 
    ignore: str
        wildcard expression to ignore for fin_and_parse 
    parser: dict/func
        dictionary of parsers for each 
    funcdict: dict/str
        dictionary {"complex_quantity": function, ...}
        or one of these strings: "ccp","qcep-ccp","qcep"
        strings recognition is case insensitive and removes non-alphabetic characters

   Returns
    -------
    pd.DataFrame
        a DataFrame with ...
    """
    ### Processing user input
    # joblist
    joblist = ut.deal_with_type(joblist,condition=str,to=get_joblist)
    # qlist
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            jsdict = ut.load_js(qlist)
            qlist = jsdict["complex_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    elif type(qlist) == list:
        if len(qlist) == 2 and type(qlist[0]) == int:
            qlist = [qlist]
    else:
        raise TypeError(""""could not process your qlist.
                        Please provide in one of the following ways:
                            - file.json (will read variable 'raw_quantities'
                            - "exc_energies"
                            - [2,"exc_energies"]
                            - ["exc_energies","osc_strength"]
                            - [[2,"exc_energies"],[2,"osc_strength"]]
                            """)
                            
    stateslist = [q[0] if type(q) in [tuple,list] else False for q in qlist]
    qlist = [q[1] if type(q) in [tuple,list] else q for q in qlist]
    # reqs
    if reqs == None:
        raise ValueError("You must specify the required raw quantities for each complex quantity as reqs")
    if type(reqs) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            if "jsdict" not in locals():
                jsdict = ut.load_js(qlist)    
            reqs = jsdict["requisites"]
        else:
            reqs = [reqs]
            print("You only gave one raw quantity as requisite. That can happen but is rare. Sure about it?")
    elif type(reqs) in [list,tuple]:
        reqs = {q: reqs for q in qlist}
    elif type(reqs) == dict:
        for k,v in reqs.items():
            reqs[k] = ut.deal_with_type(v,condition=str,to=lambda x: [x])
    else:
        raise TypeError("""reqs should be a dictionary {q1: [req1,req2,..], q2: [req1,req2,..],...} \n
                                                        otherwise a json to get it from (as jsdict["reqs"]) or a single quantity(str)""")
    assert (np.array([q in reqs.keys() for q in qlist]).all()), "Some quantity does not have its requisites specified" 
    # any raw which is a req of a complex
    all_raws = [r[1] if type(r) in [tuple,list] else r for r in set(ittl.chain.from_iterable(reqs.values()))]
    if check_parsing:
        # parser
        if parser == None:  # parser==None
            parserdict = {q: parser for q in all_raws}
        else:
            parser = ut.deal_with_type(parser,condition=tuple,to=list)
            if type(parser) == dict:
                if np.array([ i not in parser.keys() for i in all_raws]).any():
                    raise ValueError("No parser for some of your requisite raw quantities")
                try:
                    are_funcs = np.array([callable(parser[i])^(i==None) for i in all_raws]).all()
                except:
                    are_funcs = np.array([hasattr(parser[i],"__call__")^(i==None) for i in all_raws]).all()
                if not are_funcs:
                    raise TypeError("At least one of your parsers is neither a function nor None")  
                parserdict = parser
            elif type(parser) == list:
                raise TypeError(""""parser" cannot be a list. it can be None/func/dict """)
            else:
                try:  # this changes with python versions, hence try/except
                    is_func = callable(parser)  # whether it is a function
                except:
                    is_func = hasattr(parser,"__call__")  # whether it is a function
                if is_func:
                    parserdict = {q: parser for q in all_raws}
                else:
                    raise TypeError(""""Cannot process "parser", as it is none of: dict,list,func,None""")
        # parser_args
        if parser_args == None:
            parserargsdict = {q: parser_args for q in all_raws}
        elif type(parser_args) == dict:
            if np.array([i not in parser_args.keys() for i in all_raws]).any():
                raise ValueError("No parser_args for some of your requisite raw quantities")
            typelist = [type(parser_args[i]) for i in all_raws]
            if np.array([t not in [list,tuple,type(None)] for t in typelist]).any():
                raise TypeError("At least one of the elements in parser_args is none of list,tuple,None")
            parserargsdict = parser_args
        elif type(parser_args) in [list,tuple]:
            parser = ut.deal_with_type(parser,condition=tuple,to=list)
            parserargsdict = {q: parser_args for q in all_raws}
        else:
            raise ValueError("""Cannot process your "parser_args" """)
        # parser_kwargs
        if parser_kwargs == None:
            parserkwargsdict = {q: parser_kwargs for q in all_raws}
        elif type(parser_kwargs) == dict:
            if np.array([type(i) not in [dict,type(None)] for i in parser_kwargs.values()]).any():  # not a dict of dicts, using for all
                parserkwargsdict = {q: parser_kwargs for q in all_raws}
            else:
                if np.array([i not in parser_kwargs.keys() for i in all_raws]).any():
                    raise ValueError("No parser_kwargs for some of your requisite raw quantities")
                parserkwargsdict = parser_kwargs
        else:
            raise ValueError(""""parser_kwargs" can be None, dict (to use for all quantities), or dict of dict. Your type is not recognised.""")
    # levels
    if len(levels) != len(joblist[0]):
        raise TypeError("Levels and joblist do not match!")
        
    ### Let's get started
    nsl = stateslist.copy()  # new stateslist
    columns = []
    for j in joblist:
        path = os.path.join(*j)
        jsfp = os.path.join(path, "data.json")
        cnt = 0
        while cnt <= 1:  # counter to reparse at most once
            data = ut.load_js(jsfp) if os.path.exists(jsfp) else {}
            column = list(j)  
            missing = []
            for n,q in enumerate(qlist): 
                if stateslist[n]:
                    try:
                        values = [data[q][str(v)] for v in range(1,stateslist[n]+1)]
                        # if q or str[v] gives key error goes to except
                        column.extend(values)
                    except:
                        miss_state = False
                        for s in range(1,stateslist[n]+1):
                            qn = "{}_{}".format(q,s)
                            to_add = data[qn] if qn in data.keys() else np.nan
                            miss_state = miss_state if qn in data.keys() else True
                            column.append(to_add)
                        if miss_state  and cnt == 0:
                            missing.append([stateslist[n],q])
                else:  # stateslist[n] == False
                    to_add = data[q] if q in data.keys() else np.nan
                    if q not in data.keys():
                        if cnt == 0:
                            missing.append(q)
                            column.append(to_add)
                    elif type(to_add) == dict:
                        try:
                            values = [to_add[str(v)] for v in range(1,len(to_add)+1)]  # all in order
                        except KeyError:  # not ordered (1,3,4,..)
                            max_state = max([int(i) for i in  to_add.keys()])
                            values = [to_add[str(v)] if str(v) in to_add.keys() else np.nan for v in range(1,max_state+1)]
                        start = len(column)  # starting point for this set of vals
                        column.extend(values)
                        column.extend((nsl[n]-len(values))*[np.nan])
                        if len(to_add)>nsl[n]:
                            for col in columns[:-1]:  # all the previous ones
                                col = pd.concat(
                                        [col[:start+nsl[n]], pd.Series((len(values)-nsl[n])*[np.nan]),col[start+nsl[n]:]],
                                        ignore_index=True)  #  insert np.nan to match all lengths
                            nsl[n] = len(values)
                    else:  # not a dict
                        column.append(to_add)
                    if q not in data.keys() and cnt == 0:
                        missing.append(q)
            if len(missing) != 0 and cnt == 0:
#                print("launching complex")
                complex_quantities(path=path, qlist=missing, reqs=reqs, ext=ext, 
                                   ignore=ignore, parser_file=parser_file,
                                   parser=parserdict, parser_args=parserargsdict,
                                   parser_kwargs=parserkwargsdict, check_parsing=False,
                                   funcdict=funcdict)  # calculated/reparsed
                cnt += 1
            else:
                cnt = 2  # either available or already tried calculating/reparsing
        columns.append(pd.Series(column))
    xp_qlist = []  # expanded qlist
    for n,q in enumerate(qlist):  # Could not manage with list comprehension
        if nsl[n]:
            xp_qlist += ["{}_{}".format(q,s+1) for s in range(nsl[n])]
        else:
            xp_qlist.append(q)
    rows = levels + xp_qlist 
    df = pd.concat(columns, axis=1)  # every column is a j 
    df.index = rows # index
    return df.T  # Transpose: now every column is either a level spec(A,B,basis,calc) or a property, every row is a job
#

