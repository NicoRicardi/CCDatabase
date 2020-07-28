#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:51:40 2020

@author: nico
"""
import json as js
import itertools as ittl
import datetime as dt
import CCParser as ccp
import glob as gl
import numpy as np
import pandas as pd

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
        if len(data[q])<times:
            return False
        else:
            return True
        
def find_and_parse(path, ext="*.out", ignore="slurm*", parser=None, parser_args=None): 
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
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
    parser_args: list/tuple
        provide a list of the positional arguments for your parser
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
    parser, parser_args, parser_kwargs = (ccp.Parser,file,dict(software="qchem",to_file=True, to_console=False, to_json=True)  if not parser else (parser, *[file]+parser_args,**parser_kwargs)
    data = parser(parser_args,parser_kwargs)
    return data

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
        ID = max([int(i) for i in histdict.keys()] + 1  # get ID and add 1
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
                        bname = ".".join(splt[:-1])
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
    with open(histfile,"w") as f:  # updating histfile
        json.dump(histdict,f)

    with open(logfile, "w") as f:  # updating logfile
        for line in lines:
                f.write(line)
    
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

        
def raw_quantities(path="", qlist="variables.json", ext="*.out", ignore="slurm*", parser_file="CCParser.json", parser=None, parser_args=None, **kwargs):
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
    parser: None/function
        parser to call. If None, calls ccp.Parser
        If Parser != None, pass as the parser's *args as parser_args, then the parser's **kwargs
    parser_args: list/tuple
        provide a list of the positional arguments for your parser
    
    Does
    ----
    checks if the quantities are available in all "jobs". If necessary, it parses.
    Then prints what is missing
    
    Returns
    -------
    list
        the quantities missing in path
    """        
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            qlist = load_js(qlist)["raw_quantities"]
        else:  # actually a single quantity
            qlist = [qlist]
    if path == "":
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    missing = []
    jsfp = os.path.join(path, parser_file)  # path of raw quantities json file
    reparsed = {path: False}  # dict for current calculation and if needed iso, iso_g, sup, ...
    if not os.path.exists(jsfp):
        find_and_parse(path, ext=ext, ignore=ignore, parser=parser, parser_args=parser_args, **kwargs)
        reparsed [path] = True
    data = {}  # dict for current calculation and if needed iso, iso_g, sup, ...
    data[path] = load_js(jsfp)  # either already existed or we just created with find_and_parse
    for n,q in enumerate(qlist):
        times = 1
        if type(q) in [list,tuple] and len(q)==2:
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
                file = splt[0]
                type_ = file.split("."[-1])
                path_tmp = path
            else:  # cases 2,3,5
                if len(splt)==3:  # case 5
                    stdfile = False  # actually not in a standard json file
                    file = splt[1]
                    type_ = file.split("."[-1])
                else:  # cases 2,3
                    subdirs = gl.glob(os.path.join(path,"*",""))  # subdirectories of path (e.g. MP2_A for F&T)
                    paraldirs = gl.glob(os.path.join(os.path.split(path)[0],"*",""))  # directories in the parent folder (e.g. iso if path="emb")
                    if fol in subdirs and not in paraldirs:  # it is a subdir (case 2)
                        path_tmp = os.path.join(path,fol)
                    elif fol not in subdirs and in paraldirs:  # it is a parallel dir (case 3)
                        path_tmp = os.path.join(os.path.split(path)[0],fol)  
                    else:
                        raise ValueError("""the location of your quantity {} cannot be understood.
                                 Most likely it is either not present or double""".format(qlist[n]))
            q = splt[-1]
            if not stdfile:
                if type == "json":  # non-standard json file
                    fp = os.path.join(path_tmp,file)
                    if os.path.exists(fp):
                        if fp not in data.keys():
                            data[fp] = load_js(fp)
                        if q not in data[fp]:
                            missing.append(qlist[n])  # original q, not split
                        continue  # no need to run all other ifs
                    else:
                        raise FileNotFoundError("could not find {}".format(fp))
                elif type_ == "xlsx":  # excel file
                    df = pd.read_excel(file)
                    if q in not df.columns:  # times not available yet
                        missing.append(qlist[n])  # original q, not split
                    continue  # no need to run all other ifs
                elif type_ == "csv":  # csv file
                    df = pd.read_csv(file)
                    if q in not df.columns:
                        missing.append(qlist[n])  # original q, not split
                    continue  # no need to run all other ifs
                else:
                    raise TypeError("This type of non-json file is not implemented yet")
            elif path_tmp not in data.keys():
                jsfp_tmp = os.path.join(path_tmp, parser_file)
                    if not os.path.exists(jsfp_tmp):
                        find_and_parse(path_tmp, ext=ext, ignore=ignore, parser=parser, parser_args=parser_args, **kwargs)  # CCParser takes care of dumping json
                        reparsed[path_tmp] = True
                data[path_tmp] = load_js(jsfp_tmp)  # we are sure it exists    
                reparsed[path_tmp] = False
        else:  # case 1
            path_tmp = path
            jsfp_tmp = jsfp
        if not q_in_keys(data[path_tmp], q, times=times) and not reparsed[path_tmp]:  # quantity not in json, not reparsed yet
            find_and_parse(path_tmp, ext=ext, ignore=ignore, parser=parser, parser_args=parser_args, **kwargs)  # CCParser takes care of dumping json
            reparsed[path_tmp] = True
            data[path_tmp] = load_js(jsfp_tmp)  # let's read reparsed json
        if not q_in_keys(data[path_tmp], q, times=times) and reparsed[path_tmp]:
            missing.append(qlist[n])  # original q, not split
            print("{} missing {}".format(path_tmp,q))
    return missing   

def complex_quantities(path="", qlist="variables.json", reqs=None, ext="*.out", ignore="slurm*", parser=None, parser_args=None, **kwargs):
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
    if type(qlist) == list:
        if type(args[0]) != dict:
            raise TypeError("""If qlist is given as list, provide as args a dictionary:\n
                            {q1: [req1,req2,..], q2: [req1,req2,..],...}""")
        else:
            reqs = args[0]
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
            jsdict = load_js(qlist)
            qlist = jsdict["complex_quantities"]
            reqs = jsdict["requisites"]
        else:  # actually a single quantity
            qlist = [qlist]
            if len(args) == 1:
                reqs[qlist[0]] = args[0]
            else:
                reqs[qlist[0]] = args
    if path == "":
        path = os.getcwd()
    else:
        path = os.path.normpath(path)
    missing = []
    datafp = os.path.join(path, "data.json")  # path of complex quantities json file
    if os.path.exists(datafp):
        data = load_js(datafp)
    else:
        data = {}
    for q in qlist:
        times = 1
        if type(q) in [list,tuple] and len(q)==2:
            times,q = q[0],q[1]
        if not q_in_keys(data,q,times):
            miss = raw_quantities(path, qlist=reqs[q], ext=ext, ignore=ignore, parser=parser, parser_args=parser_args, **kwargs)  #Check requisite raw quantities are there (parses if necessary)
            if len(miss) == 0:
                qval = quantfuncs[q](path)
                data[q] = val
            else:
                print("cannot calculate {} because of missing raw quantity/ies".format(q))
                missing.append(q)
    dump_js(datafp,f)                        
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
            ckwargs.update("path": os.path.join(*job))  # we join the tuple elements to obtain the path
            f(**ckwargs)

def collect_data(joblist, levels=["A","B","basis","calc"], qlist="variables.json", *args):
    """

    Returns
    -------
    None.

    """
    if type(qlist) == list:
    if type(args[0]) != dict:
        raise TypeError("""If qlist is given as list, provide as args a dictionary:\n
                        {q1: [req1,req2,..], q2: [req1,req2,..],...}""")
    else:
        reqs = args[0]
    if type(qlist) == str:
        if re.match(".+\.json", qlist):  # actually a json file
        jsdict = load_js(qlist)
        qlist = jsdict["complex_quantities"]
        reqs = jsdict["requisites"]
    if len(levels) != len(joblist[0]):
        raise TypeError("Levels and joblist do not match!")
    rows = levels + qlist
    for j in joblist:
        path = os.path.join(*j)
        jsfp = os.path.join(path, "data.json")
        cnt = 0
        while cnt <= 1:  # counter to reparse at most once
            data = load_js(jsfp)
            column = levels
            for n,q in enumerate(qlist):
                to_add = data[q] if q in data.keys else np.nan
                column.append(to_add)
                if q not in data.keys() and cnt == 0:
                    missing.append(q)
            if len(missing) != 0 and cnt == 0:
                complex_quantities(path=path, qlist=missing, reqs)  # calculated/reparsed
                cnt += 1
            else:
                cnt = 2  # either available or already tried calculating/reparsing
        columns.append(column)
    df = pd.concat(columns, axis=1)
    df.index = rows
    dt = df.T
    return df

                   
        