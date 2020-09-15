#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 17:05:26 2020

@author: nico
"""
import CCDatabase.CCDatabase as ccd
import CCDatabase.utils as ut
import itertools as ittl
import pandas as pd
import CCDatabase.QCEasyParser as qcep
import os

class DataCollector():
    """
    """
    def __init__(self,joblist=None, levelnames=["A","B","basis","calc"], 
                 parsername=None, parser_file=None, parserdict=None, 
                 parserargsdict=None, parserkwargsdict=None, reqsdict=None, 
                 funcdict=None, cqlist=None, ex_qs=None, ext="*.out", 
                 ignore="slurm*", to_console=True, to_log=False,
                 logname="CCDatabase.log", printlevel=20):
        self.specs(dict(joblist=joblist,levelnames=levelnames,parser_file=parser_file,
                        parsername=parsername,parserdict=parserdict,
                        parserargsdict=parserargsdict,parserkwargsdict=parserkwargsdict,
                        reqsdict=reqsdict,funcdict=funcdict,cqlist=cqlist,
                        ex_qs=ex_qs,ext=ext,ignore=ignore,to_console=to_console,
                        to_log=to_log,logname=logname,printlevel=printlevel))
        if joblist == None:
            raise ValueError("You must specify joblist!!")
        else:
            joblist = ut.deal_with_type(joblist,condition=str,to=ccd.get_joblist)
        self.joblist = joblist
        self.levelnames = levelnames
        if parsername == None:  # process then use other args
            parsername = (''.join([i for i in parsername if i.isalpha()])).lower()  # removes non-aplhabetic
            self.reqsdict = reqsdict
            self.funcdict = funcdict
            self.all_raws = [r[1] if type(r) in [tuple,list] else r for r in set(ittl.chain.from_iterable(reqsdict.values()))]
            self.parserdict = parserdict
            self.parser_file = parser_file
            self.parserargsdict = parserargsdict
            self.parserkwargsdict = parserkwargsdict
            # TODO: process
        else:   # deduce other args
            self.reqsdict = {cq: megadict[parsername]["reqs"][cq] for cq in cqlist}
            self.funcdict = megadict[parsername]["funcdict"]
            self.all_raws = [r[1] if type(r) in [tuple,list] else r for r in set(ittl.chain.from_iterable(reqsdict.values()))]
            self.parserdict = {rq : megadict[parsername]["parser"] for rq in self.all_raws}  #TODO: change
            self.parser_file = megadict[parsername]["parser_file"]
#            self.parserargsdict = {rq: megadict[parsername]["parserargs"][rq] for rq in self.all_raws}
            ### {rq: [{rq:hook,..}]} because parser should parse all_raws at once
            self.parserargsdict = {parsername: [{rq: megadict[parsername]["parser_args"][rq] for rq in self.all_raws}]}
            self.parserkwargsdict = {parsername: megadict[parsername]["parser_kwargs"]}
        self.cqlist = cqlist 
        self.ex_qs = ex_qs if ex_qs else [cq for cq in cqlist if cq in megadict["ex_qs"]]
        self.ext = ext
        self.ignore = ignore
        self.to_console = to_console
        self.to_log = to_log
        self.logname = logname
        self.printlevel = printlevel
    
    def from_file(fname):
        specs = ut.load_js(fname)
        return DataCollector(**specs)
    
    def save_specs(self,fname):
        specs = self.specs
        ut.dump_js(specs,fname)
        
    def collect(self):
        self.df = ccd.collect_data(self.joblist, levelnames=self.levelnames, qlist=self.cqlist,
                 ex_qs=self.ex_qs, reqs=self.reqsdict, ext=self.ext, ignore=self.ignore, parser=self.parserdict,
                 parser_file=self.parser_file,parser_args=self.parserargsdict, 
                 parser_kwargs=self.parserkwargsdict, check_input=False, funcdict=self.funcdict,
                 to_console=self.to_console, to_log=self.to_log, logname=self.logname, printlevel=self.printlevel)
        
    def write_df(self,filename=None, overwrite=False):
        if filename == None:
            raise ValueError("You must specify a filename")
        if not hasattr(self,"df"):
            raise BaseException("df has not been calculated yet!")
        ext = filename.split(".")[-1]
        writers = {"csv": lambda x: self.to_csv(x), "xlsx": lambda x: self.to_excel(x)}
        if ext not in writers.keys():
            raise NotImplementedError("This extension is not implemented yet. Why don't you do it, champ?")
        if os.path.exists(filename) and not overwrite:  # we need to modify the name
            splt = filename.split(".")
            bname = ".".join(splt[:-1])
            count = 2
            while os.path.exists(filename):  # avoid overwriting it
                filename = "{}_{}.{}".format(bname,count,ext)
                count += 1
        writers[ext](filename)  # writes
        
    def get_and_write(self,filename=None, overwrite=False):
        self.collect(self)
        self.write_df(filename=filename,overwrite=overwrite)
        
        
ccpdict = {"reqs": {"ex_en": "exc_energy_rel"}, 
           "funcdict": "ccp",
           "parser_file": "CCParser.json",
           "parser": "ccp",
           "parser_args": {None},
           "parser_kwargs": {"exc_energy_rel":None}
           }

qcepccpdict = {"reqs": {"ex_en": "exc_en"}, 
           "funcdict": "qcepccp",
           "parser_file": "CCParser.json",
           "parser": "qcep",
           "parser_args": {"exc_energy": qcep.hooks["exc_energies"], "SCF": qcep.hooks["scf_energy"]},
           "parser_kwargs": {"json_file": "CCParser.json", "overwrite_vals":True}
           }

qcepdict = {"reqs": {"ex_en": "exc_en"}, 
           "funcdict": "qcep",
           "parser_file": "qcep.json",
           "parser":"qcep",
           "parser_args": {"exc_energy": qcep.hooks["exc_energies"], "SCF": qcep.hooks["scf_energy"]},
           "parser_kwargs": {"json_file": "qcep.json", "overwrite_vals":True}
           }

megadict = {"ccp": ccpdict, "qcepccp": qcepccpdict, "qcep": qcepdict, "ex_qs": ["ex_en"]}