#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:58:40 2021

@author: nico
"""
import CCDatabase.utils as ut
from CCDatabase.CCDatabase import find_and_parse
import os
import glob as gl
from CCDatabase.utils import caches
import cachetools
from CCDatabase.quantity_functions import find_emb_A, find_emb_B, deduce_expansion


def find_basfile(path, levs=4):
    cnt = 0
    while cnt < levs:
        nwchems = gl.glob(os.path.join(path, "*.nwchem"))
        if len(nwchems) == 1:
            return nwchems[0]
        path = ut.split_path(path)[0]
        cnt += 1
    raise FileNotFoundError("Could not find a unique .nwchem file going back {} levels".format(levs))
        

def find_matrix(fname, wildcardlist=[], jsonfile="CCParser.json", prop_key="ext_DM"):
    """
    """
    path = ut.split_path(fname)[0]
    if not wildcardlist:
        wildcardlist = ["Densmat_SCF.txt","Densmat_MP"]  # TODO add other options
    elif type(wildcardlist) not in [tuple, list]:
        wildcardlist = [wildcardlist]
    json_filepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    old = ut.load_js(json_filepath) if os.path.isfile(json_filepath) else {}
    matching = [ut.split_path(fp)[1] for wc in wildcardlist for fp in gl.glob(os.path.join(path,wc))]
    old[prop_key] = matching  # always overwrites values because if file is no longer there, pointer is useless
    ut.dump_js(old, json_filepath)
    return matching
    
@cachetools.cached(cache=caches["locate_iso_dmfiles"])
def locate_iso_dmfiles(path=None, filename="Densmat_SCF.txt", expansion="ME",
                       prop_key="HF_iso_{}", jsonfile="DMfinder.json",
                       ccpfile="CCParser.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if expansion == "ME":
        if os.path.isfile(os.path.join(path, "A_MP2", filename)):
            afile = os.path.join(path, "A_MP2", filename)
            afol = os.path.join(path, "A_MP2")
        elif os.path.isfile(os.path.join(path, "A_HF", filename)):
            afile = os.path.join(path, "A_HF", filename)
            afol = os.path.join(path, "A_HF")
        else:
            afile = False
        if os.path.isfile(os.path.join(path, "B_MP2", filename)):
            bfile = os.path.join(path, "B_MP2", filename)
            bfol = os.path.join(path, "B_MP2")
        elif os.path.isfile(os.path.join(path, "B_HF", filename)):
            bfile = os.path.join(path, "B_HF", filename)
            bfol = os.path.join(path, "B_HF")
        else:
            bfile = False
    elif expansion == "SE":
        if os.path.isfile(os.path.join(path, "A_MP2_gh", filename)):
            afile = os.path.join(path, "A_MP2_gh", filename)
            afol = os.path.join(path, "A_MP2_gh")
        elif os.path.isfile(os.path.join(path, "A_HF", filename)):
            afile = os.path.join(path, "A_HF_gh", filename)
            afol = os.path.join(path, "A_HF_gh")
        else:
            afile = False
        if os.path.isfile(os.path.join(path, "B_MP2_gh", filename)):
            bfile = os.path.join(path, "B_MP2_gh", filename)
            bfol = os.path.join(path, "B_MP2_gh")
        elif os.path.isfile(os.path.join(path, "B_HF_gh", filename)):
            bfile = os.path.join(path, "B_HF_gh", filename)
            bfol = os.path.join(path, "B_HF_gh")
        else:
            bfile = False
    else:
        raise NotImplementedError("Unknown expansion!! Only ME and SE so far!")
    json_filepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    old = ut.load_js(json_filepath) if os.path.isfile(json_filepath) else {}
    new = {}
#    if "{}" in prop_key:
#        prop_key = prop_key.format(expansion)
    if afile:
        if not os.path.isfile(os.path.join(afol, ccpfile)):
            find_and_parse(afol)
        raw = ut.load_js(os.path.join(afol, ccpfile))
        coords = raw["frag_xyz"][-1][0]
        if type(coords) == str:
            coords = ut.vals_from_npz(os.path.join(afol, coords), "frag_xyz")[-1][0].tolist()
        coords = [[i[0].replace("@", "X-")]+i[1:] for  i in coords]
        basA = find_basfile(afol)
        new[prop_key+"_A"] = [afile, coords, basA]
    if bfile:
        if not os.path.isfile(os.path.join(bfol, ccpfile)):
            find_and_parse(afol)
        raw = ut.load_js(os.path.join(bfol, ccpfile))
        coords = raw["frag_xyz"][-1][0]
        if type(coords) == str:
            coords = ut.vals_from_npz(os.path.join(bfol, coords), "frag_xyz")[-1][0].tolist()
        coords = [[i[0].replace("@", "X-")]+i[1:] for  i in coords]
        basB = find_basfile(bfol)
        new[prop_key+"_B"] = [bfile, coords, basB]
    if (afile or bfile) and jsonfile:
        old.update(new)
        ut.dump_js(old, json_filepath)
    return new
    
def find_iso_dmfiles(fname, filename="Densmat_SCF.txt", prop_key="HF_iso",
                     jsonfile="DMfinder.json", ccpfile="CCParser.json"):
    path = ut.split_path(fname)[0]
    fol = ut.split_path(path)[0]  # parent folder
    jsonfilepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    expansion = deduce_expansion(path=path)
    return locate_iso_dmfiles(path=fol, filename="Densmat_SCF.txt", expansion=expansion,
                       prop_key=prop_key, jsonfile=jsonfilepath if jsonfile else "", ccpfile=ccpfile)
    
@cachetools.cached(cache=caches["locate_ref_dmfile"])
def locate_ref_dmfile(path=None, filename="Densmat_SCF.txt", prop_key="HF_ref",
                      jsonfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if os.path.isfile(os.path.join(path, "AB_MP2", filename)):
        abfile = os.path.join(path, "AB_MP2", filename)
        abfol =  os.path.join(path, "AB_MP2")
    elif os.path.isfile(os.path.join(path, "AB_HF", filename)):
        abfile = os.path.join(path, "AB_HF", filename)
        abfol = os.path.join(path, "AB_HF")
    else:
        abfile = False
    json_filepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    old = ut.load_js(json_filepath) if os.path.isfile(json_filepath) else {}
    new = {}
    if abfile:
        if not os.path.isfile(os.path.join(abfol, ccpfile)):
            find_and_parse(abfol)
        raw = ut.load_js(os.path.join(abfol, ccpfile))
        coords = raw["frag_xyz"][-1][0]
        if type(coords) == str:
            coords = ut.vals_from_npz(os.path.join(abfol, coords), "frag_xyz")[-1][0].tolist()
        bas = find_basfile(abfol)
        new[prop_key] = [abfile, coords, bas]
        if jsonfile:
            old.update(new)
            ut.dump_js(old, json_filepath)
    return new
        
def find_ref_dmfile(fname, filename="Densmat_SCF.txt", prop_key="HF_ref",
                    jsonfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.split_path(fname)[0]
    fol = ut.split_path(path)[0]  # parent folder
    jsonfilepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    return locate_iso_dmfiles(path=fol, filename="Densmat_SCF.txt",
                              prop_key=prop_key, jsonfile=jsonfilepath if jsonfile else "", ccpfile=ccpfile)
    
@cachetools.cached(cache=caches["locate_fdet_dmfiles"])
def find_fdet_dmfiles(fname, filename="Densmat_SCF.txt", prop_key="HF_FDET",
                      jsonfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.split_path(fname)[0]
    expansion = deduce_expansion(path=path)
    afol, bfol = find_emb_A(path=path), find_emb_B(path=path)
    if os.path.isfile(os.path.join(afol, filename)):
        afile = os.path.join(afol, filename)
    else:
        afile = False
    if os.path.isfile(os.path.join(bfol, filename)):
        bfile = os.path.join(bfol, filename)
    else:
        bfile = False
    json_filepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    old = ut.load_js(json_filepath) if os.path.isfile(json_filepath) else {}
    new = {}
    if afile:
        if not os.path.isfile(os.path.join(afol, ccpfile)):
            find_and_parse(afol)
        raw = ut.load_js(os.path.join(afol, ccpfile))
        coords = raw["frag_xyz"][-1][0]
        if type(coords) == str:
            npz = os.path.join(afol, coords)
            coords = ut.vals_from_npz(npz, "frag_xyz")[-1][0].tolist()
            if expansion == "SE":
                ghost = ut.vals_from_npz(npz, "frag_xyz")[-1][1]
                ghost[:,0] = "X-" + ghost[:,0]
                coords.extend(ghost.tolist())
        elif expansion == "SE":
            ghost = raw["frag_xyz"][-1][1]
            coords.extend([["X+"+i[0]]+i[1:] for  i in ghost])
        basA = find_basfile(afol)
        new[prop_key+"_A"] = [afile, coords, basA]
    if bfile:
        if not os.path.isfile(os.path.join(bfol, ccpfile)):
            find_and_parse(afol)
        raw = ut.load_js(os.path.join(bfol, ccpfile))
        coords = raw["frag_xyz"][-1][0]
        if type(coords) == str:
            npz = os.path.join(bfol, coords)
            coords = ut.vals_from_npz(npz, "frag_xyz")[-1][0].tolist()
            if expansion == "SE":
                ghost = ut.vals_from_npz(npz, "frag_xyz")[-1][1]
                ghost[:,0] = "X+" + ghost[:,0]
                coords.extend(ghost.tolist())
        elif expansion == "SE":
            ghost = raw["frag_xyz"][-1][1]
            coords.extend([["X+"+i[0]]+i[1:] for  i in ghost])
        basB = find_basfile(bfol)
        new[prop_key+"_B"] = [bfile, coords, basB]
    if (afile or bfile) and jsonfile:
        old.update(new)
        ut.dump_js(old, json_filepath)
    return new
    
def get_all(fname, filename="Densmat_SCF.txt", ref="HF_ref", iso="HF_iso", FDET="HF_FDET",
                      jsonfile="DMfinder.json", ccpfile="CCParser.json"):
    path = ut.split_path(fname)[0]
    json_filepath = jsonfile if ut.split_path(jsonfile)[0] else os.path.join(path, jsonfile)
    old = ut.load_js(json_filepath) if os.path.isfile(json_filepath) else {}
    old.update(find_fdet_dmfiles(fname, filename=filename, prop_key=FDET,
                                 jsonfile="", ccpfile=ccpfile))   # jsonfile="" because we only  write at the end
    old.update(find_ref_dmfile(fname, filename=filename, prop_key=ref,
                                 jsonfile="", ccpfile=ccpfile)) 
    old.update(find_iso_dmfiles(fname, filename=filename, prop_key=iso,
                                 jsonfile="", ccpfile=ccpfile)) 
    ut.dump_js(old, json_filepath)