#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:15:05 2020

@author: nico
"""
import os
import re
import CCDatabase.utils as ut
import numpy as np
import copy as cp
import cachetools
import itertools as ittl
import subprocess as sp
import glob as gl
from pyscf import gto
from pyscf.dft import numint
from pyscf.dft.numint import eval_ao, eval_rho, NumInt
from pyscf.dft.libxc import XC
from pyscf.dft import gen_grid
import dmtools.BasisSet as bset
import logging

from CCDatabase.utils import caches

def deal_with_array(arr, to="arr"):
    """
    Note
    ----
    Deals with array-like types
    
    Parameters
    ----------
    arr: list/arr/matrix
        the array-like object to process
    to: str/type
        the desired type.
    
    Returns
    -------
    array-like
        the array like in the requested format
    """
    if type(arr) not in [list, np.ndarray, np.matrix]:
        raise ValueError("It is neither a list, nor an array, nor a matrix")
    options = {"arr": np.array, "array": np.array, "np.array": np.array,
               "np.ndarray": np.array, np.array: np.array,
               "list": list, list: list,
               "matrix": np.matrix, "np.matrix": np.matrix}
    if to not in options.keys():
        raise ValueError("Use one of the options below (str, or function when possible): \n arr, array, np.array, np.ndarray, list")
    to = options[to]
    if type(arr) == to:
        return arr
    elif to != list:
        return to(arr)  # array or matrix
    else:  # to == list
        return arr.tolist()
    
def raw_to_complex(path=None, rawfile="CCParser.json", raw_key="", n=0, 
                   linenumbers=True, first_only=False, group_values=False,
                   arr_type="arr"):
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
    if not raw_key:
        raise ValueError("you didn't specify the raw quantity's key")
    if n < 0:
        raise ValueError("negative index. Most likely you didn't specify that a quantity is in ex_qs")
    if first_only and n > 0:
        raise ValueError("You want a value only, but are asking for n > 0")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile, cached=True)
    vals = raws[raw_key]
    if group_values and np.array([type(i) in [float,np.float32,np.float64] for i in vals]).all():
        vals = group_values(raws[raw_key])
    val = vals[n][0] if linenumbers else vals[n]
    if type(val) == str and re.match(".+npz", val):
        val = ut.vals_from_npz(os.path.join(path, val), raw_key)[n]
    if type(val) not in [bool, int, float, str]:
        val = deal_with_array(val, to=arr_type)
    return val

def group_values(vals):
    """
    Note
    ----
    if values appear several times with different accuracy, replaces all appearances with the 
    highest number of decimal digits.
    
    Parameters
    ----------
    vals: list
        the values
    
    Returns
    -------
    list
        the values with the updated number of digits
    """
    decimals = lambda x: len(np.format_float_positional(x).split(".")[-1])
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

def get_index_dict(s, atomlist):
    """
    Parameters
    ----------
    s: str
        atomstring
    atomlist: list
        list of atoms in your geometry
    
    Returns
    -------
    dict
        {atomkey: [indexes]}
    """
    l=s.split(",")
    idxs = {}
    for i in l:
        if "-" in i:
            splt = i.split("-")
            al1, al2 = re.search("[A-Za-z]+", splt[0]), re.search("[A-Za-z]+", splt[1])
            al1, al2 = al1.group() if al1 else "", al2.group() if al2 else ""
            assert al1 == al2, "Inconsistent value in atomstring"
            n1, n2 = re.search("[0-9]+", splt[0]).group(), re.search("[0-9]+", splt[1]).group()
            partial = [al1+str(j) for j in range(int(n1),int(n2)+1)]
            vals = []
            for p in partial:
                if p.isalnum():
                    al, n = re.search("[A-Za-z]+", p).group(), int(re.search("[0-9]+", p).group())
                    if al in ["A","a"]:
                        vals.append(n-1)
                    else:
                        try:
                            vals.append([n for n,j in enumerate(atomlist) if j == al][n-1])
                        except:
                            pass  # perhaps logging?
                elif p.isnumeric():
                    raise ValueError("Number in atomstring. Use A1,A2,etc")
                else:
                    raise ValueError("Could not process value in atomstring")
            idxs[i] = vals
        else:
            if i.isalpha():
                idxs[i] = [n for n,j in enumerate(atomlist) if j == i]
            elif i.isalnum():
                al, n = re.search("[A-Za-z]+", i).group(), int(re.search("[0-9]+", i).group())
                if al in ["A","a"]:
                    idxs[i] = [n-1]
                else:
                    try:
                        idxs[i] = [[n for n,j in enumerate(atomlist) if j == al][n-1]]
                    except:
                        idxs[i] = []
            elif i.isnumeric():
                raise ValueError("Number in atomstring. Use A1,A2,etc")
            else:
                raise ValueError("Could not process value in atomstring")
    return idxs

def raw_atomic(path=None, atomstring="", n=0, rawfile="CCParser.json", 
               raw_key="", first_only=True, frag=0, all_frag_avail=True,
               linenumbers=True, arr_type="arr"):
    """
    Note
    ----
    Returns an atomic quantity. Uses the geometry to deduce value index.
    
    Parameters
    ----------
    path: str
        the folder path 
    atomstring: str
        a string which determines which atoms to select.
        "O" => every oxygen
        "O1" => the first oxygen
        "O2-O4" => oxygen 2,3,4
        "A1", "A2-A4": atom 1, atom 2,3,4
    n: int
        the state (0=GS,1=ES1, ..)
    rawfile: str
        the name of the raw quantity json file
    raw_key: str
        the desired property's key
    first_only: bool
        if the property is only GS
    frag: int
        what fragment the counting should be in
    all_frag_avail: bool
        whether the property is available only for frag_0 or all
    linenumbers: bool
        whether raw quantities include the linenumber or not
    arr_typ: str/type
        desired type of array-like objects
    
    Returns
    -------
    obj
        the property. if array-like, in the desired format
    """
    if not raw_key:
        raise ValueError("you didn't specify the raw quantity's key")
    if n < 0:
        raise ValueError("negative index. Most likely you didn't specify that a quantity is in ex_qs")
    if first_only and n > 0:
        raise ValueError("You want a value only, but are asking for n > 0")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile, cached=True)
    if "frag_xyz" in raws.keys():
        geomkey = "frag_xyz"
    elif "xyz" in raws.keys():
        geomkey = "xyz"
        if frag:
            raise ValueError("fragment geometries not available")
    else:
        raise ValueError("no type of geometry (frag_xyz/xyz) available")
    geoms = raws[geomkey][0][0]
    if  re.match(".+npz", geoms):
        npzfile = os.path.join(path,geoms)
        geoms = ut.vals_from_npz(npzfile, geomkey)
    else:
        if geomkey == "frag_xyz":
            geoms = [np.array(geom, dtype="object") for geom in geoms]
        else:
            geoms = np.array(geoms, dtype="object")
    if geomkey == "frag_xyz": 
        frag_atoms = geoms[frag][:,0]
        all_atoms = np.vstack(*geoms)[:,0]
        shift = sum([len(geom) for geom in geoms[:frag]])  # N_atoms in fragments before frag
    else:
        frag_atoms = geoms[:,0]
        all_atoms = frag_atoms
    all_vals = raws[raw_key]
    atomlist = all_atoms if all_frag_avail else frag_atoms
    if len(all_vals)%len(atomlist) != 0:
        raise AssertionError("The total number of values available is not a multiple of the number of atoms!")
    idict = get_index_dict(atomstring, atomlist)
    valsdict = {}
    for name,idxs in idict.items():
        idxs = [idx + n*len(atomlist) + shift for idx in idxs]  # adjust for state and previous fragments
        try:
            vals = [all_vals[idx][0] if linenumbers else all_vals[idx] for idx in idxs]
        except:
            vals = []
        for nv, val in enumerate(vals):
            if  type(val) == str and re.match(".+npz", val):
                vals[nv] = ut.vals_from_npz(os.path.join(path, val), raw_key)[idxs[nv]]
            if type(vals[nv]) not in [bool, int, float]:
                vals[nv] = deal_with_array(vals[nv], to=arr_type)
        valsdict[name] = vals
    if not list(ittl.chain.from_iterable(valsdict.values())):
        raise BaseException("No item in atomstring returned a quantity value")
    to_return = valsdict.copy()    
    for k,v in valsdict.items():
        if not v:
            del to_return[k]
    return to_return


@cachetools.cached(cache=caches["non_fdet_terms"])
def get_non_fdet_terms(path=None, rawfile="CCParser.json", linenumbers=True):  # cached=True is omitted, uses default, we can because other funcs checked
    """
    """
    ccdlog = logging.getLogger("ccd")
    ut.setupLogger()
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    E = {}
    # AB
    if os.path.isdir(os.path.join(path, "AB_MP2")):
        raw = ut.load_js(os.path.join(path, "AB_MP2", rawfile))
        E["HF_AB"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
        E["E_2_AB"] =raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    elif os.path.isdir(os.path.join(path, "AB_HF")):
        raw = ut.load_js(os.path.join(path, "AB_HF", rawfile))
        E["HF_AB"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    else:
        ccdlog.error("Could not find AB_MP2 nor AB_HF")
    # A
    if os.path.isdir(os.path.join(path, "A_MP2")):
        raw = ut.load_js(os.path.join(path, "A_MP2", rawfile))
        E["HF_A"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
        E["E_2_A"] =raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    elif os.path.isdir(os.path.join(path, "A_HF")):
        raw = ut.load_js(os.path.join(path, "A_HF", rawfile))
        E["HF_A"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    else:
        ccdlog.error("Could not find A_MP2 nor A_HF")
   # A_gh 
    if os.path.isdir(os.path.join(path, "A_MP2_gh")):
        raw = ut.load_js(os.path.join(path, "A_MP2_gh", rawfile))
        E["HF_A_gh"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
        E["E_2_A_gh"] =raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    elif os.path.isdir(os.path.join(path, "A_HF_gh")):
        raw = ut.load_js(os.path.join(path, "A_HF_gh", rawfile))
        E["HF_A_gh"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    else:
        ccdlog.error("Could not find A_MP2_gh nor A_HF_gh")
    # B
    if os.path.isdir(os.path.join(path, "B_MP2")):
        raw = ut.load_js(os.path.join(path, "B_MP2", rawfile))
        E["HF_B"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
        E["E_2_B"] =raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    elif os.path.isdir(os.path.join(path, "B_HF")):
        raw = ut.load_js(os.path.join(path, "B_HF", rawfile))
        E["HF_B"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    else:
        ccdlog.error("Could not find B_MP2 nor B_HF")
    # B_gh
    if os.path.isdir(os.path.join(path, "B_MP2_gh")):
        raw = ut.load_js(os.path.join(path, "B_MP2_gh", rawfile))
        E["HF_B_gh"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
        E["E_2_B_gh"] =raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    elif os.path.isdir(os.path.join(path, "B_HF_gh")):
        raw = ut.load_js(os.path.join(path, "B_HF_gh", rawfile))
        E["HF_B_gh"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    else:
        ccdlog.error("Could not find B_MP2_gh nor B_HF_gh")
    return E


def get_ref_terms(path=None, rawfile="CCParser.json", linenumbers=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    return get_non_fdet_terms(ut.split_path(path)[0], rawfile=rawfile, linenumbers=linenumbers)

@cachetools.cached(cache=caches["non_fdet_int"])
def get_non_fdet_int(path=None, rawfile="CCParser.json", linenumbers=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    d = get_non_fdet_terms(path, rawfile=rawfile, linenumbers=linenumbers)
    E_int = {}
    E_int["HF"] = d["HF_AB"] - d["HF_A"] - d["HF_B"]
    try:
        E_int["HF_CP"] = d["HF_AB"] - d["HF_A_gh"] - d["HF_B_gh"]
    except KeyError:
        ccdlog.info("at least some HF-ghost value unavailable")
    try:
        E_int["MP"] = E_int["HF"] + d["E_2_AB"] - d["E_2_A"] - d["E_2_B"]
    except KeyError:
        ccdlog.info("at least some MP value unavailable")
    try:
        E_int["MP_CP"] = E_int["HF_CP"] + d["E_2_AB"]  - d["E_2_A_gh"] - d["E_2_B_gh"]
    except KeyError:
        ccdlog.info("at least some MP-ghost value unavailable")
    return E_int
    
def get_ref_int(path=None, n=0, rawfile="CCParser.json", linenumbers=True):
    """
    """
    if n != 0:
        raise ValueError("Only Ground-State energy")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    return get_non_fdet_int(ut.split_path(path)[0], rawfile=rawfile, linenumbers=linenumbers)

@cachetools.cached(cache=caches["find_B"])
def find_emb_B(path=None):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if os.path.isdir(os.path.join(path, "MP2_B")):
        bfol = os.path.join(path, "MP2_B")
    elif os.path.isdir(os.path.join(path, "HF_B")):
        bfol = os.path.join(path, "HF_B")
    else:
        cyfols = [i for i in os.listdir(path) if "cy" in i]
        n_iters = len(cyfols)
        bfol = os.path.join(path, "cy{}".format(n_iters - 1 - n_iters%2))
    return bfol

@cachetools.cached(cache=caches["find_A"])
def find_emb_A(path=None):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if os.path.isdir(os.path.join(path, "MP2_A")):
        afol = os.path.join(path, "MP2_A")
    elif os.path.isdir(os.path.join(path, "HF_A")):
        afol = os.path.join(path, "HF_A")
    else:
        cyfols = [i for i in os.listdir(path) if "cy" in i]
        n_iters = len(cyfols) - 1
        afol = os.path.join(path, "cy{}".format(n_iters - n_iters%2))
    return afol

def get_energy_B(path=None, rawfile="CCParser.json", linenumbers=True, find=False, MP=True):  # cached=True is omitted, uses default, we can because other funcs checked
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if find:
        path = os.path.join(path, find_emb_B(path=path))
    raw = ut.load_js(os.path.join(path, rawfile))
    d = {}
    cycles = raw["cycle_energies"][-1][0] if linenumbers else raw["cycle_energies"][-1]
    if type(cycles) == str and re.match(".+npz", cycles):
        cycles = ut.vals_from_npz(os.path.join(path, cycles), "cycle_energies")
    d["HF_B"] = cycles[-1][-1]
    if MP:
        d["E_2_B"]= raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    return d

def get_fdet_terms_A(path=None, rawfile="CCParser.json", linenumbers=True, find=False, MP=True):  # cached=True is omitted, uses default, we can because other funcs checked
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if find:
        path = os.path.join(path, find_emb_A(path=path))
    raw = ut.load_js(os.path.join(path, rawfile))
    d = {}
    d["J"] = raw["J_int"][0][0] if linenumbers else raw["J_int"][0]
    d["V_NN"] = raw["V_AB"][0][0] if linenumbers else raw["V_AB"][0]
    d["AnucB"] = raw["AnucB"][0][0] if linenumbers else raw["AnucB"][0]
    d["BnucA"] = raw["BnucA"][0][0] if linenumbers else raw["BnucA"][0]
    d["Exc_nad_upd"] = raw["Exc_nad"][1][0] if linenumbers else raw["Exc_nad"][1]
    d["Ts_nad_upd"] = raw["Ts_nad"][1][0] if linenumbers else raw["Ts_nad"][1]
    d["Exc_nad_ref"] = raw["Exc_nad"][0][0] if linenumbers else raw["Exc_nad"][0]
    d["Ts_nad_ref"] = raw["Ts_nad"][0][0] if linenumbers else raw["Ts_nad"][0]
    d["Delta_lin"] = raw["fde_delta_lin"][0][0] if linenumbers else raw["fde_delta_lin"][0]
    d["HF_A"] = raw["scf_energy"][-1][0] if linenumbers else raw["scf_energy"][-1]
    d["expansion"] = raw["fde_expansion"][-1][0] if linenumbers else raw["fde_expansion"][-1]
    if MP:
        d["E_2_A"] = raw["mp_correction"][-1][0] if linenumbers else raw["mp_correction"][-1]
    return d


@cachetools.cached(cache=caches["fdet_terms"])
def get_fdet_terms(path=None, rawfile="CCParser.json", linenumbers=True, MP=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    afol = os.path.join(path, find_emb_A(path=path))
    bfol = os.path.join(path, find_emb_B(path=path))
    d = get_fdet_terms_A(path=afol, rawfile=rawfile, linenumbers=linenumbers, MP=MP)
    d.update(**get_energy_B(path=bfol, rawfile=rawfile, linenumbers=linenumbers, MP=MP))
    return d

@cachetools.cached(cache=caches["fdet_grouped"])
def group_fdet_terms(path=None, rawfile="CCParser.json", linenumbers=True, MP=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    fdet = get_fdet_terms(path=path, rawfile=rawfile, linenumbers=linenumbers, MP=MP)
    ref = get_ref_terms(path=path, rawfile=rawfile, linenumbers=linenumbers)
    E = {}
    E["Delta_HF_A"] = fdet["HF_A"] - ref["HF_A"] if fdet["expansion"] == "ME" else \
    fdet["HF_A"] - ref["HF_A_gh"]
    E["Delta_HF_B"] = fdet["HF_B"] - ref["HF_B"] if fdet["expansion"] == "ME" else \
    fdet["HF_B"] - ref["HF_B_gh"]
    E["Delta_HF"] = E["Delta_HF_B"] + E["Delta_HF_A"]
    E["elst_int"] = fdet["V_NN"] + fdet["J"] + fdet["AnucB"] + fdet["BnucA"]
    E["nonel,lin"] = fdet["Ts_nad_ref"] + fdet["Exc_nad_ref"] + fdet["Delta_lin"]
    E["nonel"] = fdet["Ts_nad_ref"] + fdet["Exc_nad_ref"]
    if MP:
        E["Delta_E_2_A"] = fdet["E_2_A"] - ref["E_2_A"] if fdet["expansion"] == "ME" else \
        fdet["E_2_A"] - ref["E_2_A_gh"]
        E["Delta_E_2_B"] = fdet["E_2_B"] - ref["E_2_B"] if fdet["expansion"] == "ME" else \
        fdet["E_2_B"] - ref["E_2_B_gh"]
        E["Delta_E_2"] = E["Delta_E_2_A"] + E["Delta_E_2_B"]
    return E
    
def get_fdet_int(path=None, n=0, rawfile="CCParser.json", linenumbers=True, lin=True, MP=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise ValueError("Only Ground-State energy")
    E = group_fdet_terms(path=path, rawfile=rawfile, linenumbers=linenumbers, MP=MP)
    E_int = E["Delta_HF"] + E["elst_int"]
    E_int +=  E["nonel,lin"] if lin else E["nonel"]
    if MP:
        E_int += E["Delta_E_2"]
    return E_int

def elst_hf(path=None, rawfile="CCParser.json", linenumbers=True):  # cached=True is omitted, uses default, we can because other funcs checked
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    raw = ut.load_js(os.path.join(path, rawfile))
    J = raw["tot_coulomb"][-1][0] if linenumbers else raw["tot_coulomb"][-1]
    V_ne = raw["nuc_attr"][-1][0] if linenumbers else raw["nuc_attr"][-1]
    V_nn = raw["nuc_repu"][-1][0] if linenumbers else raw["nuc_repu"][-1]
    return  J + V_ne + V_nn

def elst_fdet(path=None, rawfile="CCParser.json", linenumbers=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    afol = os.path.join(path, find_emb_A(path=path))
    bfol = os.path.join(path, find_emb_B(path=path))
    elst_A = elst_hf(path=afol, rawfile=rawfile, linenumbers=linenumbers)
    elst_B = elst_hf(path=bfol, rawfile=rawfile, linenumbers=linenumbers)
    fdet = get_fdet_terms_A(path=afol, rawfile=rawfile, linenumbers=linenumbers)
    elst_int = fdet["J"] + fdet["V_NN"] + fdet["AnucB"] + fdet["BnucA"] 
    return elst_A + elst_B + elst_int

@cachetools.cached(cache=caches["elst_ref"])
def elst_ref(path=None, rawfile="CCParser.json", linenumbers=True):
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if os.path.isdir(os.path.join(path, "AB_MP2")):
        abfol = os.path.join(path, "AB_MP2")
    elif os.path.isdir(os.path.join(path, "AB_HF")):
        abfol = os.path.join(path, "AB_HF")
    else:
        raise FileNotFoundError("Could not find AB reference folder")
    return elst_hf(path=abfol, rawfile=rawfile, linenumbers=linenumbers)

def get_elst_ref(path=None, rawfile="CCParser.json", linenumbers=True):
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    mainfol = ut.split_path(path)[0]
    return elst_ref(path=mainfol, rawfile=rawfile, linenumbers=linenumbers)

def deduce_expansion(path=None):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    basename = ut.path_basename(path).upper()
    if "ME" in basename and "SE" not in basename:
        return "ME"
    elif "SE" in basename and "ME" not in basename:
        return "SE"
    try:
        grep = str(sp.check_output("grep -i expansion {}".format(os.path.join(path, "*", "*.out")), shell=True)).lower()
    except:
        try:
            grep = str(sp.check_output("grep -i expansion {}".format(os.path.join(path, "*", "*")), shell=True)).lower()
        except:
            raise FileNotFoundError("Could not determine expansion")
    if "me" in grep and "se" not in grep:
        return "ME"
    if "se" in grep and "me" not in grep:
        return "SE"
    else:
        raise FileNotFoundError("Could not determine expansion")

@cachetools.cached(cache=caches["elst_sum_iso"])        
def elst_sum_iso(path=None, rawfile="CCParser.json", linenumbers=True, expansion=None):  # cached=True is omitted, uses default, we can because other funcs checked
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    mono = expansion == "ME"
    if expansion is None:
        raise ValueError("You must provide the expansion")
    # A
    if os.path.isdir(os.path.join(path, "A_MP2" if mono else "A_MP2_gh")):
        afol = os.path.join(path, "A_MP2" if mono else "A_MP2_gh")
    elif os.path.isdir(os.path.join(path, "A_HF" if mono else "A_HF_gh")):
        afol = os.path.join(path, "A_HF" if mono else "A_HF_gh")
    else:
        raise FileNotFoundError("Could not find isolated A")
    # B
    if os.path.isdir(os.path.join(path, "B_MP2" if mono else "B_MP2_gh")):
        bfol = os.path.join(path, "B_MP2" if mono else "B_MP2_gh")
    elif os.path.isdir(os.path.join(path, "B_HF" if mono else "B_HF_gh")):
        bfol = os.path.join(path, "B_HF" if mono else "B_HF_gh")
    else:
        raise FileNotFoundError("Could not find isolated B")
    elst_A = elst_hf(path=afol, rawfile=rawfile, linenumbers=linenumbers)
    elst_B = elst_hf(path=bfol, rawfile=rawfile, linenumbers=linenumbers)

    intfol = False
    for exp in expansion.upper(), expansion.lower():
        suspects = ["FT*{}".format(exp), "FnT*{}".format(exp)]
        for sus in suspects:
            globbed = gl.glob(os.path.join(path, sus))
            if not globbed:
                continue
            for fol in globbed:
                if os.path.isdir(os.path.join(fol, "cy0")):
                    intfol = os.path.join(fol, "cy0")
                    break  # Following part breaks all loops.
            else:
                continue
            break
        else:
            continue
        break
    if not intfol:
        raise FileNotFoundError("Could not find electrostatic interaction for sum of isolated fragments.\
                                \nGenerally this is in the 0-th cycle of freeze and thaw, please insert in \
                                for A embedded in isolated B")
    raw = ut.load_js(os.path.join(intfol, rawfile))
    J =  raw["J_sum_iso"][0] if type(raw["J_sum_iso"][0]) != list else raw["J_sum_iso"][0][0]  # not trusting linenumbers on this one
    AnucB = raw["AnucB_sum_iso"][0] if type(raw["AnucB_sum_iso"][0]) != list else raw["AnucB_sum_iso"][0][0]
    BnucA = raw["BnucA_sum_iso"][0] if type(raw["BnucA_sum_iso"][0]) != list else raw["BnucA_sum_iso"][0][0]
    V_NN = raw["V_AB"][0][0] if linenumbers else raw["V_AB"][0]
    elst_int = J + AnucB + BnucA + V_NN
    return elst_A + elst_B + elst_int

def get_elst_int_sum_iso(path=None, rawfile="CCParser.json", linenumbers=True):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    expansion = deduce_expansion(path=path)
    return elst_sum_iso(path=ut.split_path(path)[0], rawfile=rawfile, linenumbers=linenumbers, expansion=expansion)

@cachetools.cached(cache=caches["elst_change_ref"])
def elst_change_ref(path=None, n=0, rawfile="CCParser.json", linenumbers=True):
    """
    """
    if n != 0:
        raise ValueError("Only Ground-State energy")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    return get_elst_ref(path=path, rawfile=rawfile, linenumbers=linenumbers) - get_elst_int_sum_iso(path=path, rawfile=rawfile, linenumbers=linenumbers)

def elst_change_fdet(path=None, n=0, rawfile="CCParser.json", linenumbers=True):
    """/home/nico/.local/lib/python3.6/site-packages
    """
    if n != 0:
        raise ValueError("Only Ground-State energy")
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    return elst_fdet(path=path, rawfile=rawfile, linenumbers=linenumbers) - get_elst_int_sum_iso(path=path, rawfile=rawfile, linenumbers=linenumbers)
 
def has_header(fname):
    """
    Parameters
    ----------
    fname: str
        the filename to check
        
    Returns
    -------
    bool
        whether the file has a header
    """
    first_line_splt = str(sp.check_output(["head", "-n 1", fname]))[2:].split()
    header = True if len(first_line_splt) > 1 else False
    return header

def is_square(integer):
    root = np.sqrt(integer)
    if int(root + 0.5) ** 2 == integer: 
        return True
    else:
        return False
    
def read_density(dmf, header=None, not_full=True):  # TODO docstring
    """
    Parameters
    ----------
    fname: str
        the DM filename
    head: bool/None
        whether the file has a header, with None(default) that is deduced
    alpha_only: bool/None
        whether the file has only one matrix(alpha), with None(default) that is deduced
    
    Return
    ------
    np.array(Nbas, Nbas)
        alpha and beta DM
    """
    if header == None:  
        header = has_header(dmf)
    dm = np.loadtxt(dmf, dtype=np.float64, skiprows=1 if header else 0)
    nl = dm.shape[0]
    if is_square(nl):
        nbas = int(np.sqrt(nl))
        dm = dm.reshape([nbas, nbas])
        if not_full:
             dm *= 2
    if is_square(nl/2):
        nbas = int(np.sqrt(nl/2))
        dm = dm.reshape([2, nbas, nbas]).sum(axis=0)
    return dm

def read_and_reorder(dmf, coords, bas):
    """
    """
    dm = read_density(dmf)
    try:
        corr_path = os.environ.get("BAS_CORRESP")
    except:
        raise FileNotFoundError("Could not determine \"BAS_CORRESP\" environment variable.\
                                Please set it in your .bashrc and place there basis correspondence files.")
    corr_file = os.path.join(corr_path, ut.path_basename(bas).replace(".nwchem",".corr"))
    d = ut.load_js(corr_file, cached=False)
    if d["source"] == "qchem" and d["destination"] == "pyscf":
        order = bset.get_order(corr_file)
    elif d["source"] == "destination" and d["qchem"] == "pyscf":
        order = bset.get_rev_order(corr_file)
    else:
        raise ValueError("{} to {} is not a qchem <=> pyscf reordering!".format(d["source"], d["destination"]))
    atomlist = [i[0].replace("X-","") for i in coords]
    sort_arr = bset.get_sort_arr(atomlist, order)
    return bset.reorder(dm, sort_arr)
    
def dm_on_grid(mol, dm, points):
    """
    """
    ao_mol = eval_ao(mol, points, deriv=0)
    rho = eval_rho(mol, ao_mol, dm, xctype="LDA")
    return rho

def read_key(rawdict, k, b_only=False):
    if k not in rawdict.keys():
        kb = k+"_B"
        elconfb, dmfb, coordsb, basB = rawdict[kb]
        dmb = read_and_reorder(dmfb, coordsb, basB)
        molb = gto.M(atom=coordsb, basis=ut.read_file(basB),
                     charge=elconfb[0], spin=elconfb[1] - 1)
        to_return = [molb, dmb]
        ka = k+"_A"
        if ka in rawdict.keys():
            elconfa, dmfa, coordsa, basA = rawdict[ka]
            mola = gto.M(atom=coordsa, basis=ut.read_file(basA),
                         charge=elconfa[0], spin=elconfa[1] - 1)
            if b_only:
                to_return.append(mola)
            else:
                dma = read_and_reorder(dmfa, coordsa, basA)
                to_return.extend([mola, dma])
        return to_return
    else:
        elconf, dmf, coords, bas = rawdict[k]
        dm = read_and_reorder(dmf, coords, bas)
        mol = gto.M(atom=coords, basis=ut.read_file(bas),
                    charge=elconf[0], spin=elconf[1] - 1)
        return mol, dm

def get_grid(mol, gridlevel=4, obj=False):
    """
    """
    grid = gen_grid.Grids(mol)
    grid.level = gridlevel
    grid.build()
    to_return = grid if obj else (grid.coords, grid.weights)
    return to_return

def key_to_density(rawdict, k, gridpoints=False, weights=False, b_only=False, expansion="ME"):
    tmp = read_key(rawdict, k, b_only=b_only)
    if len(tmp) == 2:
        mol, dm = tmp
        case = 1
    else:
        if len(tmp) == 4:
            molb, dmb, mola, dma = tmp
            case = 2
        if len(tmp) == 3:
            molb, dmb, mola = tmp
            case = 3
        mol = mola + molb if expansion == "ME" else mola
    if gridpoints is False:
        gridpoints, weights = get_grid(mol)
    if case == 1:
        d = dm_on_grid(mol, dm, gridpoints)
    elif case == 2:
        d = dm_on_grid(mola, dma, gridpoints) + dm_on_grid(molb, dmb, gridpoints)
    else:
        d = dm_on_grid(molb, dmb, gridpoints)
    return d, gridpoints, weights
    
def densities_on_gridpoints(path=None, n=0, k1="HF_FDET", k2="HF_ref",
            b_only=False, rawfile="DMfinder.json"):
    """
    """
    ccdlog = logging.getLogger("ccd")
    ut.setupLogger()
    first_dec = lambda x: min([n for n, i in enumerate(np.format_float_positional(x).split(".")[-1]) if i != "0"])
    check_int = lambda x: first_dec(round(x) - x)
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    expansion = deduce_expansion(path=path)
    raw = ut.load_js(os.path.join(path, rawfile))
    d1, gridpoints, weights = key_to_density(raw, k1, b_only=b_only, expansion=expansion)
    int1 = np.dot(weights, d1)
    ccdlog.info("{} integrates to {}".format(k1, int1))
    assert check_int(int1) >= 3, "Non-integer integration for {}".format(k1)
    d2, *_ = key_to_density(raw, k2, gridpoints=gridpoints, weights=weights, b_only=b_only, expansion=expansion)
    int2 = np.dot(weights, d2)
    ccdlog.info("{} integrates to {}".format(k2, int2))
    assert check_int(int2) >= 3, "Non-integer integration for {}".format(k2)
    return d1, d2, gridpoints, weights

def densdiff(path=None, n=0, k1="HF_FDET", k2="HF_ref", rawfile="DMfinder.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise NotImplementedError("Only GS so far!")
    d1, d2, gridpoints, weights = densities_on_gridpoints(path=path, n=n, k1=k1,
                                                          k2=k2, rawfile=rawfile)
    return 0.5*np.dot(weights, np.absolute(d2 - d1))

def M_value(path=None, n=0, k1="HF_FDET", k2="HF_ref", rawfile="DMfinder.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise NotImplementedError("Only GS so far!")
    d1, d2, gridpoints, weights = densities_on_gridpoints(path=path, n=n, k1=k1,
                                                          k2=k2, b_only=True,
                                                          rawfile=rawfile)
    d = d2 - d1
    pos = np.where(d > 0)
    d[pos] = 0
    return -np.dot(weights, d)


def calc_kernel(func_kw, dmAvar, dmAnvar, dA, dB, grid, mola):
    """
    """
    ni = NumInt()
    vxc_tot, fxc_tot = ni.eval_xc(func_kw, dA+dB, spin=0, deriv=2)[1:3]  # on grid
    vxc_a, fxc_a = ni.eval_xc(func_kw, dA, spin=0, deriv=2)[1:3]
    vxc_nad = (vxc_tot[0] - vxc_a[0], None, None, None)  # because of restricted
    fxc_nad = (fxc_tot[0] - fxc_a[0],) + (None,)*9  # because of restricted
    fxc = numint.nr_rks_fxc(ni, mola, grid, func_kw, dmAvar, dmAvar, 0, True,
                                dA, vxc_nad, fxc_nad)
    return np.einsum('ab,ba', fxc, dmAnvar - dmAvar)

def get_xc_code(type_, val):
    """
    """
    keys = XC.keys()
    poss = []
    if val in keys:
        poss.append(val)
    elif "LDA_{}_{}".format(type_, val) in keys:
        poss.append("LDA_{}_{}".format(type_, val))
    elif "GGA_{}_{}".format(type_, val) in keys:
        "GGA_{}_{}".format(type_, val)
    if len(poss) == 0:
        raise ValueError("Could not find this functional in libxc!!")
    elif len(poss) == 1:
        return poss[0]
    else:
        value = XC[poss[0]]
        if True in [XC[i] != value for i in poss]:
            raise ValueError("Ambiguous value for this functional in libxc!!")
        else:
            return poss[0]
        
def get_kernel_terms(path=None, n=0, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise NotImplementedError("Only GS so far!")
    raw = ut.load_js(os.path.join(path, dmfindfile))
    expansion = deduce_expansion(path=path)
    molb, dmb, mola, dma = read_key(raw, kvar, b_only=False)
    mol = mola + molb if expansion == "ME" else mola
    grid = get_grid(mol, obj=True)
    dA, dB = dm_on_grid(mola, dma, grid.coords), dm_on_grid(molb, dmb, grid.coords)
    molb_,dmb_nvar,mola_,dma_nvar = read_key(raw, knvar, b_only=False)
    fols = find_emb_A(path=path), find_emb_B(path=path)
    d = {"fde_Tfunc": "K", "fde_Xfunc": "X", "fde_Cfunc": "C", "fde_XCfunc": "XC"}
    kernels = {}
    for n in range(2):
        kernel = {}
        ccpdata = ut.load_js(os.path.join(fols[n], ccpfile))   
        kw = {v: ccpdata[k][-1][0].upper() for k, v in d.items() if k in ccpdata.keys()}
        if not kw:
            continue
        assert "K" in kw.keys(), "missing kinetic"
        assert ("XC" in kw.keys()) ^ ("X" in kw.keys() and "C" in kw.keys())
        kernel.update({"{}".format(k): calc_kernel(get_xc_code(k, v), [dma, dmb][n],
              [dma_nvar, dmb_nvar][n], [dA, dB][n], [dB, dA][n], grid, [mola,molb][n]) for k, v in kw.items()})
        kernels[["A", "B"][n]] = kernel
    if not kernels:
        raise BaseException("Somehow no kernel obtained!")
    ccdlog.info("kernelsi {}".format(kernels))
    return kernels

@cachetools.cached(cache=caches["kernel"])
def kernel_sep(path=None, n=0, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise NotImplementedError("Only GS so far!")
    kernels = get_kernel_terms(path=path, n=n, kvar=kvar, knvar=knvar, dmfindfile=dmfindfile, ccpfile=ccpfile)
    kernels = {k: sum(v.values()) for k,v in kernels.items()}
    return kernels 

def kernel_tot(path=None, n=0, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json"):
    """
    """
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    if n != 0:
        raise NotImplementedError("Only GS so far!")
    kernels = kernel_sep(path=path, n=n, kvar=kvar, knvar=knvar, dmfindfile=dmfindfile, ccpfile=ccpfile)
    return sum(kernels.values()) / len(kernels.values())

def dipoles(path=None, n=0, rawfile="CCParser.json", ex_en_kw="exc_energy_rel", hf=False, linenumbers=True):
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile, cached=True)
    rdips = raws["dip_moment"]
    if n:
        if ex_en_kw not in raws.keys():
            raise ValueError("No excited states in this rawfile!")
        n_ex = len(raws[ex_en_kw])
        if n > n_ex:
            raise ValueError("The state requested is not available")
        idx = - n_ex + n -1
        return rdips[idx][0] if linenumbers else rdips[idx]
    elif not n:
        if ex_en_kw in raws.keys():
            idx = len(rdips) - len(raws[ex_en_kw]) - 1
        if not hf:
            return rdips[idx][0] if linenumbers else rdips[idx]
        else:
            return rdips[idx - 1][0] if linenumbers else rdips[idx - 1]

def tot_dipoles(path=None, n=0, rawfile="CCParser.json", ex_en_kw="exc_energy_rel", hf=False, linenumbers=True):
    path = ut.deal_with_type(path, condition=None, to=os.getcwd)
    rawfile = os.path.join(path, rawfile)
    raws = ut.load_js(rawfile, cached=True)
    rtots = raws["tot_dip"]
    if n:
        if ex_en_kw not in raws.keys():
            raise ValueError("No excited states in this rawfile!")
        n_ex = len(raws[ex_en_kw])
        if n > n_ex:
            raise ValueError("The state requested is not available")
        idx = -  + n -1
        return rtots[idx][0] if linenumbers else rtots[idx]
    elif not n:
        if ex_en_kw in raws.keys():
            idx = len(rtots) - len(raws[ex_en_kw]) - 1
        if not hf:
            return rtots[idx][0] if linenumbers else rtots[idx]
        else:
            return rtots[idx - 1][0] if linenumbers else rtots[idx - 1]
    
"""
for quantity functions it is often useful to use some general function with several parameters,
define all of them but "path" and "n" into a lambda function, which is then saved as quantity function.
e.g. 
"ex_en": lambda path, n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="exc_energies_rel")
"""
ccp_funcs = {
        "ex_en": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="exc_energy_rel"),
        "osc": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="osc_str"),
        "SCF": lambda path, n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="scf_energy", first_only=True),
        "E_linFDET_HF": lambda path, n: get_fdet_int(path=path, n=n, rawfile="CCParser.json", lin=True, MP=False),
        "E_FDET_HF": lambda path, n: get_fdet_int(path=path, n=n, rawfile="CCParser.json", lin=False, MP=False),
        "E_linFDET_MP": lambda path, n: get_fdet_int(path=path, n=n, rawfile="CCParser.json", lin=True, MP=True),
        "E_FDET_MP": lambda path, n: get_fdet_int(path=path, n=n, rawfile="CCParser.json", lin=False, MP=True),
        "E_ref_HF": lambda path, n: get_ref_int(path=path, n=n, rawfile="CCParser.json")["HF"],
        "E_ref_HF_CP": lambda path, n: get_ref_int(path=path, n=n, rawfile="CCParser.json")["HF_CP"],
        "E_ref_MP": lambda path, n: get_ref_int(path=path, n=n, rawfile="CCParser.json")["MP"],
        "E_ref_MP_CP": lambda path, n: get_ref_int(path=path, n=n, rawfile="CCParser.json")["MP_CP"],
        "elst_change_ref": lambda path, n: elst_change_ref(path=path, n=n, rawfile="CCParser.json"),
        "elst_change_FDET": lambda path, n: elst_change_fdet(path=path, n=n, rawfile="CCParser.json"),
        "densdiff_FDET_ref": lambda path, n: densdiff(path=path, n=n, k1="HF_FDET", k2="HF_ref", rawfile="DMfinder.json"),
        "densdiff_iso_ref": lambda path, n: densdiff(path=path, n=n, k1="HF_iso", k2="HF_ref", rawfile="DMfinder.json"),
        "densdiff_iso_FDET": lambda path, n: densdiff(path=path, n=n, k1="HF_iso", k2="HF_FDET", rawfile="DMfinder.json"),
        "M_value": lambda path, n: M_value(path=path, n=n, k1="HF_FDET", k2="HF_ref", rawfile="DMfinder.json"),
        "kernel_tot": lambda path, n: kernel_tot(path=path, n=n, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json"),
        "kernel_A": lambda path, n: kernel_sep(path=path, n=n, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json")["A"],
        "kernel_B": lambda path, n: kernel_sep(path=path, n=n, kvar="HF_FDET", knvar="MP_FDET", dmfindfile="DMfinder.json", ccpfile="CCParser.json")["B"],
        "tot_dip": lambda path, n: tot_dipoles(path=path, n=n, rawfile="CCParser.json", ex_en_kw="exc_energy_rel", hf=False, linenumbers=True)}

qcep_ccp_funcs = {
        "ex_en": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="exc_energy"),
        "osc": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="CCParser.json", raw_key="osc_strength"),
        "SCF": lambda path, n: raw_to_complex(path=path, n=n, rawfile="CCParser.json", raw_key="scf_energy", first_only=True),
        "EFG_e":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="CCParser.json", raw_key="EFG_tensor_e", first_only=True, arr_type="arr").items()},
         "EFG_n":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="CCParser.json", raw_key="EFG_tensor_n", first_only=True, arr_type="arr").items()},
        "EFG_t":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="CCParser.json", raw_key="EFG_tensor_t", first_only=True, arr_type="arr").items()}
        }

qcep_funcs = {
        "ex_en": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="qcep.json", raw_key="exc_energy"),
        "osc": lambda path, n: raw_to_complex(path=path, n=n-1, rawfile="qcep.json", raw_key="osc_strength"),
        "SCF": lambda path, n: raw_to_complex(path=path, n=n, rawfile="qcep.json", raw_key="scf_energy", first_only=True),
        "EFG_e":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="qcep.json", raw_key="EFG_tensor_e", first_only=True, arr_type="arr").items()},
         "EFG_n":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="qcep.json", raw_key="EFG_tensor_n", first_only=True, arr_type="arr").items()},
        "EFG_t":lambda path, atomstring, n:\
        {k: np.linalg.eigvals(v).tolist() for k,v in\
         raw_atomic(path=path, atomstring=atomstring, n=n, rawfile="qcep.json", raw_key="EFG_tensor_t", first_only=True, arr_type="arr").items()}
        }

