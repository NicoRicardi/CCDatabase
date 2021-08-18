#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 23:57:43 2021

@author: nico
"""

import os
import numpy as np
import subprocess as sp
import CCDatabase.utils as ut

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

def get_elst_int_sum(file, jsfile="CCParser.json", with_ccp=True):
    mainfol = ut.split_path(file)[0] if ut.split_path(file)[0] else os.getcwd()
    # Get potentials
    if os.path.isfile(os.path.join(mainfol,"v_coul.txt")):
        v_j_file = "v_coul.txt"
        expansion = False
    elif os.path.isfile(os.path.join(mainfol,"v_coulomb_ME.txt")):
        v_j_file = "v_coulomb_ME.txt"
        expansion = "ME"
    elif os.path.isfile(os.path.join(mainfol,"v_coulomb_SE.txt")):
        v_j_file = "v_coulomb_SE.txt"
        expansion = "SE"
    else:
        raise FileNotFoundError("Cannot find the Coulomb potential file")
    v_j = np.loadtxt(v_j_file)
    
    if os.path.isfile(os.path.join(mainfol,"v_nucA.txt")):
        v_a = np.loadtxt(os.path.join(mainfol,"v_nucA.txt"))
    elif os.path.isfile(os.path.join(mainfol,"v_nuc1.txt")):
        v_a = np.loadtxt(os.path.join(mainfol,"v_nuc1.txt"))
    else:
        raise FileNotFoundError("Cannot find file for v_A")
        
    if os.path.isfile(os.path.join(mainfol,"v_nucB.txt")):
        v_b = np.loadtxt(os.path.join(mainfol,"v_nucB.txt"))
    elif os.path.isfile(os.path.join(mainfol,"v_nuc0.txt")):
        v_b = np.loadtxt(os.path.join(mainfol,"v_nuc0.txt"))
    else:
        raise FileNotFoundError("Cannot find file for v_B")  
    
    # Get V_NANB and expansion
    if os.path.isfile(os.path.join(mainfol,jsfile)):
        jsdata = ut.load_js(os.path.join(mainfol,jsfile))
    if "V_AB" in jsdata.keys:
        V_NN = jsdata["V_AB"][0]
        if type(V_NN) == list:  # value and linenumber
            V_NN = V_NN[0]
        if expansion:
            assert expansion == jsdata["fde_expansion"][0][0] or jsdata["fde_expansion"][0]  # linenumber or not
        elif "fde_expansion" in jsdata.keys():
            expansion = jsdata["fde_expansion"][0]
            if type(expansion) == list:
                expansion = expansion[0]
    elif with_ccp:
        import CCParser as ccp
        parsed = ccp.Parser(file, to_json=True, json_file=jsfile, overwrite=False, overwrite_vals=False)
        V_NN = parsed.results.V_AB.get_last()  # check
        if expansion:
            assert expansion == parsed.results.fde_expansion # check
        else:
            expansion = parsed.results.fde_expansion # check
    else:
        try:
            s = str(sp.check_output("grep expansion {}".format(file))).upper()
            expansion = "ME" if "ME" in s else "SE"
        except:
            expansion = "ME"  # qchem default
        s = str(sp.checkout("grep \"Nuc_A <-> Nuc_B\" {}".format(file)))
        V__NN = float(re.search("-?\d+\.\d+",s).group())
        
    # Get DMs
    if os.path.isfile(os.path.join(mainfol, "Densmat_A.txt")):
        dmf_A = os.path.join(mainfol, "Densmat_A.txt")
    elif os.path.isfile(os.path.join(mainfol, "frag_1_HF_{}.txt".format(expansion))):
        dmf_A = os.path.join(mainfol, "frag_1_HF_{}.txt".format(expansion))
    if os.path.isfile(os.path.join(mainfol, "Densmat_B.txt")):
        dmf_B = os.path.join(mainfol, "Densmat_B.txt")
    elif os.path.isfile(os.path.join(mainfol, "frag_0_HF_{}.txt".format(expansion))):
        dmf_B = os.path.join(mainfol, "frag_0_HF_{}.txt".format(expansion))
    dm_A = np.loadtxt(dmf_A, dtype=np.float64, skiprows=1 if has_header(fname) else 0)
    dm_B = np.loadtxt(dmf_B, dtype=np.float64, skiprows=1 if has_header(fname) else 0)
    lA, lB = dm_A.shape[0], dm_B.shape[0]
    if is_square(lA):
        nbas = int(np.sqrt(lA))
        dm_A = 2*dm_A.reshape([nbas, nbas])  # NB supposes only alpha!!!
    elif is_square(lA/2):
        nbas = int(np.sqrt(lA/2))
        dm_A = dm_A.reshape([2,nbas, nbas]).sum(axis=0)
    if is_square(lB):
        nbas = int(np.sqrt(lB))
        dm_B = 2*dm_B.reshape([nbas, nbas])  # NB supposes only alpha!!!
    elif is_square(lB/2):
        nbas = int(np.sqrt(lB/2))
        dm_B = dm_B.reshape([2,nbas, nbas]).sum(axis=0)
    J = np.trace(np.dot(dm_A, v_J))
    AnucB = np.trace(np.dot(dm_A, v_B))  
    BnucA= np.trace(np.dot(dm_B, v_A))  
    tot = J + AnucB + BnucA + V_NN
    
    