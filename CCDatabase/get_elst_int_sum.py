#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 23:57:43 2021

@author: nico
"""

import os
import numpy as np
import subprocess as sp
#import re
import CCDatabase.utils as ut
import warnings

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

def elst_int_sum_iso(file, jsfile="CCParser.json", with_ccp=True, linenumbers=True):
    mainfol = ut.split_path(file)[0] if ut.split_path(file)[0] else os.getcwd()
    # Get potentials
    if os.path.isfile(os.path.join(mainfol,"v_coul.txt")):
        v_j_file = os.path.join(mainfol,"v_coul.txt")
        expansion = False
    elif os.path.isfile(os.path.join(mainfol,"v_coulomb_ME.txt")):
        v_j_file = os.path.join(mainfol,"v_coulomb_ME.txt")
        expansion = "ME"
    elif os.path.isfile(os.path.join(mainfol,"v_coulomb_SE.txt")):
        v_j_file = os.path.join(mainfol,"v_coulomb_SE.txt")
        expansion = "SE"
    else:
        raise FileNotFoundError("Cannot find the Coulomb potential file")
    v_j = np.loadtxt(v_j_file)
    
    if os.path.isfile(os.path.join(mainfol,"v_nucA.txt")):
        v_a = np.loadtxt(os.path.join(mainfol,"v_nucA.txt"))
    elif os.path.isfile(os.path.join(mainfol,"v_nuc_1.txt")):
        v_a = np.loadtxt(os.path.join(mainfol,"v_nuc_1.txt"))
    else:
        raise FileNotFoundError("Cannot find file for v_A")
        
    if os.path.isfile(os.path.join(mainfol,"v_nucB.txt")):
        v_b = np.loadtxt(os.path.join(mainfol,"v_nucB.txt"))
    elif os.path.isfile(os.path.join(mainfol,"v_nuc_0.txt")):
        v_b = np.loadtxt(os.path.join(mainfol,"v_nuc_0.txt"))
    else:
        raise FileNotFoundError("Cannot find file for v_B")  
    
    # Get expansion
    jsdata = ut.load_js(os.path.join(mainfol,jsfile)) if os.path.isfile(os.path.join(mainfol,jsfile)) else {}
    if "fde_expansion" in jsdata.keys():
        if expansion:
            assert expansion == jsdata["fde_expansion"][0][0] or jsdata["fde_expansion"][0]  # linenumber or not
        else:
            expansion = jsdata["fde_expansion"][0]
            if type(expansion) == list:
                if not linenumbers:
                    warnings.warn("get_elst_int_sum: You are adding values without \"linenumbers\"\
                              in a json file which has them. This can lead to issues in reading data.\
                              Consider passing \"linenumbers=True\", dummy linenumbers will be added")
                expansion = expansion[0]
            elif linenumbers:
                warnings.warn("get_elst_int_sum: You are adding values with \"linenumbers\"\
                              in a json file which does not have them. This can lead to issues in reading data.\
                              Consider passing \"linenumbers=False\"")
    elif with_ccp:
        import CCParser as ccp
        parsed = ccp.Parser(file, to_json=True, json_file=jsfile, to_console=False,
                            overwrite_file=False, overwrite_vals=False)
#        V_NN = parsed.results.V_AB.get_last()  # check
        if expansion:
            assert expansion == parsed.results.fde_expansion[0] 
        else:
            expansion = parsed.results.fde_expansion[0]
    else:
        try:
            s = str(sp.check_output("grep -i expansion {}".format(file))).upper()
            expansion = "ME" if "ME" in s else "SE"
        except:
            raise FileNotFoundError("Could not determine expansion")
#        s = str(sp.check_output("grep \"Nuc_A <-> Nuc_B\" {}".format(file)))
#        V_NN = float(re.search("-?\d+\.\d+",s).group())
#        
#    if "V_AB" in jsdata.keys:
#        V_NN = jsdata["V_AB"][0]
#        if type(V_NN) == list:  # value and linenumber
#            if not linenumbers:
#                warnings.warn("get_elst_int_sum: You are adding values without \"linenumbers\"\
#                              in a json file which has them. This can lead to issues in reading data.\
#                              Consider passing \"linenumbers=True\", dummy linenumbers will be added")
#            V_NN = V_NN[0]
#        elif linenumbers:
#            warnings.warn("get_elst_int_sum: You are adding values with \"linenumbers\"\
#                              in a json file which does not have them. This can lead to issues in reading data.\
#                              Consider passing \"linenumbers=False\"")
#        
    # Get DMs
    if os.path.isfile(os.path.join(mainfol, "Densmat_A.txt")):
        dmf_A = os.path.join(mainfol, "Densmat_A.txt")
    elif os.path.isfile(os.path.join(mainfol, "frag_1_HF_{}.txt".format(expansion))):
        dmf_A = os.path.join(mainfol, "frag_1_HF_{}.txt".format(expansion))
    if os.path.isfile(os.path.join(mainfol, "Densmat_B.txt")):
        dmf_B = os.path.join(mainfol, "Densmat_B.txt")
    elif os.path.isfile(os.path.join(mainfol, "frag_0_HF_{}.txt".format(expansion))):
        dmf_B = os.path.join(mainfol, "frag_0_HF_{}.txt".format(expansion))
    dm_A = np.loadtxt(dmf_A, dtype=np.float64, skiprows=1 if has_header(file) else 0)
    dm_B = np.loadtxt(dmf_B, dtype=np.float64, skiprows=1 if has_header(file) else 0)
    lA, lB = dm_A.shape[0], dm_B.shape[0]
    if is_square(lA):
        nbasA = int(np.sqrt(lA))
        dm_A = 2*dm_A.reshape([nbasA, nbasA])  # NB supposes only alpha!!!
    elif is_square(lA/2):
        nbasA = int(np.sqrt(lA/2))
        dm_A = dm_A.reshape([2,nbasA, nbasA]).sum(axis=0)
    if is_square(lB):
        nbasB = int(np.sqrt(lB))
        dm_B = 2*dm_B.reshape([nbasB, nbasB])  # NB supposes only alpha!!!
    elif is_square(lB/2):
        nbasB = int(np.sqrt(lB/2))
        dm_B = dm_B.reshape([2, nbasB, nbasB]).sum(axis=0)
    v_j = v_j.reshape([nbasA, nbasA])
    v_b = v_b.reshape([nbasA, nbasA])
    v_a = v_a.reshape([nbasB, nbasB])
    J = np.einsum('ab,ba', v_j, dm_A)
    AnucB = np.einsum('ab,ba',v_b, dm_A) 
    BnucA= np.einsum('ab,ba', v_a, dm_B)
    if linenumbers:
        J = [J, -1]
        AnucB = [AnucB, -1]
        BnucA = [BnucA, -1]
    jsdata.update(dict(J_sum_iso=[J], AnucB_sum_iso=[AnucB], BnucA_sum_iso=[BnucA]))  # Using same structure as ccp
    ut.dump_js(jsdata, os.path.join(mainfol,jsfile))
    
    