"""Some useful tools to comunicate with other QC software."""

import os
import re
import json
import numpy as np
from CCParser.QChem import extract_floats, parse_symmetric_matrix

def parse_inline_vec(line, asarray=True):  # from CCP, just to avoid version issues
    """ Extracts a vector of the format
    '[ 1.000, 2.000, 3.000]' from the current line"""
    pattern = r"[+-]?\d+\.\d*"
    match = re.findall(pattern, line)
    if len(match) > 0:
        match = list(map(float, match))
        if asarray:
            return np.asarray(match)
        else:
            return match
        
def parse_molecule(n, readlin, join=False):
    """Parse the geometry of all fragements in one file.

    Parameters
    ----------
    n : int
        Line number of identifier
    readlin : list
        Readlines list object

    Returns
    -------
    geos : list
        List of all atom symbols and cartesian coordinates.
    """
#    sep = []  #TODO check with Cris
    frag = 0
    geos = [[]]
    # First find the limits
    for line in readlin[n+2:]:
        if "--" in line:
            frag += 1
            geos.append([])
        elif "$end" in line:
            break
        elif "read" in line:
            geos[frag] = "read"
        else:
            data = line.split()
            if len(data) < 4:
                continue
            else:
                geos[frag].append([data[0]] + list(map(float, data[1:])))
    if join:
        from itertools import chain
        geos = list(chain.from_iterable(geos))
    return geos

def parse_elconf(n, readlin):
    """Parse the geometry of all fragements in one file.

    Parameters
    ----------
    n : int
        Line number of identifier
    readlin : list
        Readlines list object

    Returns
    -------
    geos : list
        List of all atom symbols and cartesian coordinates.
    frag_ids : dict
        Dictionary with the indices of each fragment.
    """
    elconfs, seps = [], 0
    # First find the limits
    for line in readlin[n:]:
        splt = line.split()
        if len(splt) == 2:
            elconfs.append((splt[0],splt[1]))
        if "--" in line:
            seps += 1
    if seps:
        assert len(elconfs) - seps  == 1
    else:
        assert len(elconfs) == 1
    return elconfs

def parse_simple_matrix(n, readlin, stop_signals=None, asmatrix=False):
    """Parse a symmetric matrix printed columnwise

    Parameters
    ----------
    n : int
        Line number of identifier
    readlin : list
        Readlines list object
    asmatrix : bool
        Whether to return a numpy.matrix object or not

    Returns
    -------
    numpy.matrix
        Parsed AO matrix as numpy.matrix object if asmatrix=True
    list
        Parsed AO matrix as list of lists if asmatrix=False
    """
    matrix = []
    index_line = n+1
    if stop_signals is None:
        stop_signals = ["Gap", "=", "eV", "Convergence", "criterion"]
    for iline, line in enumerate(readlin[index_line:]):
        if any(stop in line for stop in stop_signals):
            break
        matrix.append(extract_floats(line))
    if asmatrix:  # return np.matrix object
        return np.asmatrix(matrix)
    else:  # return list of lists
        return matrix


def parse_qchem(fname, hooks, to_file=True, json_file='CCParser.json', 
                overwrite_vals=True, large_fn='matrices.npz', size_thresh=3):
    """Parse the QChem file for a list of hooks

    Parameters
    ---------
    fname : str
        Name of file to be parsed
    hooks : dict
        Dictionary with the following info:
        type of data ('matrix', 'vector', 'number'),
        hook_string ('Dipone moment')
        args [optional]: (line_shift, position in line)
    to_file : bool
        Whether it needs to be writen to file.
    json_file : str
        Name of json file where to save the parsed data.
    large_fn: str
        Name of npz file to save matrices/arrays of size above thresh
    size_thresh: int
        size of the matrix above which they are not saved in json
    """
    with open(fname, 'r') as ifile:
        lines = ifile.readlines()
    if not isinstance(hooks, dict):
        raise TypeError("`hooks` should be given as a dictionary.")
    parsed, large = {}, {}
    """
    if json_file is path, json_filepath = json_file
    if fname is filename, json_filepath = json_file
    if fname is path (path/fname.out), json_path is filename, saves in folder (path/jsfile.json)
    """
    json_filepath = json_file if os.path.split(json_file)[0] else os.path.join(os.path.split(fname)[0],json_file)
    large_filepath = os.path.join(os.path.split(json_filepath)[0], large_fn)
    if os.path.split(large_fn)[0]:
        large_fn = os.path.split(large_fn)[1]
        print("\"large_fn\" must be a filename. Everything before your filename is being ignored and\
                                large_fn is always placed in the same folder as \"json_file\"")
    for n, line in enumerate(lines):
        for key in hooks:
            args = None
            if len(hooks[key]) == 2:
                otype, hook = hooks[key]
            else:
                otype, hook, args = hooks[key]
            hook = re.compile(hook)
            match = hook.search(line)
            if match:
                # Get value(s)
                if otype == 'elconfs':
                    out = parse_elconf(n, lines, join=True)
                elif otype == 'geometry':
                    out = parse_molecule(n, lines, join=True)
                    if len(out) > size_thresh:
                        print("Matrix is too large, saving %s in file." % key)
                        large[key] = large[key]+[np.array(out, dtype="object")] if key in large.keys() else [np.array(out, dtype="object")]
                        out = large_fn
                elif otype == 'frag_geoms':
                    out = parse_molecule(n, lines, join=False)
                    if max([len(i) for i in out]) > size_thresh:
                        print("Matrix is too large, saving %s matrix in file." % key)
                        large[key] = large[key]+[np.array(out, dtype="object")] if key in large.keys() else [np.array(out, dtype="object")]
                        out = large_fn
                elif otype == 'simple matrix':
                    if args:
                        out = parse_simple_matrix(n, lines, **args)
                    else:
                        out = parse_simple_matrix(n, lines)
                    if len(out) > size_thresh:  # large matrix
                        print("Matrix is too large, saving %s matrix in file." % key)
                        large[key] = large[key]+[np.array(out)] if key in large.keys() else [np.array(out)]
                        out = large_fn
                elif otype == 'symmetric matrix':
                    out = parse_symmetric_matrix(n, lines)
                    if len(out) > size_thresh:  # large matrix
                        print("Matrix is too large, saving %s matrix in file." % key)
                        large[key] = large[key]+[np.array(out)] if key in large.keys() else [np.array(out)]
                        out = large_fn
                elif otype == 'vector':
                    out = parse_inline_vec(line)
                elif otype == 'number':
                    if args:
                        out = parse_number_qchem(n, lines, **args)
                    else:
                        out = parse_number_qchem(n, lines)
                else:
                    raise NotImplementedError('The requested type ({}) is not yet implemented.'.format(otype))
                # Save them in dictionary
                parsed[key] = parsed[key] + [[out,n]] if key in parsed.keys() else [[out,n]]
    if large:
        """
        if matrix_file is path, matrix_filepath = json_file
        if fname is filename, matrix_filepath = json_file
        if fname is path (path/fname.out), matrix_path is filename, saves in folder (path/matrix_file.npz)
        """
        np.savez(large_filepath, **large)
    if to_file:
        # Check if file exists and update dictionary
        if os.path.isfile(json_filepath):
            parsed = update_json_dict(json_filepath, parsed, overwrite_vals)
        # Save json TODO: replace with dump_js
        with open(json_filepath, 'w') as ofile:
            json.dump(parsed, ofile)
    return parsed


def update_json_dict(json_file, parsed, overwrite_vals):
    """Update existing dictionary.

    Parameters
    ----------
    json_file : str
        File name of the json file.
    parsed :  dict
        Dictionary with data freshly parsed.
    overwrite_vals : bool
        Whether to overwrite existing values or keep previous.

    Returns
    -------
    updated : dict
        Updated dictionary
    """
    # Read json TODO: replace with load_js
    with open(json_file, 'r') as ifile:
        updated = json.load(ifile)
    if overwrite_vals:
        updated.update(parsed)
    else:
        for key in parsed:
            if key not in updated:
                updated[key] = parsed[key]
    return updated


def read_matrix_from_json(json_file, keys):
    """Read parsed matrices saved in json_file.

    Parameters
    ----------
    json_file : str
        Name of json file to read.
    keys : list
        List of keys to look into
    """
    # Read json TODO: replace with load_js
    with open(json_file, 'r') as ifile:
        parsed = json.load(ifile)
    matrices = {}
    files = {}
    for key in keys:
        if key not in matrices:
            matrices[key] = []
        data = parsed[key]
        for info in data:
            if isinstance(info[0][0], str):
                fname = info[0][0]
                k = info[0][1]
                if fname not in files:
                    if fname.split('.')[-1] == 'npz':
                        files[fname] = np.load(fname)
                        matrices[key].append(files[fname][k])
                    else:
                        raise ValueError('Only npz files')
                else:
                    if fname.split('.')[-1] == 'npz':
                        matrices[key].append(files[fname][k])
                    else:
                        raise ValueError('Only npz files')
            else:
                matrices[key].append(np.array(info[0]))
    return matrices


def parse_number_qchem(n, lines, line_shift=0, position=-1):
    """Parse the the number from the line.

    Parameters
    ----------
    lines : str
        Lines from text that contains the number.
    n : int
        Number of the line where the hook starts.
    line_shift : int
        Number of lines that need to be shift to start reading.
    position : int
        From the splitted line, which position occupies
        the requiered number, default assumes the last part.

    Returns
    -------
    number : int or float

    """
    number = lines[n+line_shift].split()[position]
    try:
        number = int(number)
    except ValueError:
        number = float(number)
    return number
 

###########################################################################
# The hooks
# Each hook contains: (type, hook_string, extra_args)
# type: could be 'matrix', 'vector', 'number'
hooks = {"scf_energy" : ("number", "SCF   energy in the final basis set"),
         "overlap_matrix" : ("matrix", " Overlap Matrix"),
         "mp_energy": ("number", r"(RI-)?MP\([2-3]\) Summary", dict(line_shift=3, position=2)),
         "exc_energies": ("number", "Excitation energy:", dict(position=-2)),
         "osc_strength": ("number", "Osc. strength:"),
         "total_dipole" : ("number", "Total dipole"),
         "xyz": ('geometry', r"\$molecule"),
         "frag_xyz": ('frag_geoms', r"\$molecule"),
         "elconf": ('elconf', r"\$molecule"),
         "EFG_tensor_e": ('simple matrix', r"^  Raw EFG tensor \(electronic", dict(stop_signals=['Raw', ' EFG'])),
         "EFG_tensor_n": ('simple matrix', r"^  Raw EFG tensor \(nuclear", dict(stop_signals=['Raw', ' EFG'])),
         "EFG_tensor_t": ('simple matrix', r"^  Raw EFG tensor \(total", dict(stop_signals=['Principal'])),
        }


#if __name__ == "__main__":
#    fname = 'emb1.out'
#    parse_qchem(fname, hooks, to_file=True)
