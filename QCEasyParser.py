"""Some useful tools to comunicate with other QC software."""

import os
import re
import json
import numpy as np
from CCParser.QChem import parse_symmetric_matrix, parse_inline_vec 


def parse_qchem(fname, hooks, to_file=True, output='CCParser.json', overwrite_vals=True):
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
    output : str
        Name of file where to save the parsed data.
    """
    with open(fname, 'r') as ifile:
        lines = ifile.readlines()
    if not isinstance(hooks, dict):
        raise TypeError("`hooks` should be given as a dictionary.")
    parsed = {}
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
                if otype == 'matrix':
                    out = parse_symmetric_matrix_qchem(lines, n)
                if otype == 'vector':
                    out = parse_inline_vec(line)
                elif otype == 'number':
                    if args:
                        out = parse_number_qchem(lines, n, **args)
                    else:
                        out = parse_number_qchem(lines, n)
                else:
                    raise NotImplementedError('Only matrices and numbers can be parsed.')
                # Save them in dictionary
                if key in parsed:
                    parsed[key].append([out, n])
                else:
                    parsed[key] = [[out, n]]
    if to_file:
        # Check if file exists and update dictionary
        if os.path.isfile(output):
            parsed = update_json_dict(output, parsed, overwrite_vals)
        # Save json TODO: replace with dump_js
        with open(output, 'w') as ofile:
            json.dump(parsed, ofile)
    else:
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


def parse_number_qchem(lines, n, line_shift=0, position=-1):
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
        }


if __name__ == "__main__":
    fname = 'emb1.out'
    parse_qchem(fname, hooks, to_file=True)
