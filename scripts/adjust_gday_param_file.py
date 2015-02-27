#!/usr/bin/env python

"""
Change various params in a G'DAY input file.

All that needs to be supplied is a dictionary with the variable to be changed
and the corresponding paramter value (oh and of course the param fname!).
"""

import os
import re
import tempfile
import shutil

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.06.2011)"
__email__   = "mdekauwe@gmail.com"

def adjust_param_file(fname, replacements):
    """ adjust model parameters in the file and save over the original.

    Parameters:
    ----------
    fname : string
        parameter filename to be changed.
    replacements : dictionary
        dictionary of replacement values.

    """
    fin = open(fname, 'r')
    param_str = fin.read()
    fin.close()
    new_str = replace_keys(param_str, replacements)
    fd, path = tempfile.mkstemp()
    os.write(fd, new_str)
    os.close(fd)
    shutil.copy(path, fname)
    os.remove(path)

def replace_keys(text, replacements_dict):
    """ Function expects to find GDAY input file formatted key = value.

    Parameters:
    ----------
    text : string
        input file data.
    replacements_dict : dictionary
        dictionary of replacement values.

    Returns:
    --------
    new_text : string
        input file with replacement values

    """
    lines = text.splitlines()
    for i, row in enumerate(lines):
        # skip blank lines
        if not row.strip():
            continue
        # skip .cfg section dividers
        elif not row.startswith("["):
            key, sep, val = row.split()
            lines[i] = " ".join((key, sep, replacements_dict.get(key, val)))
        elif row.startswith("[print]"):
            break


    return '\n'.join(lines) + '\n'


if __name__ == "__main__":

    cfg_fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/dk_varyco2_varyndep_grassequilib_then_forest_dukegrass_youngforest.cfg"
    replace_dict = { "albedo": "0.123"}
    adjust_param_file(cfg_fname, replace_dict)
    print cfg_fname
