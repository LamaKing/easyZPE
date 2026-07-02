#!/usr/bin/env python3

import sys, os, argparse, shutil
import re
from os.path import join as pjoin

import numpy as np

import ase
from ase.io import read as ase_read
from ase.io.espresso import write_espresso_in
from ase import Atoms
from ase.units import Bohr
from ase.io.espresso import write_espresso_in, grep_valence, kspacing_to_grid

from load_calculator import get_Espresso_input

def recursive_print(ddict, level, indent=30, fill='.'):
    """Recursive print a dictionary: if value of is a dictionary, call this function again with increased nested level for printing"""
    for key, val in ddict.items():
        if type(val) == dict:
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key),n=indent, c=fill))
            recursive_print(val, level+1)
        else:
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key),n=indent,c=fill), val)

def get_magnetisation(calculator_input):
    elements = [a for a in calculator_input['pseudopotentials'].keys()]
    magnetisations = {a: 0  for a in elements}
    for key, val in calculator_input['input_data']['system'].items():
        match = re.match('starting_magnetization\(([0-9])+\)', key)
        if match:
            ii = int(match.group(1))-1 # QE starts from 1, python starts from 0
            magnetisations[elements[ii]] = val
    return magnetisations

def del_magnetisations(original_input):
    calculator_input = original_input.copy()

    to_del = []
    for key, val in calculator_input['input_data']['system'].items():
        match = re.match('starting_magnetization\(([0-9])+\)', key)
        if match:
            to_del.append(key)
    for k in to_del:
        del calculator_input['input_data']['system'][k]

    return calculator_input

def generate_conf(indices, displ, geom, calculator_input, extra_cp, outdir, outfile, magnetisations=None, debug=False):

    print('displ %f [Ang]' %  displ, end=' ')
    print('%f [a0]' % (displ/Bohr)) # convert to units of QE XML output

    if debug:
        recursive_print(calculator_input, 0)

    # set magnetisation right! MAYBE MAKES MORE SENSE TO DO IT OUTSIDE THIS FUNCTION?
    #magnetisations = get_magnetisation(calculator_input)   # not necessary, ASE does it automatically
    #print(magnetisations)
    calculator_input = del_magnetisations(calculator_input) # delete magnetisations from input. ASE will overrite IF atoms has magnetic moments set!
    if magnetisations:
        pseudo_dir = os.environ['ESPRESSO_PSEUDO']
        geom.set_initial_magnetic_moments([magnetisations[atom.symbol]*grep_valence(pjoin(pseudo_dir,
                                                                                      calculator_input['pseudopotentials'][atom.symbol])
                                                                                )
                                           for atom in geom])

    coord_label = ['x', 'y', 'z'] # shortcut for printing

    # quick check for silly mistakes
    if 'tprnfor' not in calculator_input['input_data']['control'].keys():
        print('!'*20, 'No tprnfor in QE inputs: forgot it?')
    elif not calculator_input['input_data']['control']['tprnfor']:
        print('!'*20, 'tprnfor set to False: how get forces?')
    if calculator_input['input_data']['control']['calculation'] != 'scf':
        print('!'*20, 'calculation is "%s" instead of "scf": you sure?' % calculator_input['input_data']['control']['calculation'])

    print('creating geometries with displacements:')
    for ia, atm in enumerate(geom):
        # skip atom if not in selected indices
        if ia not in indices: continue

        print('atom %3i (%2s) original coord %s Ang' % (ia, atm.symbol, '%.2f %.2f %.2f' % tuple(atm.position)))
        for j in range(3): # x y z

            print('displace coord %s' % (coord_label[j]))

            # --- positive displ
            c_geom = geom.copy()
            if debug: print('before %.3g' % (c_geom[ia].position[j]), end=' ')
            c_geom[ia].position[j] += displ
            if debug: print('after %.3g' % c_geom[ia].position[j])
            # --- --- write
            c_dir  = 'dpos_%i_%i' % (ia, j)
            os.makedirs(pjoin(outdir, c_dir), exist_ok=True)
            with open(pjoin(outdir, c_dir, outfile), 'w') as outstr:
                write_espresso_in(outstr, c_geom, format='espresso-in', **calculator_input)
            for f in extra_cp:
                shutil.copy(f, pjoin(outdir, c_dir))

            # --- negative displ
            c_geom = geom.copy()
            if debug: print('before %.3g' % (c_geom[ia].position[j]), end=' ')
            c_geom[ia].position[j] -= displ
            if debug: print('after %.3g' % c_geom[ia].position[j])
            # --- --- write
            c_dir  = 'dneg_%i_%i' % (ia, j)
            os.makedirs(pjoin(outdir, c_dir), exist_ok=True)
            with open(pjoin(outdir, c_dir, outfile), 'w') as outstr:
                write_espresso_in(outstr, c_geom, format='espresso-in', **calculator_input)
            for f in extra_cp:
                shutil.copy(f, pjoin(outdir, c_dir))

        print('-'*48)

if __name__ == "__main__":

    # --- argument parser
    parser = argparse.ArgumentParser(description="Generate QuantumEspresso inputs with finite displacement for Zero Point Energy evaluation.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # --- --- positional arguments
    parser.add_argument('indices',
                        type=int, nargs='+',
                        help='indices of atoms displaced')
    # --- --- optional args
    parser.add_argument('--extra_cp',
                        type=str, nargs='*', default=[],
                        help='extra files to copy in directory, e.g. submission scrip for HPC')
    parser.add_argument('--outdir',
                        type=str, default='./',
                        help='where to create folders with displacements')
    parser.add_argument('--outfile',
                        type=str, default='pw.in',
                        help='filename for QE input with displacement')
    parser.add_argument('-d', '--displ',
                        type=float, default=0.01,
                        help='displacement used in finite difference [Ang]')
    parser.add_argument('-f', '--filename',
                        type=str, default='pw.in',
                        help='QE input file')
    parser.add_argument('--debug',
                        action='store_true',
                        help='print debug informations')
    # --- --- init
    args = parser.parse_args(sys.argv[1:])

    infile = args.filename # QE input file with starting structure
    displ = args.displ # Ang
    indices = args.indices # which atoms must be displaced
    geom = ase_read(infile, format='espresso-in') # read initial geometry
    calculator_input = get_Espresso_input(infile) # read QE inputs

    extra_cp = args.extra_cp
    outdir, outfile = args.outdir, args.outfile

    generate_conf(indices, displ, geom, calculator_input, extra_cp, outdir, outfile, debug=args.debug)
