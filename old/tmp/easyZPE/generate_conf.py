#!/usr/bin/env python3

import sys, os, argparse, shutil
from os.path import join as pjoin

import numpy as np

import ase
from ase.io import read as ase_read
from ase.io.espresso import write_espresso_in
from ase import Atoms

from load_calculator import get_Espresso_input

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
    debug = args.debug

    infile = args.filename # QE input file with starting structure
    displ = args.displ # Ang
    indices = args.indices # which atoms must be displaced
    geom = ase_read(infile, format='espresso-in') # read initial geometry
    calculator_input = get_Espresso_input(infile) # read QE inputs 
    coord_label = ['x', 'y', 'z'] # shortcut for printing

    # quick check for silly mistakes
    if 'tprnfor' not in calculator_input['input_data']['control'].keys():
        print('!'*20, 'No tprnfor in QE inputs: forgot it?')
    elif not calculator_input['input_data']['control']['tprnfor']:
        print('!'*20, 'tprnfor set to False: how get forces?')
    if calculator_input['input_data']['control']['calculation'] != 'scf':
        print('!'*20, 'calculation is "%s" instead of "scf": you sure?' % calculator_input['input_data']['control']['calculation'])


    for ia, atm in enumerate(geom):
        # skip atom if not in selected indices
        if ia not in indices: continue
    
        for j in range(3): # x y z
    
            print('atm %i coord %s' % (ia, coord_label[j]))
            
            # --- positive displ
            c_geom = geom.copy()
            if debug: print('before %.3g' % (c_geom[ia].position[j]), end=' ')
            c_geom[ia].position[j] += displ
            if debug: print('after %.3g' % c_geom[ia].position[j])
            # --- --- write
            c_dir  = 'dpos_%i_%i' % (ia, j)
            os.makedirs(c_dir, exist_ok=True)
            with open(pjoin(args.outdir, c_dir, args.outfile), 'w') as outstr:
                write_espresso_in(outstr, c_geom, format='espresso-in', **calculator_input)
            for f in args.extra_cp:
                shutil.copy(f, pjoin(args.outdir, c_dir)) 
            
            # --- negative displ
            c_geom = geom.copy()
            if debug: print('before %.3g' % (c_geom[ia].position[j]), end=' ')
            c_geom[ia].position[j] -= displ
            if debug: print('after %.3g' % c_geom[ia].position[j])
            # --- --- write
            c_dir  = 'dneg_%i_%i' % (ia, j)
            os.makedirs(c_dir, exist_ok=True)
            with open(pjoin(args.outdir, c_dir, args.outfile), 'w') as outstr:
                write_espresso_in(outstr, c_geom, format='espresso-in', **calculator_input)
            for f in args.extra_cp:
                shutil.copy(f, pjoin(args.outdir, c_dir)) 

