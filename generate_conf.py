#!/usr/bin/env python3

import sys, os, argparse, shutil, re
from os.path import join as pjoin

import numpy as np

from ase.io import read as ase_read
from ase.io.espresso import write_espresso_in
from ase.units import Bohr

from load_calculator import get_Espresso_input


def _recursive_print(d, level, indent=30, fill='.'):
    """Recursively pretty-print a nested dictionary."""
    for key, val in d.items():
        if isinstance(val, dict):
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key), n=indent, c=fill))
            _recursive_print(val, level+1, indent=indent, fill=fill)
        else:
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key), n=indent, c=fill), val)


def _del_magnetisations(calculator_input):
    """Return a copy of calculator_input with starting_magnetization keys removed.

    ASE will re-set magnetisations from the Atoms object's magnetic moments,
    so these must be stripped before passing to write_espresso_in to avoid
    conflicts.
    """
    calc = calculator_input.copy()
    calc['input_data'] = {k: v.copy() if isinstance(v, dict) else v
                          for k, v in calculator_input['input_data'].items()}
    to_del = [k for k in calc['input_data']['system']
              if re.match(r'starting_magnetization\([0-9]+\)', k)]
    for k in to_del:
        del calc['input_data']['system'][k]
    return calc


def generate_conf(indices, displ, geom, calculator_input, extra_cp,
                  outdir, outfile, debug=False):
    """Write QE input files with +-delta displacements for each selected atom and direction.

    For N atoms in `indices`, creates 6N directories named dpos_<i>_<j> and
    dneg_<i>_<j> (i = atom index, j = 0/1/2 for x/y/z) under `outdir`.
    Each directory contains `outfile` (a ready-to-run pw.x input) and any
    files listed in `extra_cp`.

    Magnetic moments are taken from geom.get_initial_magnetic_moments() if set
    (i.e. when the source pw.in had starting_magnetization entries that ase_read
    parsed into the Atoms object). ASE's write_espresso_in re-derives
    starting_magnetization from those moments, so no manual handling is needed.

    Parameters
    ----------
    indices : list of int
        0-based indices of atoms to displace.
    displ : float
        Displacement magnitude [Angstrom].
    geom : ase.Atoms
        Equilibrium geometry. Any initial magnetic moments set on this object
        are propagated to every displaced configuration.
    calculator_input : dict
        QE input parameters as returned by load_calculator.get_Espresso_input().
    extra_cp : list of str
        Paths to extra files to copy into each displacement folder
        (e.g. HPC submission scripts).
    outdir : str
        Root directory under which displacement folders are created.
    outfile : str
        Filename for the QE input inside each displacement folder.
    debug : bool
        Print atom positions before and after displacement.
    """
    coord_label = ['x', 'y', 'z']

    print('Displacement: %.6f Ang  (%.6f a0)' % (displ, displ / Bohr))

    if debug:
        _recursive_print(calculator_input, 0)

    # Sanity checks on the input before touching anything
    ctrl = calculator_input['input_data']['control']
    if 'tprnfor' not in ctrl:
        print('WARNING: tprnfor not found in &CONTROL -- forces will not be printed!')
    elif not ctrl['tprnfor']:
        print('WARNING: tprnfor = .false. -- forces will not be printed!')
    if ctrl.get('calculation', 'scf') != 'scf':
        print('WARNING: calculation = "%s" instead of "scf"' % ctrl['calculation'])

    # Strip existing magnetisation keys from the input dict; ASE will re-derive
    # starting_magnetization from geom.get_initial_magnetic_moments() automatically.
    calc_input = _del_magnetisations(calculator_input)

    if 'hubbardU' in calc_input:
        raise NotImplementedError(
            'HUBBARD card found in input but stock ASE write_espresso_in '
            'does not support it. Use a patched espresso.py with hubbardU= support.')

    print('Creating displaced configurations:')
    for ia, atom in enumerate(geom):
        if ia not in indices:
            continue

        print('  atom %3i (%2s) at %s Ang' % (
            ia, atom.symbol, '%.4f %.4f %.4f' % tuple(atom.position)))

        for j in range(3):
            for sign, prefix in [(+1, 'dpos'), (-1, 'dneg')]:
                c_geom = geom.copy()
                if debug:
                    print('    %s %s: %.6f -> %.6f Ang' % (
                        prefix, coord_label[j],
                        c_geom[ia].position[j],
                        c_geom[ia].position[j] + sign * displ))
                c_geom[ia].position[j] += sign * displ

                dname = '%s_%i_%i' % (prefix, ia, j)
                dest = pjoin(outdir, dname)
                os.makedirs(dest, exist_ok=True)

                with open(pjoin(dest, outfile), 'w') as f:
                    write_espresso_in(f, c_geom, format='espresso-in', **calc_input)

                for fpath in extra_cp:
                    shutil.copy(fpath, dest)

        print('    ' + '-' * 48)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate QE pw.x inputs with finite displacements for ZPE evaluation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('indices', type=int, nargs='+',
                        help='0-based indices of atoms to displace')
    parser.add_argument('-f', '--filename', type=str, default='pw.in',
                        help='equilibrium QE input file')
    parser.add_argument('-d', '--displ', type=float, default=0.01,
                        help='displacement magnitude [Ang]')
    parser.add_argument('--outdir', type=str, default='./',
                        help='root directory for displacement folders')
    parser.add_argument('--outfile', type=str, default='pw.in',
                        help='filename for the QE input in each displacement folder')
    parser.add_argument('--extra_cp', type=str, nargs='*', default=[],
                        help='extra files to copy into each folder (e.g. job submission script)')
    parser.add_argument('--debug', action='store_true',
                        help='print atom positions before and after displacement')

    args = parser.parse_args()

    geom = ase_read(args.filename, format='espresso-in')
    calculator_input = get_Espresso_input(args.filename)

    generate_conf(
        indices=args.indices,
        displ=args.displ,
        geom=geom,
        calculator_input=calculator_input,
        extra_cp=args.extra_cp,
        outdir=args.outdir,
        outfile=args.outfile,
        debug=args.debug,
    )
