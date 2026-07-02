#!/usr/bin/env python3

import sys, argparse, os
from os.path import join as pjoin

import numpy as np

import qeschema

from ase.io import read as ase_read
from ase.units import Bohr


# Physical constants
_CONV2THZ = 968.2886664991329 / (2 * np.pi)  # sqrt(Eh/(Da*a0^2)) -> THz
_THZ2CMINV = 33.356683                         # THz -> cm^-1
_HBAR = 6.582119568038699e-4                   # eV / THz  (= hbar in eV*ps)
_KB = 8.617333262e-5                           # eV / K


def qestruct2ase(xml_atomic_structure):
    """Build an ASE Atoms object from the atomic_structure block of a QE XML output.

    Parameters
    ----------
    xml_atomic_structure : dict
        The 'atomic_structure' subtree returned by qeschema's to_dict().

    Returns
    -------
    ase.Atoms
        Geometry in Angstrom with periodic boundary conditions set.
    """
    from ase import Atom, Atoms

    cell = Bohr * np.array(list(xml_atomic_structure['cell'].values()))
    atoms = [
        Atom(a['@name'], position=Bohr * np.array(a['$']), index=a['@index'])
        for a in xml_atomic_structure['atomic_positions']['atom']
    ]
    return Atoms(atoms, cell=cell, pbc=[True, True, True])


def _read_qe_forces(xml_path):
    """Read forces from a QE XML output file.

    Parameters
    ----------
    xml_path : str
        Path to pwscf.xml (or equivalent).

    Returns
    -------
    qe_root : dict
        Full qes:espresso dict (for geometry checks).
    forces : np.ndarray, shape (nat, 3)
        Forces in Hartree/Bohr.
    """
    qe = qeschema.PwDocument(xml_path).to_dict(validation='lax')['qes:espresso']
    dim = list(reversed(qe['output']['forces']['@dims']))  # QE stores transposed vs numpy
    forces = np.array(qe['output']['forces']['$']).reshape(dim)
    return qe, forces


def compute_freq(indices, displ, geom, indir, qeout, outdir, debug, temperature=298.15):
    """Build the dynamical matrix from QE finite-displacement data and diagonalise it.

    For each atom in `indices` and each Cartesian direction, reads forces from
    the +delta and -delta QE runs produced by generate_conf.py, assembles the
    3N x 3N force-constant matrix via central finite differences, mass-weights
    it to obtain the dynamical matrix at Gamma, and diagonalises with eigh.

    Skips modes with omega <= 0 (translational / soft modes due to DFT noise)
    when computing ZPE and thermal quantities.

    Parameters
    ----------
    indices : list of int
        0-based indices of the displaced atoms (must match generate_conf run).
    displ : float
        Displacement magnitude used in generate_conf [Angstrom].
    geom : ase.Atoms
        Equilibrium geometry.
    indir : str
        Directory containing the dpos_*/dneg_* displacement folders.
    qeout : str
        Name of the QE XML output file inside each displacement folder.
    outdir : str
        Directory where output files are written.
    debug : bool
        Print intermediate matrices.
    temperature : float, optional
        Temperature for thermal corrections [K]. Default 298.15 K.

    Returns
    -------
    eigval : np.ndarray, shape (3N,)
        Eigenvalues of the dynamical matrix [Eh/(Da*a0^2)], sorted ascending.
    eigvec : np.ndarray, shape (3N, 3N)
        Eigenvectors as columns.
    freqs : np.ndarray, shape (3N,)
        Frequencies in THz (negative for imaginary modes).
    """
    coord_label = ['x', 'y', 'z']

    print('displacement %.6f [Ang]  %.6f [a0]' % (displ, displ / Bohr))
    displ_au = displ / Bohr  # convert to Bohr for consistency with QE XML

    N = len(indices)
    # Walk indices in the same order they will be used to fill force matrix rows,
    # so mass_vec[ib] always corresponds to the atom at row ib.
    masses = np.array([geom.get_masses()[ii] for ii in indices])  # Da, in indices order

    print('Displaced atoms (%i):' % N,
          ' '.join(['%i (%s)' % (i, geom[i].symbol) for i in indices]))
    print('Masses:',
          ' '.join(['%s: %.4g Da' % (geom[i].symbol, m) for i, m in zip(indices, masses)]))

    # Force matrices: rows = (displaced atom, coord), cols = (force-on atom, coord)
    F_pos = np.zeros((3*N, 3*N))
    F_neg = np.zeros((3*N, 3*N))
    print('Force matrix shape:', F_pos.shape)

    for ib, ii in enumerate(indices):
        for j in range(3):
            row = ib * 3 + j
            print('  Load atom %i (%s) disp %s (row %i)' % (
                ii, geom[ii].symbol, coord_label[j], row))

            for sign, prefix, F_matr in [(+1, 'dpos', F_pos), (-1, 'dneg', F_neg)]:
                xml_path = pjoin(indir, '%s_%i_%i' % (prefix, ii, j), qeout)
                qe_root, forces = _read_qe_forces(xml_path)

                # Verify the displacement in the XML matches what was requested
                disp_geom = qestruct2ase(qe_root['output']['atomic_structure'])
                actual_disp = (disp_geom.positions - geom.positions)[ii, j] / Bohr
                if not np.isclose(actual_disp, sign * displ_au):
                    raise RuntimeError(
                        'Displacement mismatch for atom %i coord %s %s: '
                        'expected %.6f a0, got %.6f a0' % (
                            ii, coord_label[j], prefix, sign * displ_au, actual_disp))

                for ic, n in enumerate(indices):
                    for g in range(3):
                        col = ic * 3 + g
                        if debug:
                            print('    %s force on atom %i %s (mat[%i,%i]): %9.3g Eh/a0'
                                  % (prefix, n, coord_label[g], row, col, forces[n, g]))
                        F_matr[row, col] = forces[n, g]

            if debug:
                print('    ' + '-' * 40)

    print()

    if debug:
        print('F(+delta) [Eh/a0]')
        print(F_pos)
        print('F(-delta) [Eh/a0]')
        print(F_neg)

    np.savetxt(pjoin(outdir, 'force_posdx.dat'), F_pos, fmt='%25.15g', header='Eh/a0')
    np.savetxt(pjoin(outdir, 'force_negdx.dat'), F_neg, fmt='%25.15g', header='Eh/a0')

    # Central finite difference: Phi_ij = -dF_i/du_j
    force_const = -(F_pos - F_neg) / (2 * displ_au)  # Eh/a0^2

    if debug:
        print('Force constants [Eh/a0^2]')
        np.set_printoptions(precision=3)
        print(force_const)

    np.savetxt(pjoin(outdir, 'force_const.dat'), force_const, fmt='%25.15g', header='Eh/a0^2')

    # Mass-weight: D_ij = Phi_ij / sqrt(m_i * m_j)
    # Build mass vector aligned with the (atom, coord) indexing of force_const
    mass_vec = np.repeat(masses, 3)  # [m0,m0,m0, m1,m1,m1, ...]
    dyn_matr = force_const / np.sqrt(np.outer(mass_vec, mass_vec))  # Eh/(a0^2 Da)

    if debug:
        print('Dynamical matrix [Eh/(a0^2 Da)]')
        print(dyn_matr)
        print()

    # Symmetrise before eigh to suppress residual asymmetry from DFT noise
    dyn_matr_sym = (dyn_matr + dyn_matr.T) / 2
    np.savetxt(pjoin(outdir, 'dyna_matr.dat'), dyn_matr_sym, fmt='%25.15g', header='Eh/(a0^2 Da)')

    print('Diagonalising...')
    eigval, eigvec = np.linalg.eigh(dyn_matr_sym)
    print()

    if debug:
        print('Eigenvalues:', eigval)
        print('Eigenvectors:')
        print(eigvec)

    # Flag genuinely negative eigenvalues (not just numerical noise near zero)
    thold = 1e-10
    for val in eigval:
        if val < -thold:
            print('WARNING: negative eigenvalue %.6g Eh/(Da*a0^2)' % val)

    # Encode imaginary modes as negative frequencies for output clarity.
    freqs = np.sign(eigval) * np.sqrt(np.abs(eigval)) * _CONV2THZ

    # Print summary table
    print('         ' + ' '*15, ' '.join(['     l%-4s' % i for i in range(len(eigval))]))
    print('eigenvals %-15s' % '[Eh/(Da*a0^2)]', ' '.join(['%10.4f' % e for e in eigval]))
    print('freq      %-15s' % '[THz]',           ' '.join(['%10.4f' % f for f in freqs]))
    print('freq      %-15s' % '[cm^-1]',         ' '.join(['%10.4f' % f for f in freqs * _THZ2CMINV]))

    # Harmonic vibrational thermodynamics
    KBT = _KB * temperature
    ZPEs  = np.zeros(len(freqs))
    Evibs = np.zeros(len(freqs))
    TSs   = np.zeros(len(freqs))
    Avibs = np.zeros(len(freqs))

    for io, freq in enumerate(freqs):
        omega = freq * 2 * np.pi  # THz -> rad/THz (angular frequency)
        if omega <= 0:
            print('WARNING: mode %i has omega <= 0 (%.4g THz), skipping thermodynamics' % (io, freq),
                  file=sys.stderr)
            continue
        x = _HBAR * omega / KBT   # hbar*omega / kT, dimensionless
        ZPEs[io]  = _HBAR * omega / 2
        Evibs[io] = _HBAR * omega * (0.5 + 1.0 / (np.exp(x) - 1))
        Avibs[io] = _HBAR * omega / 2 + KBT * np.log(1 - np.exp(-x))
        TSs[io]   = Evibs[io] - Avibs[io]

    Tlabel = '%.0fK' % temperature
    print('ZPE         %-15s' % '[eV]', ' '.join(['%10.4f' % e for e in ZPEs]),
          'sum: %10.4g' % np.sum(ZPEs))
    print('Evib        %-15s' % '[eV]', ' '.join(['%10.4f' % e for e in Evibs]),
          'sum: %10.4g' % np.sum(Evibs))
    print('TS(%s)%s'   % (Tlabel, ' '*(9-len(Tlabel))) + '%-15s' % '[eV]',
          ' '.join(['%10.4f' % e for e in TSs]),
          'sum: %10.4g' % np.sum(TSs))
    print('Avib(%s)%s' % (Tlabel, ' '*(8-len(Tlabel))) + '%-15s' % '[eV]',
          ' '.join(['%10.4f' % e for e in Avibs]),
          'sum: %10.4g' % np.sum(Avibs))

    # Save per-mode data
    with open(pjoin(outdir, 'eigenvalues.dat'), 'w') as outf:
        header = '#%-24s ' % 'eig[Eh/(Da*a0^2)]' + ' '.join(
            '%-25s' % s for s in ['freq[THz]', 'freq[cm^-1]', 'ZPE[eV]',
                                   'TS(%s)[eV]' % Tlabel, 'Avib(%s)[eV]' % Tlabel])
        print(header, file=outf)
        for io, (eig, freq) in enumerate(zip(eigval, freqs)):
            print(' '.join('%-25.15g' % v for v in [
                eig, freq, freq * _THZ2CMINV, ZPEs[io], TSs[io], Avibs[io]
            ]), file=outf)

    # Print eigenvectors grouped by atom
    print('Eigenvectors (columns):')
    print(' '*7 + ' '.join(' v%-3s' % i for i in range(eigvec.shape[1])))
    for ir in range(eigvec.shape[0]):
        if ir % 3 == 0:
            print(' '*6, '-' * 6 * eigvec.shape[1])
        atom_sym = geom[indices[ir // 3]].symbol
        print('%3s_%s:' % (atom_sym, coord_label[ir % 3]),
              ' '.join('%5.2f' % v for v in eigvec[ir]))

    np.savetxt(pjoin(outdir, 'eigenvectors.dat'), eigvec, fmt='%10.5f',
               header='Eigenvectors as columns. Eigenvalues: %s' % str(eigval))

    return eigval, eigvec, freqs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Collect QE finite-displacement output and diagonalise the "
                    "dynamical matrix at Gamma to obtain ZPE and vibrational free energy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('indices', type=int, nargs='+',
                        help='0-based indices of the displaced atoms')
    parser.add_argument('-f', '--filename', type=str, default='pw.in',
                        help='equilibrium QE input file (for masses and geometry check)')
    parser.add_argument('-d', '--displ', type=float, default=0.01,
                        help='displacement magnitude used in generate_conf [Ang]')
    parser.add_argument('--indir', type=str, default='./',
                        help='directory containing the dpos_*/dneg_* folders')
    parser.add_argument('--outdir', type=str, default='./',
                        help='directory for output files')
    parser.add_argument('--qeout', type=str, default='pwscf.xml',
                        help='QE XML output filename (same for all displacement runs)')
    parser.add_argument('-T', '--temperature', type=float, default=298.15,
                        help='temperature for thermal corrections [K]')
    parser.add_argument('--debug', action='store_true',
                        help='print intermediate matrices')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    geom = ase_read(args.filename, format='espresso-in')

    compute_freq(
        indices=args.indices,
        displ=args.displ,
        geom=geom,
        indir=args.indir,
        qeout=args.qeout,
        outdir=args.outdir,
        debug=args.debug,
        temperature=args.temperature,
    )
