#!/usr/bin/env python3

"""Parser for Quantum Espresso input. Read a set of Namespaces and Card. Works but not all namespaces and cards are implemented"""

import sys, logging
from typing import List, Dict

def isint(string):
    try:
        int(string)
    except ValueError:
        return False
    return True

def isbool(string):
    if string.lower() in ['true', '.true.', 'false', '.false.']: return True
    return False

def str2bool(string):
    if string.lower() in ['true', '.true.']: return True
    elif string.lower() in ['false', '.false.']: return False
    else: raise ValueError('cannot convert %s to bool' % string)

def isfloat(string):
    try:
        float(string)
    except ValueError:
        return False
    return True

def get_Espresso_input(fname, debug=False):

    c_log = logging.getLogger('Espresso_inptu')
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        c_log.setLevel(logging.DEBUG)
        c_log.debug('PRINTING DEBUG INFO')
        c_log.warning('test')

    def listval2dict(lines, char='='):
        ddict = {}
        for l in lines:
            ls = l.split('=')
            key = ls[0].strip()
            val = ls[1].strip()
            if isint(val): val = int(val)
            elif isfloat(val): val = float(val)
            elif isbool(val): val= str2bool(val)
            ddict[key] = val
        return ddict

    # Functions to parse each block
    # each returns the object expected from Espresso calculator: name, dictionary if part of input_data or the relevant thing
    def get_control(lines):
        return 'control', listval2dict(lines)

    def get_system(lines):
        return 'system', listval2dict(lines)

    def get_electrons(lines):
        return 'electrons', listval2dict(lines)

    def get_ions(lines):
        return 'ions', listval2dict(lines)

    def get_cell(lines):
        return 'cell', listval2dict(lines)

    def get_atomic_species(lines):
        # what to do here? Pseudopotential? Mass?
        pseudos = {}
        masses = {}
        for l in lines[1:]:
            ls = l.split()
            symbol, mass, pseudo_name = ls[0].strip(), float(ls[1].strip()), ls[2].strip()
            pseudos[symbol] = pseudo_name
            masses[symbol] = mass
        return  'pseudopotentials', pseudos

    def get_kpoints(lines):
        ktype = lines[0].split()[1].strip()

        if ktype == 'automatic':
            return 'kpts', [int(kk) for kk in lines[1].split()[:3]]
        elif ktype == 'gamma':
            return 'kpts', None
        elif ktype == 'crystal':
            c_log.debug('assume list: ka kb kc weight. !! NOT WELL TESTED !!')
            kpts = []
            for l in lines[1:]:
                kpts.append([float(ii) for ii in l])
            return 'kpts', kpts
        elif ktype in ['tpiba', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c']:
            raise NotImplementedError('Kpoint %s not implemented' % ktype)
        else:
            raise RuntimeError('Something went wrong')

    def get_hubbard(lines):
        # !!! based on our own implementation of Espresso read/write function!
        hubbardU = {'params': []}
        hubbardU['option'] = lines[0].split()[1]
        for l in lines[1:]:
            ls = l.split()
            if ls[0] != 'U': raise NotImplemented('Only parsing U')
            elem, manifold = ls[1].split('-')
            Uval = float(ls[2])
            hubbardU['params'].append([elem, manifold, Uval])
        return 'hubbardU', hubbardU

    # Dictionary with field starting string and function to parse
    # Mandatory blocks
    name_must = {'&CONTROL': get_control,
                 '&SYSTEM': get_system,
                 '&ELECTRONS': get_electrons}
    # Optional blocks (with ending character)
    name_opt = {'&IONS': get_ions,
                '&CELL': get_cell}
    end_char = '/' # these blocks must end this this

    # These blocks (no ending char)
    card_must = {'ATOMIC_SPECIES': get_atomic_species,
                 'K_POINTS': get_kpoints}
    card_opt = {'HUBBARD': get_hubbard}

    # These blocks are read by ASE
    block_skip = ['CELL_PARAMETERS',
                 'ATOMIC_POSITIONS']

    lines = []
    with open(fname) as instr:
        for l in instr.readlines():
            l = l.strip().replace('"','').replace("'", "")
            c_log.debug(l)
            if not len(l): continue # skip empty lines
            if l[0] in ['!', '#']: # skip comments/empty lines
                c_log.debug('skip %s' % l)
            lines.append(l)

    # sanity check
    sanity_check = {key: False for key in [*name_must.keys(), *card_must.keys()]}
    for l in lines:
        lsplit = l.split()
        head = lsplit[0].strip()
        if head in sanity_check.keys(): sanity_check[head] = True
    c_log.debug(sanity_check)

    if not all(sanity_check.values()):
        raise RuntimeError('Mandatory field %s not found' % ' '.join([k for k,v in sanity_check.items() if not v]))

    flg_name = False
    flg_card = False
    flg_parse = False

    field_lines: List[str] # annotations!
    calculator_input = {'input_data': {}}

    for l in lines:
        lsplit = l.split()
        head = lsplit[0].strip()

        if head in name_must or head in name_opt:
            c_log.debug('starting namespace %s' % l)
            flg_name = True
            flg_parse = True
            field_key = head
            field_lines = []

        elif head in card_must or head in card_opt:
            if flg_card:
                c_log.debug('end field')
                flg_card = False
                c_log.debug('-- process %s' % field_key)
                c_log.debug(field_lines)
                k, v = card_must[field_key](field_lines) if field_key in card_must else card_opt[field_key](field_lines)
                c_log.debug('%s %s' % (str(k),str(v)))
                calculator_input[k] = v
                c_log.debug('---')
            c_log.debug('starting card %s' % l)
            flg_card = True
            flg_parse = True
            field_key = head
            field_lines = []
            field_lines.append(l)

        elif head in block_skip:
            if flg_card:
                c_log.debug('end field')
                flg_card = False
                flg_parse = False
                c_log.debug('-- process %s' % field_key)
                c_log.debug(field_lines)
                k, v = card_must[field_key](field_lines) if field_key in card_must else card_opt[field_key](field_lines)
                calculator_input[k] = v
                c_log.debug('---')
            c_log.debug('skip field %s' % l)
            field_key = head
            field_lines = []

        elif head == end_char and (flg_name):
            c_log.debug('End field char')
            c_log.debug('-- process %s' % field_key)
            c_log.debug(field_lines)
            flg_name = False
            flg_parse = False
            k, v = name_must[field_key](field_lines) if field_key in name_must else name_opt[field_key](field_lines)
            calculator_input['input_data'][k] = v
            c_log.debug('---')
        elif flg_parse:
            c_log.debug('-- %s: %s' % (field_key, l))
            field_lines.append(l)
        else:
            c_log.debug('-- SKIP: %s' % (l))
    c_log.debug('==========')
    return calculator_input


def recursive_print(ddict, level, indent=20, fill='.'):
    """Recursive print a dictionary: if value of is a dictionary, call this function again with increased nested level for printing"""
    for key, val in ddict.items():
        if type(val) == dict:
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key),n=indent, c=fill))
            recursive_print(val, level+1)
        else:
            print(' '*indent*level + '{s:{c}<{n}}'.format(s=str(key),n=indent,c=fill), val)

if __name__ == '__main__':
    fname = sys.argv[1]
    calculator_input = get_Espresso_input(fname, debug=False)
    recursive_print(calculator_input, 0)
