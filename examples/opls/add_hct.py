"""
Add a table with hct parameters to a .dms file.

Copyright Schrodinger, LLC. All rights reserved.
"""
# original Author: ivan.tubert-brohman@schrodinger.com
# revised Author: bzhang@brooklyn.cuny.edu
"""
usage: $SCHRODINGER/run add_hct.py [options] <dms_file>

Adds a table with hct parameters to the .dms file specified in the command
line.

positional arguments:
  dms_file              .dms file to parameterize

optional arguments:
  -v, -version          Show the program's version and exit.
  -h, -help             Show this help message and exit.
  -param_file PARAM_FILE
                        path to parameter file (default: 'agbnp2.param')
"""

import sys, os, math, sqlite3
from schrodinger import structure
from schrodinger.infra import mm
from schrodinger.utils import cmdline
from schrodinger.structutils.analyze import evaluate_smarts

TABLE_NAME = 'hct';


# these are the parameters *in the order in which they appear* in the
# hct.param file. They must also match the names of the columns in the
# SQL table.

PARAM_NAMES = \
    'radius screened_radius'.split()

CHARGE = 'r_m_charge1'
ID_PROP = 'i_dms_id'

class Agbnp2Error(Exception):
    """A class for exceptions raised by this module."""
    pass

def read_structure(conn):
    """
    Read a structure from the particle and bond tables of a database
    connection. Returns a schrodinger.structure.Structure object.
    The atoms get properties added with their database id and their
    sigma and epsilon parameters.
    """
    st = structure.create_new_structure()
    c = conn.cursor()
    atoms = {}
    for id, anum, x, y, z, charge in c.execute(
            "select particle.id, anum, x, y, z,charge "
            "from particle, nonbonded_param "
            "where particle.nbtype==nonbonded_param.id"):
        el = mm.mmat_get_element_by_atomic_number(anum)
        atom = atoms[id] = st.addAtom(el, x, y, z)
        atom.property[CHARGE] = charge
        atom.property[ID_PROP] = id             
        
    for p0, p1, order in c.execute('select p0, p1, "order" from bond'):
        st.addBond(atoms[p0], atoms[p1], int(order))

    return st


def read_hct_file(fname):
    """
    read a hct.param file. returns a list of atom types, where each atom 
    type is a tuple(smarts,params_dict).
    """
    #header and sample line:
    #type smart_pattern radius screened_radius
    #1    [#1][#6]      1.3      0.85
    atom_hct = []
    for line in open(fname):
        if line.startswith('#') or line.startswith('dielectric'): continue
        cols = line.split()
        smarts = cols[1]
        params_list = [float(s) for s in cols[2:]]
        if len(params_list) != len(PARAM_NAMES):
            print len(params_list)
            print len(PARAM_NAMES)
            raise Agbnp2Error("invalid number of parameters in line:\n %s"
                % line)
        params_dict = dict(zip(PARAM_NAMES, params_list))
        atom_hct.append((smarts, params_dict))
    return atom_hct
        

def add_HCT_table(conn,atom_hct,table_name):
    """
    create or replace a table with the HCT parameters for the atoms in the structure described by
    the particle and bond tables of the database connection.
    HCT needs three parameters:
    -charge: formal_charge for each atom
    -radius: radius-offset for each atom, generally we set offset==0.009
    -screen: screen number for each atom, initially we set screen==0.8*(radius-0.009)
    -the format for HCT table will be 
    id charge radius-offset screen
    """
    st = read_structure(conn)
    c  = conn.cursor()
    c.execute("drop table if exists %s" % table_name)
    c.execute("CREATE TABLE %s (id INTEGER NOT NULL REFERENCES particle, "
              "charge FLOAT, radius FLOAT, screened_radius FLOAT)"
              % table_name)
    atoms = {}
    for smarts, params in atom_hct:
        for match in evaluate_smarts(st,smarts):
            iatom = match[0]
            if iatom not in atoms:
                atom = atoms[iatom]=params.copy()
                atom['id'] = st.atom[iatom].property[ID_PROP]
                atom['charge'] = st.atom[iatom].property[CHARGE]

    
    for iatom in xrange(1,len(st.atom)+1):
        if iatom not in atoms:
            raise Agbnp2Error("parameters not found for atom %d" % iatom)
        cols = ['id'] + ['charge']+ PARAM_NAMES
        cmd = "INSERT into %s (%s) values(%s)" % (table_name,
               ','.join(cols),
               ','.join([':%s' % name for name in cols]))
        c.execute(cmd, atoms[iatom])
    conn.commit()


def add_hct_to_dms_file(dms_file, param_file):
    """
    Create or replace a table with the hct parameters for the atoms in the
    structure contained in dms_file, using the parameters specified in
    param_file.
    """
    if not os.path.isfile(dms_file):
        raise Agbnp2Error("file does not exist: %s" % dms_file)
    conn = sqlite3.connect(dms_file)
    atom_hct = read_hct_file(param_file)
    add_HCT_table(conn, atom_hct, TABLE_NAME)
    conn.close()

def parse_args():
    parser = cmdline.create_argument_parser(
        usage="$SCHRODINGER/run add_hct.py [options] <dms_file>",
        description="Adds a table with hct parameters to the .dms file "
            "specified in the command line.")
    parser.add_argument("dms_file", help=".dms file to parameterize")
    parser.add_argument("-param_file",
        help="path to parameter file (default: 'hct.param')",
        default="hct.param")
    return parser.parse_args()

def main():
    args = parse_args()
    try:
        add_hct_to_dms_file(args.dms_file, args.param_file)
    except (Agbnp2Error, IOError) as e:
        print >> sys.stderr, "Error: %s" % e
        sys.exit(1)

if __name__ == '__main__':
    main()

