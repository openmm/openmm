"""
This package contains code for reading CHARMM structure files for setting up a
simulation with CHARMM; specifically PSF, PAR, RTF, and STR files

    - PAR : Parameter file (PRM) -- this contains all of the force field
            parameters (e.g., bond lengths and strengths) for all of the atom
            types

    - RTF : Residue Topology File -- this file contains the residue
            connectivity tables as well as a definition of all of the atom
            types. Also contains an internal coordinate representation of the
            residues

    - PSF : Protein Structure File -- this is the main file type in CHARMM
            simulations that defines all of the residues in a system as well as
            the atom types and connectivity between the atoms

    - STR : Stream file -- Source of additional information and CHARMM commands
            that can contain RTF and PAR information. Allows users to define
            additional parameters without 'contaminating' the original force
            field parameter files
"""

__all__ = ['psf', 'parameters']
__authors__ = 'Jason Swails'
__contributors__ = ''
__license__ = 'MIT'
__copyright__ = 'Copyright 2014, Jason Swails'
__date__ = 'Apr. 18, 2014'


__private__ = ['topologyobjects', '_charmmfile']
