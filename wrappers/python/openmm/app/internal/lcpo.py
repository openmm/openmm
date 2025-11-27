"""
lcpo.py: LCPO coefficient tables and setup.

This is part of the OpenMM molecular simulation toolkit.
See https://openmm.org/development.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Evan Pretti
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import math
import openmm as mm
import openmm.unit as u
import warnings

LCPO_PARAMETERS = {
    # For H atoms, virtual sites, etc.
    'none': (0.0, 0.0, 0.0, 0.0, 0.0),

    # Weiser, Shenkin and Still, J. Comput. Chem. 20, 217-230 (1999).
    'C_sp3_1': (1.7, 0.77887, -0.28063, -0.0012968, 0.00039328),
    'C_sp3_2': (1.7, 0.56482, -0.19608, -0.0010219, 0.0002658),
    'C_sp3_3': (1.7, 0.23348, -0.072627, -0.00020079, 0.00007967),
    'C_sp3_4': (1.7, 0.0, 0.0, 0.0, 0.0),
    'C_sp2_2': (1.7, 0.51245, -0.15966, -0.00019781, 0.00016392),
    'C_sp2_3': (1.7, 0.070344, -0.019015, -0.000022009, 0.000016875),
    'O_sp3_1': (1.6, 0.77914, -0.25262, -0.0016056, 0.00035071),
    'O_sp3_2': (1.6, 0.49392, -0.16038, -0.00015512, 0.00016453),
    'O_sp2_1': (1.6, 0.68563, -0.1868, -0.00135573, 0.00023743),
    'O_carboxylate': (1.6, 0.88857, -0.33421, -0.0018683, 0.00049372),
    'N_sp3_1': (1.65, 0.078602, -0.29198, -0.0006537, 0.00036247),
    'N_sp3_2': (1.65, 0.22599, -0.036648, -0.0012297, 0.000080038),
    'N_sp3_3': (1.65, 0.051481, -0.012603, -0.00032006, 0.000024774),
    'N_sp2_1': (1.65, 0.73511, -0.22116, -0.00089148, 0.0002523),
    'N_sp2_2': (1.65, 0.41102, -0.12254, -0.000075448, 0.00011804),
    'N_sp2_3': (1.65, 0.062577, -0.017874, -0.00008312, 0.000019849),
    'S_1': (1.9, 0.7722, -0.26393, 0.0010629, 0.0002179),
    'S_2': (1.9, 0.54581, -0.19477, -0.0012873, 0.00029247),
    'P_3': (1.9, 0.3865, -0.18249, -0.0036598, 0.0004264),
    'P_4': (1.9, 0.03873, -0.0089339, 0.0000083582, 0.0000030381),
    'Cl': (1.8, 0.98318, -0.40437, 0.00011249, 0.00049901),

    # AmberTools
    'F': (1.47, 0.68563, -0.1868, -0.00135573, 0.00023743),
    'Mg': (1.18, 0.49392, -0.16038, -0.00015512, 0.00016453),
}

def addLCPOForces(system, paramsList, usePeriodic, surfaceTension=0.005*u.kilocalorie_per_mole/u.angstrom**2, probeRadius=1.4*u.angstrom):
    """
    Adds forces to an OpenMM System implementing the LCPO method for estimating
    solvent-accessible surface area of a molecule.

    Parameters
    ----------
    system : System
        The OpenMM System to add forces to.
    paramsList : list
        A list containing LCPO parameters for each atom, specifically, a sphere
        radiuis that does not include a solvent probe radius, and coefficients
        P1 through P4 in the LCPO equations.
    surfaceTension : energy/area
        The energy per area to scale the surface area from the LCPO method by.
    probeRadius : distance
        The radius of the solvent probe to use.
    """

    force = mm.LCPOForce()
    force.setSurfaceTension(surfaceTension)
    for atomRadius, p1, p2, p3, p4 in paramsList:
        force.addParticle(atomRadius + probeRadius if atomRadius else 0, p1, p2, p3, p4)
    force.setUsesPeriodicBoundaryConditions(usePeriodic)
    system.addForce(force)

def getLCPOParamsAmber(prmtop, elements):
    """
    Generates LCPO parameters for each atom in an Amber prmtop file.

    Parameters
    ----------
    prmtop : PrmtopLoader
        The PrmtopLoader object containing data about the prmtop file.
    elements : list(Element)
        The elements of the atoms in the prmtop file.

    Returns
    -------
    list
        A list containing LCPO parameters for each atom, specifically, a sphere
        radiuis that does not include a solvent probe radius, and coefficients
        P1 through P4 in the LCPO equations.
    """

    numHeavyBondsList = [0] * prmtop.getNumAtoms()
    for atom1, atom2, _, _ in prmtop.getBondsNoH():
        numHeavyBondsList[atom1] += 1
        numHeavyBondsList[atom2] += 1

    numTotalBondsList = numHeavyBondsList.copy()
    for atom1, atom2, _, _ in prmtop.getBondsWithH():
        numTotalBondsList[atom1] += 1
        numTotalBondsList[atom2] += 1

    paramsList = []
    for atom, (element, numHeavyBonds, numTotalBonds) in enumerate(zip(elements, numHeavyBondsList, numTotalBondsList)):
        # Give atoms with no element 'none' parameters.
        atomicNumber = 0 if element is None else element.atomic_number
        atomType = prmtop.getAtomType(atom)

        params = LCPO_PARAMETERS['none']
        warn = False

        # Use Amber logic for selecting parameters, and Amber default fallback
        # behavior for cases where no suitable parameters can be assigned.
        if atomicNumber == 6:
            if numTotalBonds == 4:
                if numHeavyBonds == 1:
                    params = LCPO_PARAMETERS['C_sp3_1']
                elif numHeavyBonds == 2:
                    params = LCPO_PARAMETERS['C_sp3_2']
                elif numHeavyBonds == 3:
                    params = LCPO_PARAMETERS['C_sp3_3']
                elif numHeavyBonds == 4:
                    params = LCPO_PARAMETERS['C_sp3_4']
                else:
                    warn = True
                    params = LCPO_PARAMETERS['C_sp3_1']
            else:
                if numHeavyBonds == 2:
                    params = LCPO_PARAMETERS['C_sp2_2']
                elif numHeavyBonds == 3:
                    params = LCPO_PARAMETERS['C_sp2_3']
                else:
                    warn = True
                    params = LCPO_PARAMETERS['C_sp3_1']
        elif atomicNumber == 8:
            if atomType == 'O':
                params = LCPO_PARAMETERS['O_sp2_1']
            elif atomType == 'O2':
                params = LCPO_PARAMETERS['O_carboxylate']
            else:
                if numHeavyBonds == 1:
                    params = LCPO_PARAMETERS['O_sp3_1']
                elif numHeavyBonds == 2:
                    params = LCPO_PARAMETERS['O_sp3_2']
                else:
                    warn = True
                    params = LCPO_PARAMETERS['O_sp3_1']
        elif atomicNumber == 7:
            if atomType == 'N3':
                if numHeavyBonds == 1:
                    params = LCPO_PARAMETERS['N_sp3_1']
                elif numHeavyBonds == 2:
                    params = LCPO_PARAMETERS['N_sp3_2']
                elif numHeavyBonds == 3:
                    params = LCPO_PARAMETERS['N_sp3_3']
                else:
                    warn = True
                    params = LCPO_PARAMETERS['N_sp3_1']
            else:
                if numHeavyBonds == 1:
                    params = LCPO_PARAMETERS['N_sp2_1']
                elif numHeavyBonds == 2:
                    params = LCPO_PARAMETERS['N_sp2_2']
                elif numHeavyBonds == 3:
                    params = LCPO_PARAMETERS['N_sp2_3']
                else:
                    warn = True
                    params = LCPO_PARAMETERS['N_sp3_1']
        elif atomicNumber == 16:
            if atomType == 'SH':
                params = LCPO_PARAMETERS['S_1']
            else:
                params = LCPO_PARAMETERS['S_2']
        elif atomicNumber == 15:
            if numHeavyBonds == 3:
                params = LCPO_PARAMETERS['P_3']
            elif numHeavyBonds == 4:
                params = LCPO_PARAMETERS['P_4']
            else:
                warn = True
                params = LCPO_PARAMETERS['P_3']
        elif atomType.startswith('Z') or atomicNumber <= 1:
            # Use default 'none' parameters.
            pass
        elif atomType == 'MG':
            params = LCPO_PARAMETERS['Mg']
        elif atomType == 'F':
            params = LCPO_PARAMETERS['F']
        elif atomicNumber == 17:
            # Cl is the only element in the LCPO paper not implemented in Amber.
            params = LCPO_PARAMETERS['Cl']
        else:
            warn = True
            params = LCPO_PARAMETERS['C_sp2_2']

        if warn:
            warnings.warn(
                f'Using fallback LCPO parameters for atom with atomic number {atomicNumber}, '
                f'type {atomType}, {numTotalBonds} bonds, and {numHeavyBonds} bonds excluding H'
            )

        paramsList.append((params[0] * u.angstrom, params[1], params[2], params[3], params[4] / u.angstrom ** 2))

    return paramsList