"""
compiled.pyx: Utility functions that are compiled with Cython for speed

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2018 Stanford University and the Authors.
Authors: Peter Eastman
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
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

from heapq import heappush, heappop

cdef extern from "math.h":
    double round(double x)
    double sqrt(double x)

cdef class periodicDistance:
    """This is a callable object that computes the distance between two points, taking
    periodic boundary conditions into account.  This is used heavily in Modeller.addSolvent().
    """
    cdef double vectors[3][3]
    cdef double invBoxSize[3]
    
    def __init__(self, boxVectors):
        for i in range(3):
            for j in range(3):
                self.vectors[i][j] = boxVectors[i][j]
            self.invBoxSize[i] = 1.0/boxVectors[i][i]

    def __call__(self, pos1, pos2):
        cdef double dx, dy, dz, scale1, scale2, scale3
        dx = pos1[0]-pos2[0]
        dy = pos1[1]-pos2[1]
        dz = pos1[2]-pos2[2]
        scale3 = round(dz*self.invBoxSize[2])
        dx -= scale3*self.vectors[2][0]
        dy -= scale3*self.vectors[2][1]
        dz -= scale3*self.vectors[2][2]
        scale2 = round(dy*self.invBoxSize[1])
        dx -= scale2*self.vectors[1][0]
        dy -= scale2*self.vectors[1][1]
        scale1 = round(dx*self.invBoxSize[0])
        dx -= scale1*self.vectors[0][0]
        return sqrt(dx*dx + dy*dy + dz*dz)


def matchResidueToTemplate(res, template, bondedToAtom, bint ignoreExternalBonds=False):
    """Determine whether a residue matches a template and return a list of corresponding atoms.
    This is used heavily in ForceField.

    Parameters
    ----------
    res : Residue
        The residue to check
    template : _TemplateData
        The template to compare it to
    bondedToAtom : list
        Enumerates which other atoms each atom is bonded to
    ignoreExternalBonds : bool
        If true, ignore external bonds when matching templates

    Returns
    -------
    list
        a list specifying which atom of the template each atom of the residue
        corresponds to, or None if it does not match the template
    """
    cdef int numAtoms, i, j
    atoms = list(res.atoms())
    numAtoms = len(atoms)
    if numAtoms != len(template.atoms):
        return None

    # Translate from global to local atom indices, and record the bonds for each atom.

    renumberAtoms = {}
    for i in range(numAtoms):
        renumberAtoms[atoms[i].index] = i
    bondedTo = []
    externalBonds = []
    for atom in atoms:
        bonds = [renumberAtoms[x] for x in bondedToAtom[atom.index] if x in renumberAtoms]
        bondedTo.append(bonds)
        externalBonds.append(0 if ignoreExternalBonds else len([x for x in bondedToAtom[atom.index] if x not in renumberAtoms]))

    # For each unique combination of element and number of bonds, make sure the residue and
    # template have the same number of atoms.

    residueTypeCount = {}
    for i, atom in enumerate(atoms):
        key = (atom.element, len(bondedTo[i]), externalBonds[i])
        if key not in residueTypeCount:
            residueTypeCount[key] = 1
        residueTypeCount[key] += 1
    templateTypeCount = {}
    for i, atom in enumerate(template.atoms):
        key = (atom.element, len(atom.bondedTo), 0 if ignoreExternalBonds else atom.externalBonds)
        if key not in templateTypeCount:
            templateTypeCount[key] = 1
        templateTypeCount[key] += 1
    if residueTypeCount != templateTypeCount:
        return None

    # Identify template atoms that could potentially be matches for each atom.

    candidates = [[] for i in range(numAtoms)]
    cdef bint exactNameMatch
    for i in range(numAtoms):
        exactNameMatch = (atoms[i].element is None and any(atom.element is None and atom.name == atoms[i].name for atom in template.atoms))
        for j, atom in enumerate(template.atoms):
            if (atom.element is not None and atom.element != atoms[i].element) or (exactNameMatch and atom.name != atoms[i].name):
                continue
            if len(atom.bondedTo) != len(bondedTo[i]):
                continue
            if not ignoreExternalBonds and atom.externalBonds != externalBonds[i]:
                continue
            candidates[i].append(j)

    # Find an optimal ordering for matching atoms.  This means 1) start with the one that has the fewest options,
    # and 2) follow with ones that are bonded to an already matched atom.

    searchOrder = []
    atomsToOrder = set(range(numAtoms))
    efficientAtomSet = set()
    efficientAtomHeap = []
    while len(atomsToOrder) > 0:
        if len(efficientAtomSet) == 0:
            fewestNeighbors = numAtoms+1
            for i in atomsToOrder:
                if len(candidates[i]) < fewestNeighbors:
                    nextAtom = i
                    fewestNeighbors = len(candidates[i])
        else:
            nextAtom = heappop(efficientAtomHeap)[1]
            efficientAtomSet.remove(nextAtom)
        searchOrder.append(nextAtom)
        atomsToOrder.remove(nextAtom)
        for i in bondedTo[nextAtom]:
            if i in atomsToOrder:
                if i not in efficientAtomSet:
                    efficientAtomSet.add(i)
                    heappush(efficientAtomHeap, (len(candidates[i]), i))
    inverseSearchOrder = [0]*numAtoms
    for i in range(numAtoms):
        inverseSearchOrder[searchOrder[i]] = i
    bondedTo = [[inverseSearchOrder[bondedTo[i][j]] for j in range(len(bondedTo[i]))] for i in searchOrder]
    candidates = [candidates[i] for i in searchOrder]

    # Recursively match atoms.

    matches = numAtoms*[0]
    hasMatch = numAtoms*[False]
    if _findAtomMatches(template, bondedTo, matches, hasMatch, candidates, 0):
        return [matches[inverseSearchOrder[i]] for i in range(numAtoms)]
    return None


def _getAtomMatchCandidates(template, bondedTo, matches, candidates, position):
    """Get a list of template atoms that are potential matches for the next atom."""
    for bonded in bondedTo[position]:
        if bonded < position:
            # This atom is bonded to another one for which we already have a match, so only consider
            # template atoms that *that* one is bonded to.
            return template.atoms[matches[bonded]].bondedTo
    return candidates[position]


def _findAtomMatches(template, bondedTo, matches, hasMatch, candidates, int position):
    """This is called recursively from inside matchResidueToTemplate() to identify matching atoms."""
    if position == len(matches):
        return True
    cdef int i
    for i in _getAtomMatchCandidates(template, bondedTo, matches, candidates, position):
        atom = template.atoms[i]
        if not hasMatch[i] and i in candidates[position]:
            # See if the bonds for this identification are consistent

            allBondsMatch = all((bonded > position or matches[bonded] in atom.bondedTo for bonded in bondedTo[position]))
            if allBondsMatch:
                # This is a possible match, so try matching the rest of the residue.

                matches[position] = i
                hasMatch[i] = True
                if _findAtomMatches(template, bondedTo, matches, hasMatch, candidates, position+1):
                    return True
                hasMatch[i] = False
    return False
