"""
tinkerfiles.py: Used for loading TINKER key/prm and xyz files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Joao Morado
Contributors: Peter Eastman

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

__author__ = "Joao Morado"

import datetime
import io
import math
import os
import re
import shlex
import xml.etree.ElementTree as etree
import numpy as np
from functools import wraps
from typing import Any, Dict, List, Set, Tuple, Union

from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import nanometers
from openmm.vec3 import Vec3

from . import element as elem
from . import forcefield as ff
from . import topology as top



class TinkerFiles:
    """TinkerFiles parses Tinker files (.xyz, .prm, .key), constructs a Topology, and (optionally) an OpenMM System from it."""


    RECOGNIZED_FORCES: Dict[str, Any] = {
        "bond": 1,
        "angle": 1,
        "anglep": 1,
        "strbnd": 1,
        "ureybrad": 1,
        "opbend": 1,
        "torsion": 1,
        "pitors": 1,
        "strtors": 1,
        "angtors": 1,
        "vdw": 1,
        "vdwpr": 1,
        "polarize": 1,
        "tortors": None,
        "multipole": None,
    }

    RECOGNIZED_SCALARS: Dict[str, str] = {
        "forcefield": "-2.55",
        "bond-cubic": "-2.55",
        "bond-quartic": "3.793125",
        "angle-cubic": "-0.014",
        "angle-quartic": "0.000056",
        "angle-pentic": "-0.0000007",
        "angle-sextic": "0.000000022",
        "opbendtype": "ALLINGER",
        "opbend-cubic": "-0.014",
        "opbend-quartic": "0.000056",
        "opbend-pentic": "-0.0000007",
        "opbend-sextic": "0.000000022",
        "torsionunit": "0.5",
        "vdwtype": "BUFFERED-14-7",
        "radiusrule": "CUBIC-MEAN",
        "radiustype": "R-MIN",
        "radiussize": "DIAMETER",
        "epsilonrule": "HHG",
        "dielectric": "1.0",
        "polarization": "MUTUAL",
        "vdw-13-scale": "0.0",
        "vdw-14-scale": "1.0",
        "vdw-15-scale": "1.0",
        "mpole-12-scale": "0.0",
        "mpole-13-scale": "0.0",
        "mpole-14-scale": "0.4",
        "mpole-15-scale": "0.8",
        "polar-12-scale": "0.0",
        "polar-13-scale": "0.0",
        "polar-14-scale": "1.0",
        "polar-15-scale": "1.0",
        "polar-14-intra": "0.5",
        "direct-11-scale": "0.0",
        "direct-12-scale": "1.0",
        "direct-13-scale": "1.0",
        "direct-14-scale": "1.0",
        "mutual-11-scale": "1.0",
        "mutual-12-scale": "1.0",
        "mutual-13-scale": "1.0",
        "mutual-14-scale": "1.0",
    }

    @staticmethod
    def _initialize_standard_biopolymer_and_biotypes():
        """Helper function to initialize standard biopolymer and biotype dictionaries"""
        # Amino Acid Codes
        _AMINO_ACID_CODES = {
            'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I',
            'SER': 'S', 'THR': 'T', 'CYS': 'C', 'CYX': 'C', 'CYD': 'c',
            'PRO': 'P', 'PHE': 'F', 'TYR': 'Y', 'TYD': 'y', 'TRP': 'W',
            'HIS': 'H', 'HID': 'U', 'HIE': 'Z', 'ASP': 'D', 'ASH': 'd',
            'ASN': 'N', 'GLU': 'E', 'GLH': 'e', 'GLN': 'Q', 'MET': 'M',
            'LYS': 'K', 'LYD': 'k', 'ARG': 'R', 'ORN': 'O', 'AIB': 'B',
            'PCA': 'J', 'H2N': 't', 'FOR': 'f', 'ACE': 'a', 'COH': 'x',
            'NH2': 'n', 'NME': 'm', 'UNK': 'X'
        }
        _AMINO_ACID_ORDER = [
            'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'CYS', 'CYX', 'CYD',
            'PRO', 'PHE', 'TYR', 'TYD', 'TRP', 'HIS', 'HID', 'HIE', 'ASP', 'ASH',
            'ASN', 'GLU', 'GLH', 'GLN', 'MET', 'LYS', 'LYD', 'ARG', 'ORN', 'AIB',
            'PCA', 'H2N', 'FOR', 'ACE', 'COH', 'NH2', 'NME', 'UNK'
        ]

        # Nucleotide Codes
        _NUCLEOTIDE_CODES = {
            '  A': 'A', '  G': 'G', '  C': 'C', '  U': 'U', ' DA': 'D',
            ' DG': 'B', ' DC': 'I', ' DT': 'T', ' MP': '1', ' DP': '2',
            ' TP': '3', 'UNK': 'X'
        }
        _NUCLEOTIDE_ORDER = [ 
            '  A', '  G', '  C', '  U', ' DA', ' DG', ' DC', ' DT',
            ' MP', ' DP', ' TP', 'UNK'
        ]

        _AMINO_ACID_BIOTYPES_MID = {}
        _AMINO_ACID_BIOTYPES_NTERM = {}
        _AMINO_ACID_BIOTYPES_CTERM = {}
        _NUCLEOTIDE_BIOTYPES = {}

        # Biopolymer types for mid-chain peptide backbone atoms
        nt   = [1,   7,  15,  27,  41,  55,  65,  77,  87,  96, 105, 116, 131, 147, 162, 185, 202, 218, 234, 244, 256, 268, 280, 294, 308, 321, 337, 353, 370, 384, 391,   0,   0,   0,   0,   0,   0,   1]
        cat  = [2,   8,  16,  28,  42,  56,  66,  78,  88,  97, 106, 117, 132, 148, 163, 186, 203, 219, 235, 245, 257, 269, 281, 295, 309, 322, 338, 354, 371, 385, 392,   0,   0,   0,   0,   0,   0,   2]
        ct   = [3,   9,  17,  29,  43,  57,  67,  79,  89,  98, 107, 118, 133, 149, 164, 187, 204, 220, 236, 246, 258, 270, 282, 296, 310, 323, 339, 355, 372, 386, 393,   0,   0,   0,   0,   0,   0,   3]
        hnt  = [4,  10,  18,  30,  44,  58,  68,  80,  90,  99,   0, 119, 134, 150, 165, 188, 205, 221, 237, 247, 259, 271, 283, 297, 311, 324, 340, 356, 373, 387, 394,   0,   0,   0,   0,   0,   0,   4]
        ot   = [5,  11,  19,  31,  45,  59,  69,  81,  91, 100, 108, 120, 135, 151, 166, 189, 206, 222, 238, 248, 260, 272, 284, 298, 312, 325, 341, 357, 374, 388, 395,   0,   0,   0,   0,   0,   0,   5]
        hat  = [6,  12,  20,  32,  46,  60,  70,  82,  92, 101, 109, 121, 136, 152, 167, 190, 207, 223, 239, 249, 261, 273, 285, 299, 313, 326, 342, 358, 375,   0, 396,   0,   0,   0,   0,   0,   0,   6]
        cbt  = [0,  13,  21,  33,  47,  61,  71,  83,  93, 102, 110, 122, 137, 153, 168, 191, 208, 224, 240, 250, 262, 274, 286, 300, 314, 327, 343, 359, 376, 389, 397,   0,   0,   0,   0,   0,   0,   0]
        for i, res in enumerate(_AMINO_ACID_ORDER):
            _AMINO_ACID_BIOTYPES_MID[res] = {
                'N': nt[i], 'CA': cat[i], 'C': ct[i], 'H': hnt[i], 'O': ot[i], 'HA': hat[i], 'CB': cbt[i]
            }

        # Biopolymer types for N-terminal peptide backbone atoms
        nn   = [403, 409, 415, 421, 427, 433, 439, 445, 451, 457, 463, 471, 477, 483, 489, 495, 501, 507, 513, 519, 525, 531, 537, 543, 549, 555, 561, 567, 573, 579, 391, 762,   0,   0,   0,   0,   0, 403]
        can  = [404, 410, 416, 422, 428, 434, 440, 446, 452, 458, 464, 472, 478, 484, 490, 496, 502, 508, 514, 520, 526, 532, 538, 544, 550, 556, 562, 568, 574, 580, 392,   0,   0, 767,   0,   0,   0, 404]
        cn   = [405, 411, 417, 423, 429, 435, 441, 447, 453, 459, 465, 473, 479, 485, 491, 497, 503, 509, 515, 521, 527, 533, 539, 545, 551, 557, 563, 569, 575, 581, 393,   0, 764, 769,   0,   0,   0, 405]
        hnn  = [406, 412, 418, 424, 430, 436, 442, 448, 454, 460, 466, 474, 480, 486, 492, 498, 504, 510, 516, 522, 528, 534, 540, 546, 552, 558, 564, 570, 576, 582, 394, 763,   0,   0,   0,   0,   0, 406] # HN1, HN2, HN3 usually
        on   = [407, 413, 419, 425, 431, 437, 443, 449, 455, 461, 467, 475, 481, 487, 493, 499, 505, 511, 517, 523, 529, 535, 541, 547, 553, 559, 565, 571, 577, 583, 395,   0, 766, 770,   0,   0,   0, 407]
        han  = [408, 414, 420, 426, 432, 438, 444, 450, 456, 462, 468, 476, 482, 488, 494, 500, 506, 512, 518, 524, 530, 536, 542, 548, 554, 560, 566, 572, 578,   0, 396,   0, 765, 768,   0,   0,   0, 408]
        for i, res in enumerate(_AMINO_ACID_ORDER):
            _AMINO_ACID_BIOTYPES_NTERM[res] = {
                'N': nn[i], 'CA': can[i], 'C': cn[i], 'H': hnn[i], 'O': on[i], 'HA': han[i]
                # Note: CB is usually the same as mid-chain, sidechain atoms need separate lookup
            }

        # Biopolymer types for C-terminal peptide backbone atoms
        nc   = [584, 590, 596, 602, 608, 614, 620, 626, 632, 638, 644, 649, 655, 661, 667, 673, 679, 685, 691, 697, 703, 709, 715, 721, 727, 733, 739, 745, 751, 757,   0,   0,   0,   0, 773, 775, 777, 584]
        cac  = [585, 591, 597, 603, 609, 615, 621, 627, 633, 639, 645, 650, 656, 662, 668, 674, 680, 686, 692, 698, 704, 710, 716, 722, 728, 734, 740, 746, 752, 758,   0,   0,   0,   0,   0,   0, 779, 585]
        cc   = [586, 592, 598, 604, 610, 616, 622, 628, 634, 640, 646, 651, 657, 663, 669, 675, 681, 687, 693, 699, 705, 711, 717, 723, 729, 735, 741, 747, 753, 759,   0,   0,   0,   0, 771,   0,   0, 586]
        hnc  = [587, 593, 599, 605, 611, 617, 623, 629, 635, 641,   0, 652, 658, 664, 670, 676, 682, 688, 694, 700, 706, 712, 718, 724, 730, 736, 742, 748, 754, 760,   0,   0,   0,   0, 774, 776, 778, 587]
        oc   = [588, 594, 600, 606, 612, 618, 624, 630, 636, 642, 647, 653, 659, 665, 671, 677, 683, 689, 695, 701, 707, 713, 719, 725, 731, 737, 743, 749, 755, 761,   0,   0,   0,   0, 772,   0,   0, 588] # OXT, O'
        hac  = [589, 595, 601, 607, 613, 619, 625, 631, 637, 643, 648, 654, 660, 666, 672, 678, 684, 690, 696, 702, 708, 714, 720, 726, 732, 738, 744, 750, 756,   0,   0,   0,   0,   0,   0,   0, 780, 589]
        for i, res in enumerate(_AMINO_ACID_ORDER):
            _AMINO_ACID_BIOTYPES_CTERM[res] = {
                'N': nc[i], 'CA': cac[i], 'C': cc[i], 'H': hnc[i], 'OXT': oc[i], 'HA': hac[i]
                # Note: CB is usually the same as mid-chain, sidechain atoms need separate lookup
            }

        # Biopolymer types for nucleotide phosphate and sugar atoms
        o5t   = [1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203, 0, 0, 0, 0]
        c5t   = [1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204, 0, 0, 0, 0]
        h51t  = [1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205, 0, 0, 0, 0] # H5'
        h52t  = [1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206, 0, 0, 0, 0] # H5''
        c4t   = [1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207, 0, 0, 0, 0]
        h4t   = [1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208, 0, 0, 0, 0] # H4'
        o4t   = [1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209, 0, 0, 0, 0]
        c1t   = [1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210, 0, 0, 0, 0]
        h1t   = [1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211, 0, 0, 0, 0] # H1'
        c3t   = [1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212, 0, 0, 0, 0]
        h3t   = [1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213, 0, 0, 0, 0] # H3'
        c2t   = [1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214, 0, 0, 0, 0]
        h21t  = [1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215, 0, 0, 0, 0] # H2'
        o2t   = [1014, 1044, 1075, 1103,    0,    0,    0,    0, 0, 0, 0, 0] # O2' (RNA only)
        h22t  = [1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216, 0, 0, 0, 0] # H2'' (often same type as H2')
        o3t   = [1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217, 0, 0, 0, 0]
        pt    = [1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242, 0, 0, 0, 0]
        opt   = [1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243, 0, 0, 0, 0] # OP1, OP2
        h5tt  = [1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245, 0, 0, 0, 0] # HO5' (5' terminal H)
        h3tt  = [1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250, 0, 0, 0, 0] # HO3' (3' terminal H)
        for i, res in enumerate(_NUCLEOTIDE_ORDER):
            _NUCLEOTIDE_BIOTYPES[res] = {
                 "O5'": o5t[i], "C5'": c5t[i], "H5'": h51t[i], "H5''": h52t[i],
                 "C4'": c4t[i], "H4'": h4t[i], "O4'": o4t[i], "C1'": c1t[i],
                 "H1'": h1t[i], "C3'": c3t[i], "H3'": h3t[i], "C2'": c2t[i],
                 "H2'": h21t[i], "O2'": o2t[i], "H2''": h22t[i], "O3'": o3t[i],
                 "P": pt[i], "OP1": opt[i], "OP2": opt[i], # OP1/OP2 often share a type
                 "HO5'": h5tt[i], "HO3'": h3tt[i]
                 # Base atoms need separate lookup based on A/G/C/T/U
             }

        return (
            _AMINO_ACID_CODES,
            _AMINO_ACID_ORDER,
            _NUCLEOTIDE_CODES,
            _NUCLEOTIDE_ORDER,
            _AMINO_ACID_BIOTYPES_MID,
            _AMINO_ACID_BIOTYPES_NTERM,
            _AMINO_ACID_BIOTYPES_CTERM,
            _NUCLEOTIDE_BIOTYPES
        )

    # Call the helper function to populate the dictionaries when the class is defined
    (_AMINO_ACID_CODES, 
     _AMINO_ACID_ORDER,
     _NUCLEOTIDE_CODES, 
     _NUCLEOTIDE_ORDER,
     _AMINO_ACID_BIOTYPES_MID, 
     _AMINO_ACID_BIOTYPES_NTERM, 
     _AMINO_ACID_BIOTYPES_CTERM, 
     _NUCLEOTIDE_BIOTYPES) = _initialize_standard_biopolymer_and_biotypes()

    def __init__(
        self,
        xyz: str,
        key: str,
        seq: str = None,
        periodicBoxVectors: Tuple[Vec3, Vec3, Vec3] = None,
        unitCellDimensions: Vec3 = None,
        writeXmlFiles: bool = False,
    ):
        """
        Load exactly one .xyz file and one or more .key or .prm files.

        Parameters
        ----------
        xyz : str
            The path to the xyz file to load.
        key : str or list of str
            The path(s) to the key/prm file(s) to load.
        seq : str, optional
            The path to the .seq file to load.
        periodicBoxVectors : tuple of Vec3
            The vectors defining the periodic box.
            If provided, this overwrites the box information from the xyz file.
        unitCellDimensions : Vec3, optional
            The dimensions of the crystallographic unit cell.
            For non-rectangular unit cells, specify periodicBoxVectors instead.
            If provided, this overwrites the box information from the xyz file.
        writeXmlFiles : bool, optional, default=False
            If True, the residue and force field XML files are written to disk.
        """
        # ----------------------- INTERNAL VARIABLES -----------------------
        # Populate parser functions for the recognized forces
        self.RECOGNIZED_FORCES["tortors"] = TinkerFiles.__addTorTor
        self.RECOGNIZED_FORCES["multipole"] = TinkerFiles.__addMultipole

        # Store the input parameters
        self._writeXmlFiles = writeXmlFiles
        self._XmlFilesList = None

        # Topology
        self.topology = None

        # Position and box information
        self.positions = None
        self.boxVectors = None
        self._numpyPositions = None
        self._numpyBoxVectors = None

        # Internal variable to store the data from the xyz file
        self._xyzDict = None

        # Internal variables to store the data from the key file(s)
        self._atomTypes, self._bioTypes, self._forces, self._scalars = [], [], [], []

        # Internal variable to store the combined atom data

        # ----------------------- LOAD FILES -----------------------
        # Load the .seq file
        if seq is not None:
            self._seqDict = self.readSeqFile(seq)
        else:
            self._seqDict = None


        # Load the .xyz file
        self._xyzDict, self.boxVectors, self.positions = self._loadSeqFile(xyz)
        self.positions = self.positions * nanometers
    
        # Load the .key or .prm file(s)
        key = key if isinstance(key, list) else [key]
        for keyFile in key:
            atomTypes, bioTypes, forces, scalars = self._loadKeyFile(keyFile)
            self._atomTypes.append(atomTypes)
            self._bioTypes.append(bioTypes)
            self._forces.append(forces)
            self._scalars.append(scalars)

        # Load the .xyz file
        self._xyzDict, self.boxVectors, self.positions = self._loadXyzFile(xyz)
        self.positions = self.positions * nanometers
        
        # Combine the data from the .xyz and .key files
        self._atomData = TinkerFiles._combineXyzAndKeyData(
            self._xyzDict, self._atomTypes, self._bioTypes
        )

        # ----------------------- CREATE TOPOLOGY -----------------------
        # Create the topology
        self.topology = self._createTopology(self._atomData)

        # Set the periodic box vectors as specified in the xyz file
        if self.boxVectors is not None:
            self.topology.setPeriodicBoxVectors(self.boxVectors)

        # If provided, set the periodic box vectors or unit cell dimensions in the topology
        # This overwrites the box information from the xyz file
        if periodicBoxVectors is not None and unitCellDimensions is not None:
            raise ValueError(
                "Specify either periodicBoxVectors or unitCellDimensions, but not both"
            )
        if periodicBoxVectors is not None:
            self.topology.setPeriodicBoxVectors(periodicBoxVectors)
        elif unitCellDimensions is not None:
            self.topology.setUnitCellDimensions(unitCellDimensions)

        # ----------------------- CREATE XML FILES -----------------------
        self.forceFields = []
        self.implictForceFields = []
        for keyFile, atomTypes, bioTypes, forces, scalars in zip(
            key, self._atomTypes, self._bioTypes, self._forces, self._scalars
        ):
            xmlFile, implicitXmlFile = self._createXmlFile(
                keyFile, atomTypes, bioTypes, forces, scalars, self._atomData
            )

            # Reset the file pointers
            xmlFile.seek(0)
            implicitXmlFile.seek(0)

            # Store the XML files
            self.forceFields.append(xmlFile)
            self.implictForceFields.append(implicitXmlFile)

    def createSystem(
        self,
        nonbondedMethod = ff.PME,
        nonbondedCutoff = 1.0 * nanometers,
        constraints = None,
        rigidWater: bool = False,
        removeCMMotion: bool = True,
        hydrogenMass = None,
        polarization: str = "mutual",
        mutualInducedTargetEpsilon: float = 0.00001,
        implicitSolvent: bool = False,
        *args,
        **kwargs,
    ) -> Any:
        """
        Create an OpenMM System from the parsed Tinker files.

        Parameters
        ----------
        nonbondedMethod : object=PME
            The method to use for nonbonded interactions.
            Allowed values are NoCutoff, and PME.
        nonbondedCutoff : distance=1.0*nanometers
            The cutoff distance to use for nonbonded interactions.
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : bool, optional, default=True
            If true, water molecules will be fully rigid regardless of the value passed for the constraints argument.
            Note that AMOEBA waters are flexible.
        removeCMMotion : bool, optional, default=True
            If True, center of mass motion will be removed.
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.
            Any mass added to a hydrogen is subtracted from the heavy atom to keep their total mass the same.
        polarization : str, optional, default="mutual"
            The method to use for calculating induced dipoles.
            Allowed values are "mutual", "direct", or "extrapolated".
        mutualInducedTargetEpsilon : float, optional, default=0.00001
            The target epsilon for mutual induced dipoles.
            Only used if polarization="mutual".
        implicitSolvent : bool, optional, default=False
            If True, solvent will be modeled implicitly.

        Returns
        -------
        openmm.System
            The created OpenMM System.
        """
        xmlFilesList = self.forceFields
        if implicitSolvent:
            xmlFilesList += self.implictForceFields

        # Reset the file pointers
        for f in xmlFilesList:
            f.seek(0)

        forcefield = ff.ForceField(*xmlFilesList)

        system = forcefield.createSystem(
            self.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            removeCMMotion=removeCMMotion,
            hydrogenMass=hydrogenMass,
            polarization=polarization,
            mutualInducedTargetEpsilon=mutualInducedTargetEpsilon,
            *args,
            **kwargs,
        )

        return system

    # ------------------------------------------------------------------------------------------ #
    #                                     POSITIONS AND BOX                                      #
    # ------------------------------------------------------------------------------------------ #
    def getPositions(self, asNumpy: bool = False) -> Union[List[Vec3], np.ndarray]:
        """
        Get the atomic positions.

        Parameters
        ----------
        asNumpy : bool=False
            If true, the values are returned as a numpy array instead of a list of Vec3s.

        Returns
        -------
        list of Vec3 or np.ndarray
            The atomic positions
        """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = np.array(self.positions.value_in_unit(nanometers)) * nanometers
            return self._numpyPositions
        return self.positions

    def getBoxVectors(self, asNumpy: bool = False) -> Union[List[Vec3], np.ndarray]:
        """
        Get the periodic box vectors.

        Parameters
        ----------
        asNumpy : bool, optional, default=False
            If true, the values are returned as a numpy array instead of a list of Vec3s.

        Returns
        -------
        list of Vec3 or np.ndarray
            The periodic box vectors.
        """
        if self.boxVectors is None:
            raise AttributeError(f"Box information not found in {self.file}")
        if asNumpy:
            if self._numpyBoxVectors is None:
                self._numpyBoxVectors = []
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[0].value_in_unit(nanometers)) * nanometers
                    )
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[1].value_in_unit(nanometers)) * nanometers
                    )
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[2].value_in_unit(nanometers)) * nanometers
                )
            return self._numpyBoxVectors
        return self.boxVectors
    

    def _loadSeqFile(self, file: str) -> Dict[str, List[str]]:
        """
        Load the biopolymer sequence from a TINKER .seq file. 

        Parameters
        ----------
        file : str
            The path to the .seq file.

        Notes
        -----
        This function is heavily based on the TINKER readseq.f routine.

        Tinker .seq files list residues belonging to one or more chains.
        Each line typically starts with a residue number, optionally preceded
        by a single-character chain ID. Lines starting with residue number 1
        indicate the beginning of a new chain. If no chain ID is provided
        for the first residue, a default ID (A, B, C, ...) is assigned.

        With Chain ID:
            A     1  CYS PHE GLU PRO PRO PRO ALA THR THR THR GLN THR GLY PHE ARG
            A    16  LEU ...

        Without Chain ID (Chain 'A' assumed for first block, 'B' for second, etc.):
            1  CYS PHE GLU PRO PRO PRO ALA THR THR THR GLN THR GLY PHE ARG
           16  LEU ...
            1  MET VAL ...  (Starts a new chain, 'B')


        Returns
        -------
        dict[str, list[str]]
            A dictionary where keys are chain IDs and values are lists of
            three-letter residue codes for that chain.
        """
        seqDict = {}
        defaultChainIndex = 0
        currentChainId = None 

        if os.path.exists(file):
            with open(file, "r") as f:
                seqContent = f.readlines()
        else:
            raise FileNotFoundError(f"File {file} not found")

        for line_num, line in enumerate(seqContent, 1):
            lineSplit = line.split()
            if not lineSplit:
                continue  

            parsedChainId = None
            startResidue = None
            residuesStartIndex = -1

            # Determine line format: Chain ID + Res Num or just Res Num
            if lineSplit[0].isalpha() and len(lineSplit[0]) == 1 and lineSplit[1].isdigit():
                parsedChainId = lineSplit[0]
                startResidue = int(lineSplit[1])
                residuesStartIndex = 2
            elif lineSplit[0].isdigit():
                startResidue = int(lineSplit[0])
                residuesStartIndex = 1
            else:
                raise ValueError(
                    f"Line {line_num}: Does not have the expected format "
                    f"(ChainID ResNum ... or ResNum ...): {line.strip()}"
                )
            
            # Force 3-letter residues, otherwise add whitespaces to the left
            residues = [residue.rjust(3) for residue in lineSplit[residuesStartIndex:]]
            
            if not residues:
                raise ValueError(f"Line {line_num}: No residues found after parsing!")

            if startResidue == 1:
                if parsedChainId:
                    if parsedChainId in seqDict:
                        raise ValueError(
                            f"Line {line_num}: Duplicate start definition (residue 1) "
                            f"for chain ID '{parsedChainId}'"
                        )
                    currentChainId = parsedChainId
                else:
                    currentChainId = chr(65 + (defaultChainIndex % 26))
                    defaultChainIndex += 1

                seqDict[currentChainId] = {
                    "residues": [],
                    "chainType": "GENERIC"
                }

            else: 
                if currentChainId is None:
                     raise ValueError(
                        f"Line {line_num}: Continuation line encountered before "
                        f"any chain was started (residue 1)."
                    )

                if parsedChainId is not None and parsedChainId != currentChainId:
                    raise ValueError(
                        f"Line {line_num}: Chain ID '{parsedChainId}' provided for residue "
                        f"{startResidue}, but expected continuation of chain '{currentChainId}'"
                    )

                expectedResidueNum = len(seqDict[currentChainId]["residues"]) + 1
                if startResidue != expectedResidueNum:
                    raise ValueError(
                        f"Line {line_num}: Residue numbering inconsistency for chain "
                        f"'{currentChainId}'. Expected residue {expectedResidueNum}, "
                        f"but line starts with {startResidue}."
                    )

            # Add residues to the current chain
            seqDict[currentChainId]["residues"].extend(residues)

        # Determine chain types
        for chainId, chainData in seqDict.items():
            residues = chainData["residues"]
            aminoAcidIntersect = set(residues) & set(self._AMINO_ACID_ORDER)
            if aminoAcidIntersect:
                chainData["chainType"] = "PEPTIDE"
                continue

            nucleotideIntersect = set(residues) & set(self._NUCLEOTIDE_ORDER)
            if nucleotideIntersect:
                chainData["chainType"] = "NUCLEIC"
                continue

        return seqDict

    # ------------------------------------------------------------------------------------------ #
    #                                   XYZ FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _parseAndStoreXyzLine(line: str, xyzDict: dict) -> None:
        """
        Parse a line of an .xyz file and store the data in a dictionary.

        Parameters
        ----------
        line : str
            The line containing atom data.
        xyzDict : dict
            The dictionary to store parsed atom data.
        """
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                "Each line in the TINKER .xyz file must have at least 6 fields"
            )

        index = int(fields[0]) - 1
        symbol = str(fields[1])
        x = float(fields[2]) * 0.1
        y = float(fields[3]) * 0.1
        z = float(fields[4]) * 0.1
        position = Vec3(x, y, z)
        atomType = str(fields[5])
        bonds = [int(bond) - 1 for bond in fields[6:]]

        xyzDict[index] = {
            "symbol": symbol,
            "positions": position,
            "atomType": atomType,
            "bonds": bonds,
        }

    def _loadXyzFile(self, file: str) -> Tuple[Dict, List[Vec3], List[Vec3]]:
        """
        Load a TINKER .xyz file.

        Parameters
        ----------
        file : str
            The name of the .xyz file to load.

        Returns
        -------
        self._xyzDict, self.boxVectors, self.positions : dict, list of Vec3, list of Vec3
            The dictionary with the atom data, the box vectors, and the atomic positions.
        """
        try:
            self._xyzDict = {}

            with open(file, "r") as f:
                # Read number of atoms
                nAtoms = int(f.readline().split()[0])
                linesLeft = nAtoms

                # Read the second line
                secondLine = f.readline().strip()
                secondLineSplit = secondLine.split()

                # Check for box information
                if len(secondLineSplit) == 6 and secondLineSplit[0] != "1":
                    # Read box unit vectors and angles from the second line
                    box = [
                        float(val) * (0.1 if i < 3 else math.pi / 180.0)
                        for i, val in enumerate(secondLineSplit)
                    ]
                    self.boxVectors = computePeriodicBoxVectors(*box)
                else:
                    # No box information, so treat the second line as atom positions
                    TinkerFiles._parseAndStoreXyzLine(secondLine, self._xyzDict)
                    linesLeft -= 1

                # Process the remaining atom lines
                for _ in range(linesLeft):
                    atomLine = f.readline().strip()
                    TinkerFiles._parseAndStoreXyzLine(atomLine, self._xyzDict)

            # Store the positions
            self.positions = [self._xyzDict[i]["positions"] for i in range(nAtoms)]

        except FileNotFoundError:
            raise IOError(f"Could not find file {file}")
        except Exception as e:
            raise ValueError(f"Error parsing {file}: {e}")

        return self._xyzDict, self.boxVectors, self.positions

    # ------------------------------------------------------------------------------------------ #
    #                                   KEY FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _addAtomType(
        atomTypeDict: Dict,
        atomType: str,
        atomClass: str,
        nameShort: str,
        nameLong: str,
        atomicNumber: int,
        mass: float,
        valence: int,
        element: str,
        residue: str,
    ) -> None:
        """
        Helper function to validate atom type information.

        Parameters
        ----------
        atomTypeDict : dict
            The dictionary of atom type data.
        atomType : str
            The atom type.
        atomClass : str
            The atom class.
        nameShort : str
            The short name of the atom type.
        nameLong : str
            The long name of the atom type.
        atomicNumber : int
            The atomic number of the atom type.
        mass : float
            The mass of the atom type.
        valence : int
            The valence of the atom type.
        element : str
            The element of the atom type.
        residue : str
            The residue of the atom type.
        """
        if atomType in atomTypeDict:
            # Validate against existing atom type data
            stored = atomTypeDict[atomType]
            mismatches = {
                "atomClass": atomClass,
                "nameShort": nameShort,
                "nameLong": nameLong,
                "element": element,
                "atomicNumber": atomicNumber,
                "mass": mass,
                "valence": valence,
                "residue": residue,
            }

            for key, new_value in mismatches.items():
                if stored[key] != new_value:
                    raise ValueError(
                        f"Inconsistent {key} for atom type '{atomType}': "
                        f"expected '{stored[key]}', got '{new_value}'."
                    )

        # Add new atom type to the dictionary
        atomTypeDict[atomType] = {
            "atomClass": atomClass,
            "nameShort": nameShort,
            "nameLong": nameLong,
            "element": element,
            "atomicNumber": atomicNumber,
            "mass": mass,
            "valence": valence,
            "residue": residue,
        }

    @staticmethod
    def _addBioType(
        bioTypeDict: Dict,
        bioType: str,
        nameShort: str,
        nameLong: str,
        atomType: str,
        element: str,
        residue: str,
    ) -> None:
        """
        Helper function to validate biotype information.

        Parameters
        ----------
        bioTypeDict : dict
            The dictionary of biotype data.
        bioType : str
            The bio type.
        nameShort : str
            The short name of the biotype.
        nameLong : str
            The long name of the biotype.
        atomType : str
            The atom type counterpart of the biotype.
        element : str
            The element of the biotype.
        residue : str
            The residue of the biotype
        """
        if bioType in bioTypeDict:
            # Validate against existing atom type data
            stored = bioTypeDict[bioType]
            mismatches = {
                "nameShort": nameShort,
                "nameLong": nameLong,
                "atomType": atomType,
                "element": element,
                "residue": residue,
            }

            for key, new_value in mismatches.items():
                if stored[key] != new_value:
                    raise ValueError(
                        f"Inconsistent {key} for biotype '{bioType}': "
                        f"expected '{stored[key]}', got '{new_value}'."
                    )

        # Add new biotype to the dictionary
        bioTypeDict[bioType] = {
            "nameShort": nameShort,
            "nameLong": nameLong,
            "atomType": atomType,
            "element": element,
            "residue": residue,
        }

    @staticmethod
    def getResidueAbbreviation(residueName: str) -> str:
        """
        Get the residue abbreviation for a given residue name.

        Parameters
        ----------
        residueName : str
            The name of the residue.

        Returns
        -------
        str
            The residue abbreviation.
        """
        for abbr, data in TinkerFiles.RESIDUE_MAPPING.items():
            if data["tinkerLookupName"] in residueName:
                if data["type"] == "ion":
                    return abbr
                elif data["type"] == "AmoebaWater":
                    return "HOH"

        return residueName[:3].upper()

    def _loadKeyFile(self, keyFile: str) -> Tuple[Dict, Dict, Dict, Dict]:
        """
        Load a TINKER .key or .prm file.

        Parameters
        ----------
        keyFile : str
            The name of the file.

        Returns
        -------
        atomTypes, bioTypes, forces, scalars : dict, dict, dict, dict
            The atom types, bio types, forces, and scalars.
        """
        # Get all interesting lines from the file
        try:
            allLines = []
            with open(keyFile) as file:
                for line in file:
                    try:
                        if line.count('"') % 2 != 0:
                            # Skip lines with an odd number of quotes to avoid parsing errors
                            # with citations or other non-essential information
                            continue

                        fields = shlex.split(line)

                        if not fields or fields[0].startswith("#"):
                            # Skip empty lines and comments
                            continue
                        allLines.append(fields)
                    except Exception as e:
                        raise ValueError(f"Error parsing line in {keyFile}: {e}")
        except FileNotFoundError:
            raise IOError(f"Could not find file {keyFile}")
        except Exception as e:
            raise ValueError(f"Error reading {keyFile}: {e}")

        atomTypesDict = dict()
        bioTypesDict = dict()
        forcesDict = dict()
        # We use the default values for the scalars as a starting point
        # and update them with the values from the key file
        scalarsDict = self.RECOGNIZED_SCALARS.copy()

        lineIndex = 0
        while lineIndex < len(allLines):
            fields = allLines[lineIndex]
            if fields[0] == "atom":
                if len(fields) != 8:
                    raise ValueError(
                        f"Invalid atom line: Expected 8 fields, got {len(fields)}. Fields: {fields}"
                    )
                # Atom type information
                # atom atomType atomClass nameShort nameLong atomicNumber mass valence
                (
                    atomType,
                    atomClass,
                    nameShort,
                    nameLong,
                    atomicNumber,
                    mass,
                    valence,
                ) = fields[1:]
                element = elem.Element.getByAtomicNumber(int(atomicNumber)).symbol
                nameLong = re.sub(r"\s+", " ", nameLong.strip())
                resAbbr = TinkerFiles.getResidueAbbreviation(nameLong)
                TinkerFiles._addAtomType(
                    atomTypesDict,
                    atomType,
                    atomClass,
                    nameShort,
                    nameLong,
                    atomicNumber,
                    mass,
                    valence,
                    element,
                    resAbbr,
                )
                lineIndex += 1
            elif fields[0] == "biotype":
                # Biotype information
                if len(fields) != 5:
                    raise ValueError(
                        f"Invalid biotype line: Expected 5 fields, got {len(fields)}. Fields: {fields}"
                    )
                bioType, nameShort, nameLong, atomType = fields[1:]
                element = elem.Element.getByAtomicNumber(int(atomicNumber)).symbol
                nameLong = nameLong.replace('"', "")
                lookUp = f"{nameShort}_{nameLong}"
                if lookUp in bioTypesDict:
                    # Workaround for Tinker using the same name but different types for H2', H2'', and for H5', H5''
                    lookUp = f"{nameShort}*_{nameLong}"
                resAbbr = TinkerFiles.getResidueAbbreviation(nameLong)
                TinkerFiles._addBioType(
                    bioTypesDict,
                    bioType,
                    nameShort,
                    nameLong,
                    atomType,
                    element,
                    resAbbr,
                )
                lineIndex += 1
            elif fields[0] in self.RECOGNIZED_FORCES:
                if self.RECOGNIZED_FORCES[fields[0]] == 1:
                    if fields[0] not in forcesDict:
                        forcesDict[fields[0]] = []
                    forcesDict[fields[0]].append(fields[1:])
                    lineIndex += 1
                else:
                    # Call the function to parse the specific force
                    lineIndex = self.RECOGNIZED_FORCES[fields[0]](
                        lineIndex, allLines, forcesDict
                    )
            elif fields[0] in self.RECOGNIZED_SCALARS:
                scalar, value = fields
                scalarsDict[scalar] = value
                lineIndex += 1
            else:
                # Skip unrecognized fields
                lineIndex += 1

        return atomTypesDict, bioTypesDict, forcesDict, scalarsDict

    @staticmethod
    def __addMultipole(
        lineIndex: int, allLines: List[List[str]], forces: Dict[str, Any]
    ) -> int:
        """
        Parse and store multipole force data from the key file.

        Parameters
        ----------
        lineIndex : int
            The current line index in the key file.
        allLines : list of list of str
            All lines from the key file, split into fields.
        forces : dict
            The forces dictionary to store the parsed data.

        Returns
        -------
        int
            The updated line index after parsing the multipole force data.
        """
        if "multipole" not in forces:
            forces["multipole"] = []
        fields = allLines[lineIndex]
        multipoles = [fields[-1]]
        axisInfo = fields[1:-1]
        lineIndex += 1
        fields = allLines[lineIndex]
        multipoles.append(fields[0])
        multipoles.append(fields[1])
        multipoles.append(fields[2])
        lineIndex += 1
        fields = allLines[lineIndex]
        multipoles.append(fields[0])
        lineIndex += 1
        fields = allLines[lineIndex]
        multipoles.append(fields[0])
        multipoles.append(fields[1])
        lineIndex += 1
        fields = allLines[lineIndex]
        multipoles.append(fields[0])
        multipoles.append(fields[1])
        multipoles.append(fields[2])
        lineIndex += 1
        multipoleInfo = [axisInfo, multipoles]
        forces["multipole"].append(multipoleInfo)
        return lineIndex

    @staticmethod
    def __addTorTor(
        lineIndex: int, allLines: List[List[str]], forces: Dict[str, Any]
    ) -> int:
        """
        Parse and store torsion-torsion force data from the key file.

        Parameters
        ----------
        lineIndex : int
            The current line index in the key file.
        allLines : list of list of str
            All lines from the key file, split into fields.
        forces : dict
            The forces dictionary to store the parsed data.

        Returns
        -------
        int
            The updated line index after parsing the torsion-torsion force data.
        """
        if "tortors" not in forces:
            forces["tortors"] = []
        fields = allLines[lineIndex]
        tortorInfo = fields[1:]
        lastGridLine = lineIndex + int(fields[6]) * int(fields[7])
        grid = []
        while lineIndex < lastGridLine:
            lineIndex += 1
            grid.append(allLines[lineIndex])
        forces["tortors"].append([tortorInfo, grid])
        return lineIndex

    @staticmethod
    def _combineXyzAndKeyData(xyzDict: Dict, atomTypes: List, bioTypes: List) -> Dict:
        """
        Combine the data from the .xyz and .key files into one dictionary.

        Parameters
        ----------
        xyzDict : dict
            The data from the .xyz file.
        atomTypes : dict
            The atom types from the .key file.
        bioTypes : dict
            The bio types from the .key file.

        Returns
        -------
        dict
            The combined data.

        Notes
        -----
        atomData[index]                          index = 0, 1, 2, ...
        atomData[index]['symbol']                atom symbol
        atomData[index]['positions']             atom position
        atomData[index]['bonds']                 list of bonded atom indices
        atomData[index]['residue']               residue name

        # if atomType is present in the .key file
        atomData[index]['atomType']              atom type
        atomData[index]['atomClass']             atom class
        atomData[index]['nameShort']             short name
        atomData[index]['nameLong']              long name
        atomData[index]['element']               element
        atomData[index]['atomicNumber']          atomic number
        atomData[index]['mass']                  mass
        atomData[index]['valence']               valence

        # if bioType is present in the .key file
        atomData[index]['bioType']               bio type
        atomData[index]['nameShort']             short name
        atomData[index]['nameLong']              long name
        atomData[index]['atomType']              atom type
        atomData[index]['element']               element
        """
        atomData = dict()
        for atomIndex in xyzDict:
            atomData[atomIndex] = xyzDict[atomIndex]

            if "atomType" in atomData[atomIndex]:
                # Add all the atom type data to the atom data
                for atomTypeDict in atomTypes:
                    if atomData[atomIndex]["atomType"] in atomTypeDict:
                        atomData[atomIndex].update(
                            atomTypeDict[atomData[atomIndex]["atomType"]]
                        )
                        break
            elif "bioType" in atomData[atomIndex]:
                # Add all the biotype data to the atom data
                for bioTypeDict in bioTypes:
                    if atomData[atomIndex]["bioType"] in bioTypeDict:
                        atomData[atomIndex].update(
                            bioTypeDict[atomData[atomIndex]["bioType"]]
                        )
                        break

        return atomData

    # ------------------------------------------------------------------------------------------ #
    #                                        TOPOLOGY                                            #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _getResiduesFromBonds(atomData: Dict) -> Tuple[List[List[int]], List[str]]:
        """
        Form whole residues by traversing bonded atoms.

        Notes
        -----
        Residues are defined as connected atoms in the graph formed by the bonds between atoms.
        For example, connected amino acids in a protein are considered to be part of the same residue.

        Returns
        -------
        list of list of int
            A list of residues, where each residue is a list of atom indices.
        """

        def _getResidueName(atomIndex):
            # Find atom type in atomTypesDict
            atomType = atomData[atomIndex]["atomType"]
            residueName = atomData[atomIndex]["residue"]

            nameParts = residueName.split()
            if len(nameParts) == 2:
                residueName = nameParts[0]
            else:
                residueName = " ".join(nameParts[:2])
            return residueName

        residues = []
        residueNames = []
        seen = set()

        for atom1 in range(len(atomData)):
            if atom1 not in seen:
                # Start a new residue
                residue = []

                # Get the residue name
                residueName = _getResidueName(atom1)

                # Iterative way of traversing a graph-like structure
                # Recursive approach led to stack overflow
                stack = [atom1]
                while stack:
                    atom = stack.pop()
                    if atom not in seen:
                        seen.add(atom)
                        residue.append(atom)
                        stack.extend(atomData[atom]["bonds"])

                residues.append(sorted(residue))
                residueNames.append(residueName)

        return residues, residueNames

    def _createTopology(self, atomData: Dict) -> top.Topology:
        """
        Build the topology from the data parsed from the .xyz file and the .key file.

        Parameters
        ----------
        atomData : dict
            The atom data parsed from the .xyz and .key files.

        Returns
        -------
        openmm.app.topology.Topology
            The topology object.
        """
        self.topology = top.Topology()

        # Infer the residues from the bonds
        residues, residueNames = self._getResiduesFromBonds(atomData)

        # Add chain to the topology
        chain = self.topology.addChain()
        for residueName, residueAtoms in zip(residueNames, residues):
            # Add residues to the topology
            residue = self.topology.addResidue(residueName, chain)
            for atomIndex in residueAtoms:
                name = atomData[atomIndex]["nameShort"]
                element = elem.Element.getByAtomicNumber(
                    int(atomData[atomIndex]["atomicNumber"])
                )

                # Add atoms to the topology
                self.topology.addAtom(name, element, residue)

        # Add bonds to the topology
        atoms = list(self.topology.atoms())
        seenBonds = set()
        for atomIndex in range(len(atomData)):
            for bondedAtomIndex in atomData[atomIndex]["bonds"]:
                if (bondedAtomIndex, atomIndex) not in seenBonds:
                    self.topology.addBond(atoms[atomIndex], atoms[bondedAtomIndex])
                    seenBonds.add((atomIndex, bondedAtomIndex))

        return self.topology

    # ---------------------------------------------------------------------------------------- #
    #                                      WRITE XML FILE                                      #
    # ---------------------------------------------------------------------------------------- #
    @staticmethod
    def _writeAtomTypes(
        root: etree.Element, atomTypes: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Write the atom types to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomTypes : dict
            The atom types dictionary.
        """
        atomTypesElement = etree.SubElement(root, "AtomTypes")
        for atomType in sorted(atomTypes, key=int):
            atomTypeElement = etree.SubElement(atomTypesElement, "Type")
            atomTypeElement.attrib["name"] = atomType
            atomTypeElement.attrib["class"] = atomTypes[atomType]["atomClass"]
            atomTypeElement.attrib["element"] = atomTypes[atomType]["element"]
            atomTypeElement.attrib["mass"] = atomTypes[atomType]["mass"]

    @staticmethod
    def _writeResidues(
        root: etree.Element,
        atomData: Dict[int, Dict[str, Any]],
        topology: top.Topology,
        residuesSet: Set[str],
    ) -> None:
        """
        Write the residues to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomData : dict
            The atom data dictionary.
        topology : openmm.app.topology.Topology
            The topology object.
        """
        residuesElement = etree.SubElement(root, "Residues")
        residuesSeen = set()
        for residue in topology.residues():
            if residue.name in residuesSet and residue.name not in residuesSeen:
                residueElement = etree.SubElement(residuesElement, "Residue")
                residueElement.attrib["name"] = residue.name

                if residue.name == "HOH":
                    # AMOEBA water is flexible
                    residueElement.attrib["rigidWater"] = "false"

                for i, atom in enumerate(residue.atoms()):
                    atomElement = etree.SubElement(residueElement, "Atom")
                    atomElement.attrib["name"] = atomData[atom.index][
                        "nameShort"
                    ] + str(i)
                    atomElement.attrib["type"] = atomData[atom.index]["atomType"]

                baseIndex = next(atom.index for atom in residue.atoms())
                for bond in residue.bonds():
                    bondElement = etree.SubElement(residueElement, "Bond")
                    bondElement.attrib["from"] = str(bond[0].index - baseIndex)
                    bondElement.attrib["to"] = str(bond[1].index - baseIndex)

                residuesSeen.add(residue.name)

    @staticmethod
    def _writeAmoebaBondForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaBondForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "bond" not in forces:
            return

        cubic = 10.0 * float(scalars["bond-cubic"])
        quartic = 100.0 * float(scalars["bond-quartic"])
        bondForceElement = etree.SubElement(root, "AmoebaBondForce")
        bondForceElement.attrib["bond-cubic"] = str(cubic)
        bondForceElement.attrib["bond-quartic"] = str(quartic)

        for bond in forces["bond"]:
            length = float(bond[3]) * 0.1
            k = float(bond[2]) * 100.0 * 4.184
            bondElement = etree.SubElement(bondForceElement, "Bond")
            bondElement.attrib["class1"] = bond[0]
            bondElement.attrib["class2"] = bond[1]
            bondElement.attrib["length"] = str(length)
            bondElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaAngleForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaAngleForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "angle" not in forces and "anglep" not in forces:
            return

        cubic = float(scalars["angle-cubic"])
        quartic = float(scalars["angle-quartic"])
        pentic = float(scalars["angle-pentic"])
        sextic = float(scalars["angle-sextic"])
        angleForceElement = etree.SubElement(root, "AmoebaAngleForce")
        angleForceElement.attrib["angle-cubic"] = str(cubic)
        angleForceElement.attrib["angle-quartic"] = str(quartic)
        angleForceElement.attrib["angle-pentic"] = str(pentic)
        angleForceElement.attrib["angle-sextic"] = str(sextic)

        radian = 57.2957795130
        radian2 = 4.184 / (radian * radian)

        for angleSet in ["angle", "anglep"]:
            if angleSet not in forces:
                continue
            for angle in forces[angleSet]:
                k = float(angle[3]) * radian2
                angleElement = etree.SubElement(angleForceElement, "Angle")
                angleElement.attrib["class1"] = angle[0]
                angleElement.attrib["class2"] = angle[1]
                angleElement.attrib["class3"] = angle[2]
                angleElement.attrib["k"] = str(k)
                angleElement.attrib["angle1"] = angle[4]
                if len(angle) > 5:
                    angleElement.attrib["angle2"] = angle[5]
                if len(angle) > 6:
                    angleElement.attrib["angle3"] = angle[6]
                angleElement.attrib["inPlane"] = (
                    "true" if angleSet == "anglep" else "false"
                )

    @staticmethod
    def _writeAmoebaOutOfPlaneBendForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaOutOfPlaneBendForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "opbend" not in forces:
            return

        cubic = float(scalars["opbend-cubic"])
        quartic = float(scalars["opbend-quartic"])
        pentic = float(scalars["opbend-pentic"])
        sextic = float(scalars["opbend-sextic"])
        opbendType = scalars["opbendtype"]

        opbendForceElement = etree.SubElement(root, "AmoebaOutOfPlaneBendForce")
        opbendForceElement.attrib["type"] = opbendType
        opbendForceElement.attrib["opbend-cubic"] = str(cubic)
        opbendForceElement.attrib["opbend-quartic"] = str(quartic)
        opbendForceElement.attrib["opbend-pentic"] = str(pentic)
        opbendForceElement.attrib["opbend-sextic"] = str(sextic)

        radian = 57.2957795130
        radian2 = 4.184 / (radian * radian)

        for opbend in forces["opbend"]:
            k = float(opbend[4]) * radian2
            angleElement = etree.SubElement(opbendForceElement, "Angle")
            angleElement.attrib["class1"] = opbend[0]
            angleElement.attrib["class2"] = opbend[1]
            angleElement.attrib["class3"] = opbend[2]
            angleElement.attrib["class4"] = opbend[3]
            angleElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaTorsionForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "torsion" not in forces:
            return

        torsionUnit = float(scalars["torsionunit"])
        conversion = 4.184 * torsionUnit
        radian = 57.2957795130

        torsionForceElement = etree.SubElement(root, "PeriodicTorsionForce")
        for torsion in forces["torsion"]:
            torsionElement = etree.SubElement(torsionForceElement, "Proper")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            startIndex = 4
            for ii in range(0, 3):
                torsionSuffix = str(ii + 1)
                amplitudeAttributeName = f"k{torsionSuffix}"
                angleAttributeName = f"phase{torsionSuffix}"
                periodicityAttributeName = f"periodicity{torsionSuffix}"
                amplitude = float(torsion[startIndex]) * conversion
                angle = float(torsion[startIndex + 1]) / radian
                periodicity = int(torsion[startIndex + 2])
                torsionElement.attrib[amplitudeAttributeName] = str(amplitude)
                torsionElement.attrib[angleAttributeName] = str(angle)
                torsionElement.attrib[periodicityAttributeName] = str(periodicity)
                startIndex += 3

    @staticmethod
    def _writeAmoebaPiTorsionForce(root: etree.Element, forces: Dict[str, Any]) -> None:
        """
        Write the AmoebaPiTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "pitors" not in forces:
            return

        piTorsionUnit = 1.0
        conversion = 4.184 * piTorsionUnit
        piTorsionForceElement = etree.SubElement(root, "AmoebaPiTorsionForce")
        piTorsionForceElement.attrib["piTorsionUnit"] = str(piTorsionUnit)
        for piTorsion in forces["pitors"]:
            k = float(piTorsion[2]) * conversion
            piTorsionElement = etree.SubElement(piTorsionForceElement, "PiTorsion")
            piTorsionElement.attrib["class1"] = piTorsion[0]
            piTorsionElement.attrib["class2"] = piTorsion[1]
            piTorsionElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaStretchTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaStretchTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "strtors" not in forces:
            return

        stretchTorsionForceElement = etree.SubElement(root, "AmoebaStretchTorsionForce")
        for torsion in forces["strtors"]:
            v = [float(x) * 10 * 4.184 for x in torsion[4:]]
            torsionElement = etree.SubElement(stretchTorsionForceElement, "Torsion")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            torsionElement.attrib["v11"] = str(v[0])
            torsionElement.attrib["v12"] = str(v[1])
            torsionElement.attrib["v13"] = str(v[2])
            torsionElement.attrib["v21"] = str(v[3])
            torsionElement.attrib["v22"] = str(v[4])
            torsionElement.attrib["v23"] = str(v[5])
            torsionElement.attrib["v31"] = str(v[6])
            torsionElement.attrib["v32"] = str(v[7])
            torsionElement.attrib["v33"] = str(v[8])

    @staticmethod
    def _writeAmoebaAngleTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaAngleTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "angtors" not in forces:
            return

        angleTorsionForceElement = etree.SubElement(root, "AmoebaAngleTorsionForce")
        for torsion in forces["angtors"]:
            v = [float(x) * 4.184 for x in torsion[4:]]
            torsionElement = etree.SubElement(angleTorsionForceElement, "Torsion")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            torsionElement.attrib["v11"] = str(v[0])
            torsionElement.attrib["v12"] = str(v[1])
            torsionElement.attrib["v13"] = str(v[2])
            torsionElement.attrib["v21"] = str(v[3])
            torsionElement.attrib["v22"] = str(v[4])
            torsionElement.attrib["v23"] = str(v[5])

    @staticmethod
    def _writeAmoebaStretchBendForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaStretchBendForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "strbnd" not in forces:
            return

        stretchBendUnit = 1.0
        radian = 57.2957795130
        conversion = 41.84 / radian
        stretchBendForceElement = etree.SubElement(root, "AmoebaStretchBendForce")
        stretchBendForceElement.attrib["stretchBendUnit"] = str(stretchBendUnit)

        for stretchBend in forces["strbnd"]:
            k1 = float(stretchBend[3]) * conversion
            k2 = float(stretchBend[4]) * conversion
            stretchBendElement = etree.SubElement(
                stretchBendForceElement, "StretchBend"
            )
            stretchBendElement.attrib["class1"] = stretchBend[0]
            stretchBendElement.attrib["class2"] = stretchBend[1]
            stretchBendElement.attrib["class3"] = stretchBend[2]
            stretchBendElement.attrib["k1"] = str(k1)
            stretchBendElement.attrib["k2"] = str(k2)

    @staticmethod
    def _writeAmoebaTorsionTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaTorsionTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "tortors" not in forces:
            return

        torsionTorsionUnit = 1.0

        torsionTorsionForceElement = etree.SubElement(root, "AmoebaTorsionTorsionForce")
        for index, torsionTorsion in enumerate(forces["tortors"]):
            torsionTorsionElement = etree.SubElement(
                torsionTorsionForceElement, "TorsionTorsion"
            )
            torsionTorsionElement.attrib["class1"] = torsionTorsion[0][0]
            torsionTorsionElement.attrib["class2"] = torsionTorsion[0][1]
            torsionTorsionElement.attrib["class3"] = torsionTorsion[0][2]
            torsionTorsionElement.attrib["class4"] = torsionTorsion[0][3]
            torsionTorsionElement.attrib["class5"] = torsionTorsion[0][4]
            torsionTorsionElement.attrib["grid"] = str(index)
            torsionTorsionElement.attrib["nx"] = str(torsionTorsion[0][5])
            torsionTorsionElement.attrib["ny"] = str(torsionTorsion[0][6])

        for index, torsionTorsion in enumerate(forces["tortors"]):
            gridElement = etree.SubElement(
                torsionTorsionForceElement, "TorsionTorsionGrid"
            )
            gridElement.attrib["grid"] = str(index)
            gridElement.attrib["nx"] = str(torsionTorsion[0][5])
            gridElement.attrib["ny"] = str(torsionTorsion[0][6])
            for gridIndex, gridEntry in enumerate(torsionTorsion[1]):
                if len(gridEntry) > 5:
                    f = float(gridEntry[2]) * 4.184
                    fx = float(gridEntry[3]) * 4.184
                    fy = float(gridEntry[4]) * 4.184
                    fxy = float(gridEntry[5]) * 4.184
                    gridEntryElement = etree.SubElement(gridElement, "Grid")
                    gridEntryElement.attrib["angle1"] = gridEntry[0]
                    gridEntryElement.attrib["angle2"] = gridEntry[1]
                    gridEntryElement.attrib["f"] = str(f)
                    gridEntryElement.attrib["fx"] = str(fx)
                    gridEntryElement.attrib["fy"] = str(fy)
                    gridEntryElement.attrib["fxy"] = str(fxy)
                elif len(gridEntry) > 2:
                    f = float(gridEntry[2]) * 4.184
                    gridEntryElement = etree.SubElement(gridElement, "Grid")
                    gridEntryElement.attrib["angle1"] = gridEntry[0]
                    gridEntryElement.attrib["angle2"] = gridEntry[1]
                    gridEntryElement.attrib["f"] = str(f)

    @staticmethod
    def _writeAmoebaVdwForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaVdwForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "vdw" not in forces and "vdwpr" not in forces:
            return

        vdwForceElement = etree.SubElement(root, "AmoebaVdwForce")
        vdwForceElement.attrib["type"] = scalars["vdwtype"]
        vdwForceElement.attrib["radiusrule"] = scalars["radiusrule"]
        vdwForceElement.attrib["radiustype"] = scalars["radiustype"]
        vdwForceElement.attrib["radiussize"] = scalars["radiussize"]
        vdwForceElement.attrib["epsilonrule"] = scalars["epsilonrule"]
        vdwForceElement.attrib["vdw-13-scale"] = scalars["vdw-13-scale"]
        vdwForceElement.attrib["vdw-14-scale"] = scalars["vdw-14-scale"]
        vdwForceElement.attrib["vdw-15-scale"] = scalars["vdw-15-scale"]

        if "vdw" in forces:
            for vdw in forces["vdw"]:
                sigma = float(vdw[1]) * 0.1
                epsilon = float(vdw[2]) * 4.184
                reduction = vdw[3] if len(vdw) > 3 else 1.0
                vdwElement = etree.SubElement(vdwForceElement, "Vdw")
                vdwElement.attrib["class"] = vdw[0]
                vdwElement.attrib["sigma"] = str(sigma)
                vdwElement.attrib["epsilon"] = str(epsilon)
                vdwElement.attrib["reduction"] = str(reduction)

        if "vdwpr" in forces:
            for pair in forces["vdwpr"]:
                sigma = float(pair[2]) * 0.1
                epsilon = float(pair[3]) * 4.184
                pairElement = etree.SubElement(vdwForceElement, "Pair")
                pairElement.attrib["class1"] = pair[0]
                pairElement.attrib["class2"] = pair[1]
                pairElement.attrib["sigma"] = str(sigma)
                pairElement.attrib["epsilon"] = str(epsilon)

    @staticmethod
    def _writeAmoebaMultipoleForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaMultipoleForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "multipole" not in forces and "polarize" not in forces:
            return

        amoebaMultipoleForceElement = etree.SubElement(root, "AmoebaMultipoleForce")
        amoebaMultipoleForceElement.attrib["mpole12Scale"] = scalars["mpole-12-scale"]
        amoebaMultipoleForceElement.attrib["mpole13Scale"] = scalars["mpole-13-scale"]
        amoebaMultipoleForceElement.attrib["mpole14Scale"] = scalars["mpole-14-scale"]
        amoebaMultipoleForceElement.attrib["mpole15Scale"] = scalars["mpole-15-scale"]
        amoebaMultipoleForceElement.attrib["polar12Scale"] = scalars["polar-12-scale"]
        amoebaMultipoleForceElement.attrib["polar13Scale"] = scalars["polar-13-scale"]
        amoebaMultipoleForceElement.attrib["polar14Scale"] = scalars["polar-14-scale"]
        amoebaMultipoleForceElement.attrib["polar15Scale"] = scalars["polar-15-scale"]
        amoebaMultipoleForceElement.attrib["polar14Intra"] = scalars["polar-14-intra"]
        amoebaMultipoleForceElement.attrib["direct11Scale"] = scalars["direct-11-scale"]
        amoebaMultipoleForceElement.attrib["direct12Scale"] = scalars["direct-12-scale"]
        amoebaMultipoleForceElement.attrib["direct13Scale"] = scalars["direct-13-scale"]
        amoebaMultipoleForceElement.attrib["direct14Scale"] = scalars["direct-14-scale"]
        amoebaMultipoleForceElement.attrib["mutual11Scale"] = scalars["mutual-11-scale"]
        amoebaMultipoleForceElement.attrib["mutual12Scale"] = scalars["mutual-12-scale"]
        amoebaMultipoleForceElement.attrib["mutual13Scale"] = scalars["mutual-13-scale"]
        amoebaMultipoleForceElement.attrib["mutual14Scale"] = scalars["mutual-14-scale"]

        bohr = 0.52917720859
        dipoleConversion = 0.1 * bohr
        quadrupoleConversion = 0.01 * bohr * bohr / 3.0
        for multipoleInfo in forces["multipole"]:
            axisInfo = multipoleInfo[0]

            if int(axisInfo[0]) < 0:
                # Some Tinker files include a multipole parameter definition
                # for atom types with negative indices with no multipole parameters (all zeros).
                # We skip these entries as they cause an error in OpenMM.
                continue

            multipoles = multipoleInfo[1]
            multipoleElement = etree.SubElement(
                amoebaMultipoleForceElement, "Multipole"
            )
            multipoleElement.attrib["type"] = axisInfo[0]

            axisInfoLen = len(axisInfo)

            if axisInfoLen > 1:
                multipoleElement.attrib["kz"] = axisInfo[1]
            if axisInfoLen > 2:
                multipoleElement.attrib["kx"] = axisInfo[2]
            if axisInfoLen > 3:
                multipoleElement.attrib["ky"] = axisInfo[3]

            multipoleElement.attrib["c0"] = multipoles[0]
            multipoleElement.attrib["d1"] = str(dipoleConversion * float(multipoles[1]))
            multipoleElement.attrib["d2"] = str(dipoleConversion * float(multipoles[2]))
            multipoleElement.attrib["d3"] = str(dipoleConversion * float(multipoles[3]))
            multipoleElement.attrib["q11"] = str(
                quadrupoleConversion * float(multipoles[4])
            )
            multipoleElement.attrib["q21"] = str(
                quadrupoleConversion * float(multipoles[5])
            )
            multipoleElement.attrib["q22"] = str(
                quadrupoleConversion * float(multipoles[6])
            )
            multipoleElement.attrib["q31"] = str(
                quadrupoleConversion * float(multipoles[7])
            )
            multipoleElement.attrib["q32"] = str(
                quadrupoleConversion * float(multipoles[8])
            )
            multipoleElement.attrib["q33"] = str(
                quadrupoleConversion * float(multipoles[9])
            )

        polarityConversion = 0.001
        for polarize in forces["polarize"]:
            polarizeElement = etree.SubElement(amoebaMultipoleForceElement, "Polarize")
            polarizeElement.attrib["type"] = polarize[0]
            polarizeElement.attrib["polarizability"] = str(
                polarityConversion * float(polarize[1])
            )
            polarizeElement.attrib["thole"] = polarize[2]
            for ii in range(3, len(polarize)):
                polarizeElement.attrib[f"pgrp{ii - 2}"] = polarize[ii]

    @staticmethod
    def _writeAmoebaGeneralizedKirkwoodForce(
        root: etree.Element,
        forces: Dict[str, Any],
        atomTypes: Dict[str, Dict[str, Any]],
    ) -> None:
        """
        Write the AmoebaGeneralizedKirkwoodForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        atomTypes : dict
            The atom types dictionary.
        bioTypes : dict
            The bio types dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "multipole" not in forces:
            return

        solventDielectric = 78.3
        soluteDielectric = 1.0
        includeCavityTerm = 1
        probeRadius = 0.14
        surfaceAreaFactor = -6.0 * 3.1415926535 * 0.0216 * 1000.0 * 0.4184

        gkForceElement = etree.SubElement(root, "AmoebaGeneralizedKirkwoodForce")
        gkForceElement.attrib["solventDielectric"] = str(solventDielectric)
        gkForceElement.attrib["soluteDielectric"] = str(soluteDielectric)
        gkForceElement.attrib["includeCavityTerm"] = str(includeCavityTerm)
        gkForceElement.attrib["probeRadius"] = str(probeRadius)
        gkForceElement.attrib["surfaceAreaFactor"] = str(surfaceAreaFactor)

        for multipoleInfo in forces["multipole"]:
            axisInfo = multipoleInfo[0]
            multipoles = multipoleInfo[1]
            atomType = axisInfo[0]

            shct = 0.8  # Default value for overlap scale factor
            if atomType in atomTypes:
                element = atomTypes[atomType]["element"]
                if element == "H":
                    shct = 0.85
                elif element == "C":
                    shct = 0.72
                elif element == "N":
                    shct = 0.79
                elif element == "O":
                    shct = 0.85
                elif element == "F":
                    shct = 0.88
                elif element == "P":
                    shct = 0.86
                elif element == "S":
                    shct = 0.96
                elif element == "Fe":
                    shct = 0.88
                else:
                    print(
                        f"Warning no overlap scale factor for type={atomType} element={element}"
                    )
            else:
                print(f"Warning no overlap scale factor for type={atomType} ")

            gkElement = etree.SubElement(gkForceElement, "GeneralizedKirkwood")
            gkElement.attrib["type"] = axisInfo[0]
            gkElement.attrib["charge"] = multipoles[0]
            gkElement.attrib["shct"] = str(shct)

    @staticmethod
    def _writeAmoebaWcaDispersionForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaWcaDispersionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "vdw" not in forces:
            return

        epso = 0.1100
        epsh = 0.0135
        rmino = 1.7025
        rminh = 1.3275
        awater = 0.033428
        slevy = 1.0
        dispoff = 0.26
        shctd = 0.81

        root = etree.SubElement(root, "AmoebaWcaDispersionForce")
        root.attrib["epso"] = str(epso * 4.184)
        root.attrib["epsh"] = str(epsh * 4.184)
        root.attrib["rmino"] = str(rmino * 0.1)
        root.attrib["rminh"] = str(rminh * 0.1)
        root.attrib["awater"] = str(1000.0 * awater)
        root.attrib["slevy"] = str(slevy)
        root.attrib["dispoff"] = str(0.1 * dispoff)
        root.attrib["shctd"] = str(shctd)
        convert = 0.1

        if scalars["radiustype"] == "SIGMA":
            convert *= 1.122462048309372

        if scalars["radiussize"] == "DIAMETER":
            convert *= 0.5

        for vdw in forces["vdw"]:
            sigma = float(vdw[1])
            sigma *= convert
            epsilon = float(vdw[2]) * 4.184
            vdwElement = etree.SubElement(root, "WcaDispersion")
            vdwElement.attrib["class"] = vdw[0]
            vdwElement.attrib["radius"] = str(sigma)
            vdwElement.attrib["epsilon"] = str(epsilon)

    @staticmethod
    def _writeAmoebaUreyBradleyForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaUreyBradleyForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "ureybrad" not in forces:
            return

        cubic = 0.0
        quartic = 0.0
        ureyBradleyForceElement = etree.SubElement(root, "AmoebaUreyBradleyForce")
        ureyBradleyForceElement.attrib["cubic"] = str(cubic)
        ureyBradleyForceElement.attrib["quartic"] = str(quartic)

        for ub in forces["ureybrad"]:
            k = float(ub[3]) * 4.184 * 100.0
            d = float(ub[4]) * 0.1
            ureyBradleyElement = etree.SubElement(
                ureyBradleyForceElement, "UreyBradley"
            )
            ureyBradleyElement.attrib["class1"] = ub[0]
            ureyBradleyElement.attrib["class2"] = ub[1]
            ureyBradleyElement.attrib["class3"] = ub[2]
            ureyBradleyElement.attrib["k"] = str(k)
            ureyBradleyElement.attrib["d"] = str(d)

    def _createXmlFile(self, filename, atomTypes, bioTypes, forces, scalars, atomData):
        """
        Create the XML file.

        Parameters
        ----------
        filename : str
            The name of .key file.
        atomTypes : dict
            The atom types dictionary.
        bioTypes : dict
            The bio types dictionary.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        atomData : dict
            The atom data dictionary.

        Returns
        -------
        tinkXmlFile, gkXmlFile : tuple(file, file)
            The XML files.
        """
        root = etree.Element("ForceField")

        # Write the header
        info = etree.SubElement(root, "Info")
        etree.SubElement(info, "Source").text = os.path.basename(filename)
        etree.SubElement(info, "DateGenerated").text = datetime.date.today().isoformat()
        etree.SubElement(info, "Reference")

        # Copy the root element for the Generalized Kirkwood force
        gkRoot = root.__copy__()

        # Write the AtomTypes
        TinkerFiles._writeAtomTypes(root, atomTypes)

        # Extract residues present in this file
        residuesSet = {
            atomInfo["residue"]
            for atomInfo in atomData.values()
            if atomInfo["atomType"] in atomTypes
        }

        # Write the Residues
        TinkerFiles._writeResidues(root, atomData, self.topology, residuesSet)

        # Write the Forces
        # Bonded forces
        TinkerFiles._writeAmoebaBondForce(root, forces, scalars)
        TinkerFiles._writeAmoebaAngleForce(root, forces, scalars)
        TinkerFiles._writeAmoebaOutOfPlaneBendForce(root, forces, scalars)
        TinkerFiles._writeAmoebaTorsionForce(root, forces, scalars)
        TinkerFiles._writeAmoebaPiTorsionForce(root, forces)
        TinkerFiles._writeAmoebaStretchTorsionForce(root, forces)
        TinkerFiles._writeAmoebaAngleTorsionForce(root, forces)
        TinkerFiles._writeAmoebaStretchBendForce(root, forces)
        TinkerFiles._writeAmoebaTorsionTorsionForce(root, forces)

        # Non-bonded forces
        TinkerFiles._writeAmoebaVdwForce(root, forces, scalars)
        TinkerFiles._writeAmoebaMultipoleForce(root, forces, scalars)
        TinkerFiles._writeAmoebaUreyBradleyForce(root, forces)

        # Write the Generalized Kirkwood Force
        TinkerFiles._writeAmoebaGeneralizedKirkwoodForce(gkRoot, forces, atomTypes)
        TinkerFiles._writeAmoebaWcaDispersionForce(gkRoot, forces, scalars)

        # Write footer
        tree = etree.ElementTree(root)
        etree.indent(root, "  ")
        gkTree = etree.ElementTree(gkRoot)
        etree.indent(gkRoot, "  ")

        # Store the XML file in memory
        tinkerXmlFile = io.StringIO()
        tree.write(
            tinkerXmlFile, xml_declaration=True, encoding="unicode", method="xml"
        )

        gkTinkerXmlFile = io.StringIO()
        gkTree.write(
            gkTinkerXmlFile, xml_declaration=True, encoding="unicode", method="xml"
        )

        if self._writeXmlFiles:
            if scalars["forcefield"] == "-2.55":
                scalars["forcefield"] = filename.rsplit(".", 1)[0]

            tinkerXmlFileName = scalars["forcefield"] + ".xml"
            with open(tinkerXmlFileName, "w") as f:
                f.write(tinkerXmlFile.getvalue())

            gkTinkerXmlFileName = scalars["forcefield"] + "_gk.xml"
            with open(gkTinkerXmlFileName, "w") as f:
                f.write(gkTinkerXmlFile.getvalue())

        return tinkerXmlFile, gkTinkerXmlFile
