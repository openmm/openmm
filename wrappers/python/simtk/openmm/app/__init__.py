"""
OpenMM Application
"""
from __future__ import absolute_import

__docformat__ = "epytext en"

__author__ = "Peter Eastman"
__copyright__ = "Copyright 2016-2019, Stanford University and Peter Eastman"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Peter Eastman"
__email__ = "peastman@stanford.edu"

from .topology import Topology, Chain, Residue, Atom
from .pdbfile import PDBFile
from .pdbxfile import PDBxFile
from .forcefield import ForceField
from .simulation import Simulation
from .pdbreporter import PDBReporter, PDBxReporter
from .amberprmtopfile import AmberPrmtopFile, HCT, OBC1, OBC2, GBn, GBn2
from .amberinpcrdfile import AmberInpcrdFile
from .dcdfile import DCDFile
from .gromacsgrofile import GromacsGroFile
from .gromacstopfile import GromacsTopFile
from .dcdreporter import DCDReporter
from .modeller import Modeller
from .statedatareporter import StateDataReporter
from .element import Element
from .desmonddmsfile import DesmondDMSFile
from .checkpointreporter import CheckpointReporter
from .charmmcrdfiles import CharmmCrdFile, CharmmRstFile
from .charmmparameterset import CharmmParameterSet
from .charmmpsffile import CharmmPsfFile, CharmmPSFWarning
from .simulatedtempering import SimulatedTempering
from .metadynamics import Metadynamics, BiasVariable

# Enumerated values

NoCutoff = forcefield.NoCutoff
CutoffNonPeriodic = forcefield.CutoffNonPeriodic
CutoffPeriodic = forcefield.CutoffPeriodic
Ewald = forcefield.Ewald
PME = forcefield.PME
LJPME = forcefield.LJPME

HBonds = forcefield.HBonds
AllBonds = forcefield.AllBonds
HAngles = forcefield.HAngles

Single = topology.Single
Double = topology.Double
Triple = topology.Triple
Aromatic = topology.Aromatic
Amide = topology.Amide

