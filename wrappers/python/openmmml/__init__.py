from .mlpotential import MLPotential
from . import models

# Register CACE potential
from openmmml.models.cacepotential import CACEPotentialImplFactory
MLPotential.registerImplFactory('cace-lr', CACEPotentialImplFactory())

# Register UMA batched RPMD potentials
from openmmml.models.umapotential_pythonforce_batch import UMAPotentialPythonForceBatchedImplFactory

# UMA models
MLPotential.registerImplFactory('uma-s-1-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
MLPotential.registerImplFactory('uma-s-1p1-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
MLPotential.registerImplFactory('uma-m-1p1-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())

# OMol25 eSEN models
MLPotential.registerImplFactory('esen-md-direct-all-omol-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
MLPotential.registerImplFactory('esen-sm-conserving-all-omol-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
MLPotential.registerImplFactory('esen-sm-direct-all-omol-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())

# OC25 eSEN models
MLPotential.registerImplFactory('esen-sm-conserving-all-oc25-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
MLPotential.registerImplFactory('esen-md-direct-all-oc25-pythonforce-batch', UMAPotentialPythonForceBatchedImplFactory())
